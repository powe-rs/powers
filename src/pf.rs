use crate::dc::{dc_pf, make_b_dc};
use crate::fd;
use crate::gauss;
use crate::mpc::MPC;
use crate::mpopt::{Alg, BusVoltage, GenQLimits, MPOpt, NodalBalance};
use crate::newton::*;
use crate::order::{ext_to_int, int_to_ext};
use crate::powers::{bus_types, make_d_sbus_d_vm, make_sbus, make_ybus, SBus};
use crate::radial::radial_pf;
use crate::total_load::{total_load, LoadType, LoadZone};

use crate::debug::format_polar_vec;
use anyhow::Result;
use casecsv::*;
use num_complex::Complex64;
use sparsetools::csr::CSR;
use spsolve::Solver;
use std::collections::{HashMap, HashSet};
use std::f64::consts::PI;
use std::ops::Deref;
use std::time::Instant;

struct MakeSBus<'a> {
    base_mva: f64,
    bus: &'a [Bus],
    gen: &'a [Gen],
    mpopt: &'a MPOpt,
}

impl<'a> SBus for MakeSBus<'a> {
    fn s_bus(&self, v_m: &[f64]) -> Vec<Complex64> {
        make_sbus(
            self.base_mva,
            self.bus,
            self.gen,
            self.mpopt,
            Some(v_m),
            None,
        )
    }

    fn d_sbus_d_vm(&self, v_m: &[f64]) -> CSR<usize, Complex64> {
        make_d_sbus_d_vm(
            self.base_mva,
            self.bus,
            self.gen,
            self.mpopt,
            Some(v_m),
            None,
        )
    }
}

pub fn runpf(
    casedata: &MPC,
    mpopt: &MPOpt,
    solver: &dyn Solver<usize, f64>,
) -> Result<(MPC, bool)> {
    // options
    let qlim = mpopt.pf.enforce_q_limits != GenQLimits::IgnoreLimits; // enforce Q limits on gens?
    let dc = mpopt.dc; // use DC formulation?

    // read data
    let mpc = casedata;

    // convert to internal indexing
    let mut mpc = ext_to_int(mpc);
    let (base_mva, mut bus, mut gen, mut branch) = (mpc.base_mva, mpc.bus, mpc.gen, mpc.branch);

    let (_t0, success, its) = if !bus.is_empty() {
        // get bus index lists of each type of bus
        let (ref_, pv, pq) = bus_types(&bus, &gen);

        //-----  run the power flow  -----
        let t0 = Instant::now();
        let mut success = false;
        let mut its = 0; // total iterations

        // if mpopt.verbose > 0 {
        //     v = mpver('all');
        //     fprintf('\nPowers Version % s, %s', v.Version, v.Date);
        // }
        if dc {
            // DC formulation
            // if mpopt.verbose > 0 {
            log::info!("DC Power Flow");
            // }
            // initial state
            let v_a0: Vec<f64> = bus.iter().map(|b| b.va * PI / 180.0).collect();

            // build B matrices and phase shift injections
            let (b_dc, b_f, p_businj, p_finj) = make_b_dc(base_mva, &bus, &branch);

            // compute complex bus power injections (generation - load)
            // adjusted for phase shifters and real shunts
            let s_bus = make_sbus(base_mva, &bus, &gen, &mpopt, None, None);
            // let gs = bus.iter().map(|b| b.gs).collect_vec();
            // let Pbus: Vec<f64> = s_bus.real() - Pbusinj - gs / baseMVA;
            let p_bus: Vec<f64> = (0..bus.len())
                .map(|i| s_bus[i].re - p_businj[i] - bus[i].gs / base_mva)
                .collect();

            // "run" the power flow
            let (v_a, succ): (Vec<f64>, bool) =
                dc_pf(&b_dc, &p_bus, &v_a0, &ref_, &pv, &pq, solver)?;
            success = succ;

            its = 1;

            // update data matrices with solution
            let pf: Vec<f64> = (b_f * &v_a)
                .iter()
                .enumerate()
                .map(|(i, pf)| (pf + p_finj[i]) * base_mva)
                .collect();
            for (i, br) in branch.iter_mut().enumerate() {
                br.qf = Some(0.0);
                br.qt = Some(0.0);
                br.pf = Some(pf[i]);
                br.pt = Some(-pf[i]);
            }
            for (i, b) in bus.iter_mut().enumerate() {
                b.vm = 1.0;
                b.va = v_a[i] * 180.0 / PI;
            }
            // update Pg for slack generator (1st gen at ref bus)
            // (note: other gens at ref bus are accounted for in Pbus)
            //      Pg = Pinj + Pload + Gs
            //      newPg = oldPg + newPinj - oldPinj
            let b_ref = b_dc.select(Some(&ref_), None)?;
            let p_ref = b_ref * &v_a;
            for r in ref_ {
                for g in gen.iter_mut() {
                    if g.gen_bus == r {
                        g.pg += (p_ref[r] - p_bus[r]) * base_mva;
                        break;
                    }
                }
            }
        } else {
            let alg = mpopt.pf.algorithm;

            // initial state
            // let v0 = Arr::ones(bus.len()); // flat start
            let mut v0 = Vec::with_capacity(bus.len());
            for b in bus.iter() {
                v0.push(Complex64::from_polar(b.vm, b.va * PI / 180.0));
            }
            let pq_bus: HashSet<usize> = HashSet::from_iter(pq.clone().into_iter()); // exclude PQ buses
            for g in gen.iter() {
                if g.is_on() && !pq_bus.contains(&g.gen_bus) {
                    v0[g.gen_bus] =
                        Complex64::new(g.vg / (v0[g.gen_bus] * v0[g.gen_bus]).norm(), 0.0);
                }
            }
            log::debug!("V0: {}", format_polar_vec(&v0));

            if qlim {
                // let ref0 = ref_.clone(); // save index and angle of
                // let Varef0 = ref0.iter().map(|&r0| bus[r0].va).collect::<Vec<f64>>(); //   original reference bus(es)
                // let mut limited = vec![]; // list of indices of gens @ Q lims
                // let mut fixedQg = Arr::zeros(gen.len()); // Qg of gens at Q limits
            }

            // build admittance matrices
            let (y_bus, y_f, y_t) = make_ybus(base_mva, &bus, &branch, true);
            let (y_bus, y_f, y_t) = (y_bus.to_csr(), y_f.unwrap().to_csr(), y_t.unwrap().to_csr());
            log::trace!("Ybus:\n{}", y_bus.to_table());

            let mut repeat = true;
            while repeat {
                // function for computing V dependent complex bus power injections
                // (generation - load)
                // let s_bus: SBus = |v_m: &Arr<f64>| {
                //     make_sbus(baseMVA, &bus, &gen, mpopt, Some(v_m), None, false).0;
                // };
                let s_bus = MakeSBus {
                    base_mva,
                    bus: &bus,
                    gen: &gen,
                    mpopt: &mpopt,
                };

                let (v, succ, iterations) = match alg {
                    Alg::NR => {
                        let newtonpf_fcn = match mpopt.pf.current_balance {
                            NodalBalance::CURRENT => {
                                match mpopt.pf.v_cartesian {
                                    BusVoltage::POLAR => {
                                        newtonpf_i_polar // current, polar
                                    }
                                    BusVoltage::CARTESIAN => {
                                        newtonpf_i_cart // current, cartesian
                                    }
                                    BusVoltage::HYBRID => {
                                        newtonpf_i_hybrid // current, hybrid
                                    }
                                }
                            }
                            NodalBalance::POWER => {
                                match mpopt.pf.v_cartesian {
                                    BusVoltage::POLAR => {
                                        newtonpf_s_polar // default - power, polar
                                    }
                                    BusVoltage::CARTESIAN => {
                                        newtonpf_s_cart // power, cartesian
                                    }
                                    BusVoltage::HYBRID => {
                                        newtonpf_s_hybrid // power, hybrid
                                    }
                                }
                            }
                        };
                        newtonpf_fcn(&y_bus, &s_bus, &v0, &ref_, &pv, &pq, solver, &mpopt, None)?
                    }
                    Alg::FDBX | Alg::FDXB => {
                        let (b_p, b_pp) = fd::make_b(base_mva, &bus, &branch, alg, true);
                        let (b_p, b_pp) = (b_p.to_csr(), b_pp.unwrap().to_csr());
                        let progress = fd::PrintProgress {};
                        fd::fdpf(
                            &y_bus,
                            &s_bus,
                            &v0,
                            &b_p,
                            &b_pp,
                            &ref_,
                            &pv,
                            &pq,
                            solver,
                            &mpopt,
                            Some(&progress),
                        )?
                    }
                    Alg::GS => {
                        let progress = gauss::PrintProgress {};
                        gauss::gausspf(
                            &y_bus,
                            &s_bus,
                            &v0,
                            &ref_,
                            &pv,
                            &pq,
                            solver,
                            &mpopt,
                            Some(&progress),
                        )?
                    }
                    Alg::SUM => {
                        let (_mpc, success, iterations) =
                            radial_pf(base_mva, &bus, &gen, &branch, mpopt)?;
                        (Vec::new(), success, iterations)
                    }
                };
                success = succ;
                its = its + iterations;

                // update data matrices with solution
                match alg {
                    Alg::NR | Alg::FDBX | Alg::FDXB | Alg::GS => {
                        let mpc_s = pfsoln(
                            base_mva, &bus, &gen, &branch, &y_bus, &y_f, &y_t, &v, &ref_, &pv, &pq,
                            &mpopt,
                        );
                        bus = mpc_s.0;
                        gen = mpc_s.1;
                        branch = mpc_s.2;
                    }
                    _ => {}
                }

                if success && qlim {
                    // enforce generator Q limits
                    unimplemented!("generator Q limits");
                } else {
                    repeat = false; // don't enforce generator Q limits, once is enough
                }
            }
            // TODO: adjust voltage angles to make original ref bus correct
        }
        (t0, success, its)
    } else {
        // if mpopt.verbose {
        log::error!("Power flow not valid: Case contains no connected buses");
        // }
        (Instant::now(), false, 0)
    };
    // mpc.et = t0;
    mpc.success = Some(success);
    mpc.iterations = Some(its);

    // -----  output results  ----- //
    // convert back to original bus numbering & print results
    mpc.bus = bus;
    mpc.gen = gen;
    mpc.branch = branch;
    let mut results = int_to_ext(&mpc).unwrap();

    let order = results.order.take().unwrap();

    // zero out result fields of out-of-service gens & branches
    let off = order.gen.status.off.clone();
    for i in off {
        results.gen[i].pg = 0.0;
        results.gen[i].qg = 0.0;
    }
    let off = order.branch.status.off.clone();
    for i in off {
        results.branch[i].pf = Some(0.0);
        results.branch[i].qf = Some(0.0);
        results.branch[i].pt = Some(0.0);
        results.branch[i].qt = Some(0.0);
    }

    // printpf(&results, 1, mpopt);

    Ok((results, success))
}

fn _have_zip_loads(mpopt: &MPOpt) -> bool {
    if let Some(pw) = mpopt.exp.sys_wide_zip_loads.pw {
        // TODO: check indexing
        if pw[1..2].iter().any(|&v| v != 0.0) {
            return true;
        }
    }
    if let Some(qw) = mpopt.exp.sys_wide_zip_loads.qw {
        if qw[1..2].iter().any(|&v| v != 0.0) {
            return true;
        }
    }
    false
    // (mpopt.exp.sys_wide_zip_loads.pw.is_some()
    //     && any(&mpopt.exp.sys_wide_zip_loads.pw.unwrap()[1..2]))
    //     || (mpopt.exp.sys_wide_zip_loads.qw.is_some()
    //         && any(&mpopt.exp.sys_wide_zip_loads.qw.unwrap()[1..2])) // TODO: check indexing
}

pub(crate) fn pfsoln(
    base_mva: f64,
    bus0: &[Bus],
    gen0: &[Gen],
    branch0: &[Branch],
    y_bus: &CSR<usize, Complex64>,
    y_f: &CSR<usize, Complex64>,
    y_t: &CSR<usize, Complex64>,
    v: &[Complex64],
    refbus: &[usize],
    _pv: &[usize],
    _pq: &[usize],
    mpopt: &MPOpt,
) -> (Vec<Bus>, Vec<Gen>, Vec<Branch>) {
    // let (mut bus, mut gen, mut branch) = (bus0.clone(), gen0.clone(), branch0.clone());
    let (mut bus, mut gen, mut branch) = (Vec::from(bus0), Vec::from(gen0), Vec::from(branch0));
    // let mut bus: Vec<Bus> = Vec::from(bus0);

    for (i, b) in bus.iter_mut().enumerate() {
        b.vm = v[i].norm();
        b.va = v[i].arg() * 180.0 / PI;
    }

    let i_bus = y_bus * v;
    // let s_bus = gen.iter().map(|g| V[g.bus] * i_bus[g.bus].conj()).collect();

    let (_pd_gbus, qd_gbus) = total_load(
        &bus, /*(gbus, :)*/
        None,
        LoadZone::Bus,
        None,
        LoadType::Fixed,
        false,
        mpopt,
        true,
    );
    let qd_gbus = qd_gbus.as_ref().unwrap();

    let is_on = |g: &Gen| g.is_on() && bus[g.gen_bus].is_pq();

    let nb = bus.len();
    let mut ngb = vec![0; nb];
    for g in gen.iter_mut() {
        let s_bus = v[g.gen_bus] * i_bus[g.gen_bus].conj(); // compute total injected bus power
        if is_on(&g.deref()) {
            g.qg = s_bus.im * base_mva + qd_gbus[g.gen_bus]; // inj Q + local Qd

            ngb[g.gen_bus] += 1;
        } else if g.is_off() {
            g.qg = 0.0;
        }
    }
    let ngg = gen
        .iter()
        .filter(|&g| is_on(g))
        .map(|g| ngb[g.gen_bus])
        .collect::<Vec<usize>>(); // number of gens at this gen's bus

    // ...at this point any buses with more than one generator will have
    // the total Q dispatch for the bus assigned to each generator. This
    // must be split between them. We do it first equally, then in proportion
    // to the reactive range of the generator.
    for (i, g) in gen.iter_mut().filter(|g| is_on(&g.deref())).enumerate() {
        g.qg = g.qg / ngg[i] as f64;

        let mut m = g.qg.abs();
        if !g.qmax.is_infinite() {
            m = m + g.qmax.abs();
        }
        if !g.qmin.is_infinite() {
            m = m + g.qmin.abs();
        }
        // TODO: each gen gets sum over all gens at same bus

        // replace +/- Inf limits with proxy +/- M
        if g.qmin.is_infinite() {
            g.qmin = if g.qmin.is_sign_positive() { m } else { -m };
        }
        if g.qmin.is_infinite() {
            g.qmax = if g.qmax.is_sign_positive() { m } else { -m };
        }
    }

    // Buses with more than one generator have the total reactive power
    // dispatch for the bus divided equally among each online generator.
    // The reactive power is divided in proportion to the reactive range
    // of each generator, according to the logic in the pfsoln function
    // by Ray Zimmerman from MATPOWER v7.
    let mut cg: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, g) in gen.iter().enumerate() {
        if is_on(g) {
            // cg[g.bus] = append(cg[g.bus], gen)
            if cg.contains_key(&g.gen_bus) {
                cg.get_mut(&g.gen_bus).unwrap().push(i)
            } else {
                cg.insert(g.gen_bus, vec![i]);
            }
            // match cg.get_mut(&g.bus) {
            //     Some(l) => l.push(g),
            //     None => cg.insert(g.bus, vec![g]),
            // }
        }
    }

    for l in cg.values() {
        if l.len() < 2 {
            continue;
        }
        let mut qg_tot = 0.0; // Total Qg at the bus.
        for &i in l.iter() {
            let g: &mut Gen = gen.get_mut(i).unwrap();
            qg_tot += g.qg;
        }

        // The sum of absolute Qg, Qmax and Qmin for all generators
        // at the bus (if Qmax/Qmin is non-infinite). Used as limit
        // when Qmax/Qmin is infinite.
        let mut m = 0.0;
        for &i in l.iter() {
            let g: &Gen = gen.get(i).unwrap();

            let mut mg = g.qg.abs();
            if !g.qmax.is_infinite() {
                mg = mg + g.qmax.abs();
            }
            if !g.qmin.is_infinite() {
                mg = mg + g.qmin.abs();
            }

            m += mg
        }
        let mut qmin = vec![0.0; l.len()];
        let mut qmax = vec![0.0; l.len()];

        for (i, &j) in l.iter().enumerate() {
            let g: &Gen = gen.get(j).unwrap();

            qmin[i] = g.qmin;
            qmax[i] = g.qmax;

            // replace +/- Inf limits with proxy +/- M
            if g.qmin.is_infinite() {
                qmin[i] = if g.qmin.is_sign_positive() { m } else { -m };
            }
            if g.qmax.is_infinite() {
                qmax[i] = if g.qmax.is_sign_positive() { m } else { -m };
            }
        }
        let qg_min: f64 = qmin.iter().sum(); // Minimum total Qg at the bus.
        let qg_max: f64 = qmax.iter().sum(); // Maximum total Qg at the bus.

        if (qg_min - qg_max).abs() > 1e-13 {
            let q = (qg_tot - qg_min) / (qg_max - qg_min);
            // for (i, mut pv) in l.iter_mut().enumerate() {
            for &i in l {
                let g: &mut Gen = gen.get_mut(i).unwrap();

                //pv.Qg = pv.QMin + (((QgTot - QgMin) / (QgMax - QgMin)) / (pv.QMax - pv.QMin))
                g.qg = g.qmin + (q * (g.qmax - g.qmin));
            }
        } else {
            // Zero Qg range at bus. Qg set such that all generators
            // at the bus violate their limits by the same amount.

            // Total mismatch at bus divided by number of online generators.
            let mis = (qg_tot - qg_min) / (l.len() as f64);
            for (i, &j) in l.iter().enumerate() {
                let g: &mut Gen = gen.get_mut(j).unwrap();

                g.qg = qmin[i] + mis;
            }
        }
    }

    // Update Pg for slack gen(s).
    for &r in refbus {
        let (pd_refk, _) = total_load(
            &vec![bus[r].clone()],
            None,
            LoadZone::Bus,
            None,
            LoadType::Fixed,
            false,
            mpopt,
            false,
        );
        // let mut refgen = 0;
        let mut refgen0: Option<usize> = None;
        let mut sum = 0.0;
        for (i, g) in gen.iter_mut().enumerate() {
            if g.gen_bus == r {
                let s_bus = v[g.gen_bus] * i_bus[g.gen_bus].conj();

                g.pg = s_bus.re * base_mva + pd_refk[0]; // inj P + local Pd

                if refgen0.is_some() {
                    // more than one generator at this ref bus
                    sum += g.pg;
                }
                refgen0 = Some(i);
            }
        }
        // subtract off what is generated by other gens at this bus
        if let Some(i) = refgen0 {
            gen[i].pg -= sum;
        }
    }

    // Update/compute branch power flows.
    let i_fr_bus = y_f * v;
    let i_to_bus = y_t * v;
    for (i, br) in branch.iter_mut().enumerate() {
        if br.is_on() {
            let s_f: Complex64 = v[i] * i_fr_bus[i].conj() * base_mva; // complex power at "from" bus
            let s_t: Complex64 = v[i] * i_to_bus[i].conj() * base_mva; // complex power injected at "to" bus

            br.pf = Some(s_f.re);
            br.qf = Some(s_f.im);
            br.pt = Some(s_t.re);
            br.qt = Some(s_t.im);
        } else {
            br.pf = Some(0.0);
            br.qf = Some(0.0);
            br.pt = Some(0.0);
            br.qt = Some(0.0);
        }
    }

    return (Vec::from(bus), Vec::from(gen), Vec::from(branch));
}

// fn printpf<W: Write>(_results: &MPC, _fd: W, _mpopt: &MPOpt) {}
