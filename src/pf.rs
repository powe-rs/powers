use crate::mpc::*;
use crate::mpopt::MPOpt;
use crate::total_load::{total_load, LoadType, LoadZone};
use ndarray::Array1;
use num_complex::Complex64;
use sprs::CsMatView;
use std::collections::HashMap;
use std::f64::consts::PI;
use std::ops::Deref;

pub(crate) fn pfsoln(
    base_mva: f64,
    bus0: &[Bus],
    gen0: &[Gen],
    branch0: &[Branch],
    Ybus: CsMatView<Complex64>,
    Yf: CsMatView<Complex64>,
    Yt: CsMatView<Complex64>,
    V: Array1<Complex64>,
    refbus: &[usize],
    pv: &[usize],
    pq: &[usize],
    mpopt: &MPOpt,
) -> (Vec<Bus>, Vec<Gen>, Vec<Branch>) {
    // let (mut bus, mut gen, mut branch) = (bus0.clone(), gen0.clone(), branch0.clone());
    let (mut bus, mut gen, mut branch) = (Vec::from(bus0), Vec::from(gen0), Vec::from(branch0));
    // let mut bus: Vec<Bus> = Vec::from(bus0);

    for (i, b) in bus.iter_mut().enumerate() {
        b.vm = V[i].norm();
        b.va = V[i].arg() * 180.0 / PI;
    }

    let i_bus = &Ybus * &V;
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

    let is_on = |g: &Gen| g.status && bus[g.bus].bus_type != BusType::PQ;

    let nb = bus.len();
    let mut ngb = vec![0; nb];
    for (i, g) in gen.iter_mut().enumerate() {
        let s_bus = V[g.bus] * i_bus[g.bus].conj(); // compute total injected bus power
        if is_on(&g.deref()) {
            g.qg = s_bus.im * base_mva + qd_gbus[g.bus]; // inj Q + local Qd

            ngb[g.bus] += 1;
        } else if !g.status {
            g.qg = 0.0;
        }
    }
    let ngg = gen
        .iter()
        .filter(|&g| is_on(g))
        .map(|g| ngb[g.bus])
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
    let mut cg = HashMap::<usize, Vec<usize>>::new();
    for (i, g) in gen.iter().enumerate() {
        if is_on(g) {
            // cg[g.bus] = append(cg[g.bus], gen)
            if cg.contains_key(&g.bus) {
                cg.get_mut(&g.bus).unwrap().push(i)
            } else {
                cg.insert(g.bus, vec![i]);
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
            if g.bus == r {
                let s_bus = V[g.bus] * i_bus[g.bus].conj();

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
    let i_fr_bus = &Yf * &V;
    let i_to_bus = &Yf * &V;
    for (i, br) in branch.iter_mut().enumerate() {
        if br.status {
            let s_f: Complex64 = V[i] * i_fr_bus[i].conj() * base_mva; // complex power at "from" bus
            let s_t: Complex64 = V[i] * i_to_bus[i].conj() * base_mva; // complex power injected at "to" bus

            br.pf = s_f.re;
            br.qf = s_f.im;
            br.pt = s_t.re;
            br.qt = s_t.im;
        } else {
            br.pf = 0.0;
            br.qf = 0.0;
            br.pt = 0.0;
            br.qt = 0.0;
        }
    }

    return (Vec::from(bus), Vec::from(gen), Vec::from(branch));
}
