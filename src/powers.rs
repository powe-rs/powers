// Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
// by Ray Zimmerman, PSERC Cornell
// Copyright (c) 2022-2023, Richard Lincoln. All rights reserved.

use crate::mpopt::MPOpt;

use anyhow::Result;
use casecsv::*;
use num_complex::Complex64;
use sparsetools::coo::Coo;
use sparsetools::csr::CCSR;
use sparsetools::csr::CSR;
use std::collections::HashSet;
use std::f64::consts::PI;
use std::iter::zip;

pub(crate) const J: Complex64 = Complex64 { re: 0.0, im: 1.0 };

#[macro_export]
macro_rules! cmplx {
    () => {
        num_complex::Complex64::new(0.0, 0.0)
    };
    ($arg1:expr) => {
        num_complex::Complex64::new($arg1, 0.0)
    };
    ($arg1:expr, $arg2:expr) => {
        num_complex::Complex64::new($arg1, $arg2)
    };
}

/// Builds index lists for each type of bus (REF, PV, PQ).
///
/// Generators with "out-of-service" status are treated as PQ buses with
/// zero generation (regardless of Pg/Qg values in gen). Expects `BUS` and
/// `GEN` have been converted to use internal consecutive bus numbering.
pub(crate) fn bus_types(bus: &[Bus], gen: &[Gen]) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
    // Buses with generators that are ON.
    // let mut bus_gen_status = vec![false; bus.len()];
    // for g in gen {
    //     if g.is_on() {
    //         bus_gen_status[g.bus] = true
    //     }
    // }
    let bus_gen_status = gen
        .iter()
        .filter(|g| g.is_on())
        .map(|g| g.bus)
        .collect::<HashSet<usize>>();

    // Form index lists for slack, PV, and PQ buses.
    // ref = find(bus(:, BUS_TYPE) == REF & bus_gen_status);   /// reference bus index
    // pv  = find(bus(:, BUS_TYPE) == PV  & bus_gen_status);   /// PV bus indices
    // pq  = find(bus(:, BUS_TYPE) == PQ | ~bus_gen_status);   /// PQ bus indices
    let refbus = bus
        .iter()
        .filter(|b| b.is_ref() && bus_gen_status.contains(&b.bus_i))
        .map(|b| b.bus_i)
        .collect::<Vec<usize>>();
    let pv = bus
        .iter()
        .filter(|b| b.is_pv() && bus_gen_status.contains(&b.bus_i))
        .map(|b| b.bus_i)
        .collect::<Vec<usize>>();
    let pq = bus
        .iter()
        .filter(|b| b.is_pq() || !bus_gen_status.contains(&b.bus_i))
        .map(|b| b.bus_i)
        .collect::<Vec<usize>>();

    (refbus, pv, pq)
}

// pub type SBus = fn(v_m: &[f64]) -> (Arr<Complex64>, Option<CSR<usize, Complex64>>);

pub trait SBus {
    fn s_bus(&self, v_m: &[f64]) -> Vec<Complex64>;
    fn d_sbus_d_vm(&self, v_m: &[f64]) -> CSR<usize, Complex64>;
}

/// Builds the vector of complex bus power injections.
///
/// Returns the vector of complex bus power injections, that is, generation
/// minus load. Power is expressed in per unit. If the MPOPT and VM arguments
/// are present it evaluates any ZIP loads based on the provided voltage
/// magnitude vector. If VM is empty, it assumes nominal voltage. If SG is
/// provided, it is a complex ng x 1 vector of generator power injections in
/// p.u., and overrides the PG and QG columns in GEN, using GEN only for
/// connectivity information.
///
///   [SBUS, DSBUS_DVM] = MAKESBUS(BASEMVA, BUS, GEN, MPOPT, VM)
///
/// See also MAKEYBUS.
pub(crate) fn make_sbus(
    base_mva: f64,
    bus: &[Bus],
    gen: &[Gen],
    mpopt: &MPOpt,
    vm: Option<&[f64]>,
    sg: Option<&[Complex64]>,
) -> Vec<Complex64> {
    let nb = bus.len();
    let base_mva = Complex64::new(base_mva, 0.0);

    // Form net complex bus power injection vector
    // (power injected by generators + power injected by loads).
    let mut s_bus = vec![Complex64::default(); nb];

    if let Some(sg) = sg {
        gen.iter()
            .zip(sg)
            .filter(|(g, _)| g.is_on())
            .for_each(|(g, s_pu)| {
                s_bus[g.bus] += s_pu;
            });
    } else {
        gen.iter().filter(|g| g.is_on()).for_each(|g| {
            s_bus[g.bus] += Complex64::new(g.pg, g.qg) / base_mva;
        });
    }

    // get load parameters
    let pw = mpopt.exp.sys_wide_zip_loads.pw.unwrap_or([1.0, 0.0, 0.0]);
    let qw = mpopt.exp.sys_wide_zip_loads.qw.unwrap_or(pw);

    bus.iter()
        .filter(|b| b.pd != 0.0 || b.qd != 0.0)
        .for_each(|b| {
            // Compute per-bus loads in p.u.
            let sd_z = cmplx!(b.pd * pw[2], b.qd * qw[2]) / base_mva;
            let sd_i = cmplx!(b.pd * pw[1], b.qd * qw[1]) / base_mva;
            let sd_p = cmplx!(b.pd * pw[0], b.qd * qw[0]) / base_mva;

            let vm_i = match vm {
                Some(vm) => Complex64::new(vm[b.bus_i], 0.0),
                None => Complex64::new(1.0, 0.0),
            };

            let sd = sd_p + (sd_i * vm_i) + (sd_z * (vm_i * vm_i));

            s_bus[b.bus_i] -= sd;
        });

    s_bus
}

/// Computes the partial derivative of the bus injections with respect to
/// voltage magnitude. If `vm` is empty, it assumes no voltage dependence
/// and returns a sparse zero matrix.
pub(crate) fn make_d_sbus_d_vm(
    base_mva: f64,
    bus: &[Bus],
    _gen: &[Gen],
    mpopt: &MPOpt,
    vm: Option<&[f64]>,
    _sg: Option<&[Complex64]>,
) -> CSR<usize, Complex64> {
    let nb = bus.len();
    let base_mva = Complex64::new(base_mva, 0.0);

    // get load parameters
    let pw = mpopt.exp.sys_wide_zip_loads.pw.unwrap_or([1.0, 0.0, 0.0]);
    let qw = mpopt.exp.sys_wide_zip_loads.qw.unwrap_or(pw);

    if let Some(vm) = vm {
        const TWO: Complex64 = Complex64 { re: 2.0, im: 0.0 };
        let diag: Vec<Complex64> = bus
            .iter()
            .map(|b| {
                if b.pd != 0.0 || b.qd != 0.0 {
                    let sd_z = Complex64::new(b.pd * pw[2], b.qd * qw[2]) / base_mva;
                    let sd_i = Complex64::new(b.pd * pw[1], b.qd * qw[1]) / base_mva;

                    let vm_i = Complex64::new(vm[b.bus_i], 0.0);

                    -(sd_i + TWO * vm_i * sd_z)
                } else {
                    Complex64::default()
                }
            })
            .collect();
        CSR::with_diagonal(diag)
    } else {
        CSR::with_size(nb, nb)
    }
}

/// Builds vectors of nominal complex bus power demands for ZIP loads.
///
/// Returns a struct with three fields, each an nb x 1 vectors. The fields
/// 'z', 'i' and 'p' correspond to the nominal p.u. complex power
/// (at 1 p.u. voltage magnitude) of the constant impedance, constant current,
/// and constant power portions, respectively of the ZIP load model.
pub(crate) fn make_sdzip(
    base_mva: f64,
    bus: &[Bus],
    mpopt: &MPOpt,
) -> (Vec<Complex64>, Vec<Complex64>, Vec<Complex64>) {
    let pw = mpopt.exp.sys_wide_zip_loads.pw.unwrap_or([1.0, 0.0, 0.0]);
    let qw = mpopt.exp.sys_wide_zip_loads.qw.unwrap_or(pw);

    let base_mva = cmplx!(base_mva);

    let mut sd_z = Vec::with_capacity(bus.len());
    let mut sd_i = Vec::with_capacity(bus.len());
    let mut sd_p = Vec::with_capacity(bus.len());
    for b in bus {
        sd_z.push(cmplx!(b.pd * pw[2], b.qd * qw[2]) / base_mva);
        sd_i.push(cmplx!(b.pd * pw[1], b.qd * qw[1]) / base_mva);
        sd_p.push(cmplx!(b.pd * pw[0], b.qd * qw[0]) / base_mva);
    }
    (sd_z, sd_i, sd_p)
}

/// Builds the bus admittance matrix and branch admittance matrices.
pub(crate) fn make_ybus(
    base_mva: f64,
    bus: &[Bus],
    branch: &[Branch],
    yf_yt: bool,
) -> (
    Coo<usize, Complex64>,
    Option<Coo<usize, Complex64>>,
    Option<Coo<usize, Complex64>>,
) {
    let nb = bus.len();
    let nl = branch.len();

    // For each branch, compute the elements of the branch admittance matrix where:
    //
    //      | If |   | Yff  Yft |   | Vf |
    //      |    | = |          | * |    |
    //      | It |   | Ytf  Ytt |   | Vt |
    let mut y_bus = Coo::with_size(nb, nb);
    let mut y_f = if yf_yt {
        Some(Coo::with_size(nl, nb))
    } else {
        None
    };
    let mut y_t = if yf_yt {
        Some(Coo::with_size(nl, nb))
    } else {
        None
    };

    for (i, br) in branch.iter().enumerate() {
        let y_s = if br.is_on() {
            Complex64::new(1.0, 0.0) / Complex64::new(br.r, br.x)
        } else {
            Complex64::default()
        }; // series admittance
        let b_c = if br.is_on() { br.b } else { 0.0 }; // line charging susceptance
        let t = if br.tap == 0.0 { 1.0 } else { br.tap }; // default tap ratio = 1
        let tap = Complex64::from_polar(t, br.shift * PI / 180.0); // add phase shifters

        let y_tt = y_s + Complex64::new(0.0, b_c / 2.0);
        let y_ff = y_tt / (tap * tap.conj());
        let y_ft = -y_s / tap.conj();
        let y_tf = -y_s / tap;

        let (f, t) = (br.f_bus, br.t_bus);

        if yf_yt {
            let y_f = y_f.as_mut().unwrap();
            let y_t = y_t.as_mut().unwrap();

            y_f.push(i, f, y_ff);
            y_f.push(i, t, y_ft);

            y_t.push(i, f, y_tf);
            y_t.push(i, t, y_tt);
        }

        y_bus.push(f, f, y_ff);
        y_bus.push(f, t, y_ft);
        y_bus.push(t, f, y_tf);
        y_bus.push(t, t, y_tt);
    }

    let y_sh = bus
        .iter()
        .map(|b| Complex64::new(b.gs, b.bs) / Complex64::new(base_mva, 0.0))
        .collect::<Vec<Complex64>>();

    for (i, _) in bus.iter().enumerate() {
        y_bus.push(i, i, y_sh[i]);
    }
    (y_bus, y_f, y_t)

    /*
    // series admittance
    let y_s = branch.iter().map(|br| br.y_s()).collect::<Vec<Complex64>>();
    // line charging susceptance
    let b_c = branch
        .iter()
        .map(|br| if br.status { br.b } else { 0.0 })
        .collect::<Vec<f64>>();
    let tap = branch
        .iter()
        .map(|br| {
            let t = if br.tap == 0.0 { 1.0 } else { br.tap }; // default tap ratio = 1
            Complex64::from_polar(t, br.shift * PI / 180.0) // add phase shifters
        })
        .collect::<Vec<Complex64>>();

    // Ytt = Ys + 1j*Bc/2;
    let y_tt = y_s
        .iter()
        .zip(&b_c)
        .map(|(ys, bc)| ys + cmplx!(0.0, bc / 2.0))
        .collect::<Vec<Complex64>>();
    // Yff = Ytt ./ (tap .* conj(tap));
    let y_ff = y_tt
        .iter()
        .zip(&tap)
        .map(|(ytt, t)| ytt / (t * t.conj()))
        .collect::<Vec<Complex64>>();
    // Yft = - Ys ./ conj(tap);
    let y_ft = y_s
        .iter()
        .zip(&tap)
        .map(|(ys, t)| -ys / t.conj())
        .collect::<Vec<Complex64>>();
    // Ytf = - Ys ./ tap;
    let y_tf = y_s
        .iter()
        .zip(tap)
        .map(|(ys, t)| -ys / t)
        .collect::<Vec<Complex64>>();

    // Compute shunt admittance.
    //
    // If Psh is the real power consumed by the shunt at V = 1.0 p.u.
    // and Qsh is the reactive power injected by the shunt at V = 1.0 p.u.
    // then Psh - j Qsh = V * conj(Ysh * V) = conj(Ysh) = Gs - j Bs,
    // i.e. Ysh = Psh + j Qsh, so ...
    //Ysh = (bus(:, GS) + 1j * bus(:, BS)) / baseMVA; %% vector of shunt admittances
    let y_sh = bus
        .iter()
        .map(|b| b.y_sh(base_mva))
        .collect::<Vec<Complex64>>();

    let f = branch.iter().map(|br| br.from_bus).collect::<Vec<usize>>();
    let t = branch.iter().map(|br| br.to_bus).collect::<Vec<usize>>();

    {
        let lr = Vec::from_iter(0..nl);
        let l1 = vec![cmplx!(1.0); nl];

        // Build connection matrices.
        // let c_f = TriMat::from_triplets((nl, nb), lr.clone(), f, l1.clone()).to_csr();
        let c_f = sparse((nl, nb), &lr, &f, &l1);
        // let c_t = TriMat::from_triplets((nl, nb), lr.clone(), t, l1.clone()).to_csr();
        let c_t = sparse((nl, nb), &lr, &t, &l1);

        // Build Yf and Yt such that Yf * V is the vector of complex branch currents injected
        // at each branch's "from" bus, and Yt is the same for the "to" bus end
        // let y_f = TriMat::from_triplets((nl, nl), lr.clone(), lr.clone(), y_ff).to_csr() * &c_f
        //     + TriMat::from_triplets((nl, nl), lr.clone(), lr.clone(), y_ft).to_csr() * &c_t;
        let y_f = &(&spdiag(&y_ff) * &c_f) + &(&spdiag(&y_ft) * &c_t);
        // let y_t = TriMat::from_triplets((nl, nl), lr.clone(), lr.clone(), y_tf).to_csr() * c_f
        //     + TriMat::from_triplets((nl, nl), lr.clone(), lr.clone(), y_tt).to_csr() * c_t;
        let y_t = &(&spdiag(&y_tf) * &c_f) + &(&spdiag(&y_tt) * &c_t);

        // Branch admittances.
        let y_branch = &(&c_f.transpose_view() * &y_f) + &(&c_t.transpose_view() * &y_t);
        // Shunt admittance.
        let y_shunt = spdiag(&y_sh);

        let y_bus = &y_branch + &y_shunt;

        if !yf_yt {
            (y_bus, None, None)
        } else {
            (y_bus, Some(y_f), Some(y_t))
        }
    }
    */
    /*{
        let y_bus = TriMat::from_triplets(
            (nb, nb),
            [&f, &f, &t, &t].concat(),
            [&f, &t, &f, &t].concat(),
            [&y_ff, &y_ft, &y_tf, &y_tt].concat(),
        )
        .to_csr();

        if !yf_yt {
            (y_bus, None, None)
        } else {
            let i = [0..nl, 0..nl].concat();
            let y_f =
                TriMat::from_triplets((nl, nl), i, [&f, &t].concat(), [&y_ff, &y_ft].concat()).to_csr();
            let y_t =
                TriMat::from_triplets((nl, nl), i, [&f, &t].concat(), [&y_tf, &y_tt].concat()).to_csr();

            (y_bus, Some(y_f), Some(y_t))
        }

        (y_bus, y_f, y_t)
    }*/
}

/// Computes partial derivatives of power injection w.r.t. voltage.
///
/// The derivatives can be take with respect to polar or cartesian coordinates
/// of voltage, depending on the 3rd argument.
pub(crate) fn d_sbus_d_v(
    y_bus: &CSR<usize, Complex64>,
    v: &[Complex64],
    cartesian: bool,
) -> (CSR<usize, Complex64>, CSR<usize, Complex64>) {
    let i_bus = y_bus * v;

    let diag_v = CSR::<usize, Complex64>::with_diagonal(v.to_vec());
    let diag_i_bus = CSR::<usize, Complex64>::with_diagonal(i_bus);

    if cartesian {
        // dSbus/dVr = conj(diagIbus) + diagV * conj(Ybus)
        // dSbus/dVi = 1j * (conj(diagIbus) - diagV * conj(Ybus))

        let d_sbus_d_vr = diag_i_bus.conj() + &diag_v * y_bus.conj();
        let d_sbus_d_vi = (diag_i_bus.conj() - &diag_v * y_bus.conj()) * J;

        (d_sbus_d_vr, d_sbus_d_vi)
    } else {
        let v_norm = v
            .iter()
            .map(|v| v / Complex64::new(v.norm(), 0.0))
            .collect();
        let diag_v_norm = CSR::<usize, Complex64>::with_diagonal(v_norm);

        // dSbus/dVa = 1j * diagV * conj(diagIbus - Ybus * diagV)
        // dSbus/dVm = diagV * conj(Ybus * diagVnorm) + conj(diagIbus) * diagVnorm

        let mut d_sbus_d_va = &diag_v * (&diag_i_bus - y_bus * &diag_v).conj() * J;
        let d_sbus_d_vm =
            &diag_v * (y_bus * &diag_v_norm).conj() + diag_i_bus.conj() * &diag_v_norm;

        d_sbus_d_va.sort_indexes();

        (d_sbus_d_va, d_sbus_d_vm)
    }
}

/// Computes partial derivatives of current balance w.r.t. voltage.
///
/// The derivatives can be take with respect to polar or cartesian coordinates
/// of voltage, depending on the `vcart` argument.
/// Returns two matrices containing partial derivatives of the complex bus
/// current balance w.r.t voltage magnitude and voltage angle, respectively
/// (for all buses).
pub(crate) fn d_imis_d_v(
    s_bus: &Vec<Complex64>,
    y_bus: &CSR<usize, Complex64>,
    v: &[Complex64],
    vcart: bool,
) -> Result<(CSR<usize, Complex64>, CSR<usize, Complex64>)> {
    let n = v.len();

    let (d_imis_d_v1, d_imis_d_v2) = if vcart {
        /*
        diagSV2c = sparse(1:n, 1:n, conj(Sbus./(V.^2)), n, n);
        */
        let diag_sv2c = Coo::new(
            n,
            n,
            (0..n).collect(),
            (0..n).collect(),
            // (0..n).map(|i| (s_bus[i] / (v[i] * v[i]).conj())).collect(),
            zip(s_bus, v)
                .map(|(s_bus, v)| s_bus / (v * v).conj())
                .collect(),
        )?
        .to_csr();
        let d_imis_d_v1 = y_bus + &diag_sv2c; // dImis/dVr
        let d_imis_d_v2 = J * (y_bus - diag_sv2c); // dImis/dVi

        (d_imis_d_v1, d_imis_d_v2)
    } else {
        // let v_m = Vec::from_real(&v.norm());
        // let v_norm: Vec<Complex64> = (0..n).map(|i| v[i] / cmplx!(v[i].norm())).collect();
        let v_norm = v
            .iter()
            .map(|v| v / Complex64::new(v.norm(), 0.0))
            .collect();
        // let i_bus: Vec<Complex64> = (0..n).map(|i| s_bus[i] / v[i]).collect();
        let i_bus = zip(s_bus, v).map(|(s_bus, v)| s_bus / v).collect();
        let i_bus_vm: Vec<Complex64> = zip(&i_bus, v)
            .map(|(i_bus, v)| i_bus / Complex64::new(v.norm(), 0.0))
            .collect();
        /*
        diagV       = sparse(1:n, 1:n, V, n, n);
        diagIbus    = sparse(1:n, 1:n, Ibus, n, n);
        diagIbusVm  = sparse(1:n, 1:n, Ibus./Vm, n, n);
        diagVnorm   = sparse(1:n, 1:n, V./abs(V), n, n);
        */
        let diag_v = CSR::with_diagonal(v.to_vec());
        let diag_ibus = CSR::with_diagonal(i_bus);
        let diag_ibus_vm = CSR::with_diagonal(i_bus_vm);
        // let diag_v_norm = CSR::with_diag((v / v.norm()).to_vec());
        let diag_v_norm = CSR::with_diagonal(v_norm);

        let d_imis_d_v1 = J * (y_bus * diag_v - diag_ibus); // dImis/dVa
        let d_imis_d_v2 = y_bus * diag_v_norm + diag_ibus_vm; // dImis/dVm

        (d_imis_d_v1, d_imis_d_v2)
    };

    Ok((d_imis_d_v1, d_imis_d_v2))
}

/// Computes the infinity norm: `max(abs(a))`
pub(crate) fn norm_inf(a: &[f64]) -> f64 {
    let mut max = f64::NEG_INFINITY;
    a.iter().for_each(|v| {
        let absvi = v.abs();
        if absvi > max {
            max = absvi
        }
    });
    max
}
