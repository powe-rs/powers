/// Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
/// by Ray Zimmerman, PSERC Cornell
/// Copyright (c) 2022, Richard Lincoln. All rights reserved.
use crate::mpc::*;
use crate::mpopt::MPOpt;
use densetools::arr::{Arr, CArr};
use num_complex::Complex64;
use sparsetools::coo::Coo;
use sparsetools::csr::{conj, CSR};
use std::f64::consts::PI;
use std::vec;

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
fn bus_types(bus: &[Bus], gen: &[Gen]) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
    // Buses with generators that are ON.
    let mut bus_gen_status = vec![false; bus.len()];
    for g in gen {
        if g.status {
            bus_gen_status[g.bus] = true
        }
    }

    // Form index lists for slack, PV, and PQ buses.
    // ref = find(bus(:, BUS_TYPE) == REF & bus_gen_status);   /// reference bus index
    // pv  = find(bus(:, BUS_TYPE) == PV  & bus_gen_status);   /// PV bus indices
    // pq  = find(bus(:, BUS_TYPE) == PQ | ~bus_gen_status);   /// PQ bus indices
    let mut refbus = Vec::new();
    let mut pv = Vec::new();
    let mut pq = Vec::new();
    for (i, b) in bus.iter().enumerate() {
        match &b.bus_type {
            BusType::REF => {
                if bus_gen_status[i] {
                    refbus.push(b.i);
                }
            }
            BusType::PV => {
                if bus_gen_status[i] {
                    pv.push(b.i)
                }
            }
            &bus_type => {
                if bus_type == BusType::PQ || !bus_gen_status[i] {
                    pq.push(b.i);
                }
            }
        }
    }

    (refbus, pv, pq)
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
/// With `derivative` true, it computes the partial derivative of the
/// bus injections with respect to voltage magnitude, leaving the first
/// return value SBUS empty. If VM is empty, it assumes no voltage dependence
/// and returns a sparse zero matrix.
///
/// See also MAKEYBUS.
pub(crate) fn make_sbus(
    base_mva: f64,
    bus: &[Bus],
    gen: &[Gen],
    mpopt: &MPOpt,
    vm: Option<&[f64]>,
    sg: Option<&[Complex64]>,
    derivative: bool,
) -> (Arr<Complex64>, Option<Coo<usize, Complex64>>) {
    let nb = bus.len();

    // get load parameters
    let (sd_z, sd_i, sd_p) = make_sdzip(base_mva, bus, mpopt);

    let mut s_bus_g = Arr::<Complex64>::zeros(nb);
    // for (gen_bus, sg_i) in sg_on {
    //     sbusg[gen_bus] += sg_i;
    // }
    if let Some(sg) = sg {
        for (g, s_pu) in gen.iter().zip(sg) {
            if g.status {
                s_bus_g[g.bus] += s_pu
            }
        }
    } else {
        for g in gen {
            if g.status {
                let s_pu = cmplx![g.pg, g.qg] / cmplx![base_mva];
                s_bus_g[g.bus] += s_pu;
            }
        }
    }

    // Compute per-bus loads in p.u.
    let v_mag: Arr<Complex64> = match vm {
        Some(vm) => Arr::from_real(vm),
        None => Arr::zeros(nb),
    };

    let s_bus_d = &sd_p + &sd_i * &v_mag + &sd_z * (&v_mag * &v_mag);

    // Form net complex bus power injection vector
    // (power injected by generators + power injected by loads).
    let s_bus = s_bus_g - s_bus_d;

    if !derivative {
        return (s_bus, None);
    }

    let d_sbus_d_vm = if let Some(_) = vm {
        Coo::with_diagonal(&(sd_i + cmplx!(2.0) * v_mag * sd_z))
    } else {
        // CsMatBase::zero((nb, nb))
        Coo::empty(nb, nb, 0)
    };
    (s_bus, Some(d_sbus_d_vm))
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
) -> (Arr<Complex64>, Arr<Complex64>, Arr<Complex64>) {
    let pw = mpopt.exp.sys_wide_zip_loads.pw.unwrap_or([1.0, 0.0, 0.0]);
    let qw = mpopt.exp.sys_wide_zip_loads.qw.unwrap_or(pw);

    let mut sd_z = Arr::new();
    let mut sd_i = Arr::new();
    let mut sd_p = Arr::new();
    for (i, b) in bus.iter().enumerate() {
        sd_z[i] = cmplx!(b.pd * pw[2], b.qd * qw[2]) / cmplx!(base_mva);
        sd_i[i] = cmplx!(b.pd * pw[1], b.qd * qw[1]) / cmplx!(base_mva);
        sd_p[i] = cmplx!(b.pd * pw[0], b.qd * qw[0]) / cmplx!(base_mva);
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
    let mut y_bus = Coo::empty(nb, nb, 0);
    let (mut y_f, mut y_t) = if yf_yt {
        (Some(Coo::empty(nl, nb, 0)), Some(Coo::empty(nl, nb, 0)))
    } else {
        (None, None)
    };

    for (i, br) in branch.iter().enumerate() {
        let y_s = br.y_s(); // series admittance
        let b_c = if br.status { br.b } else { 0.0 }; // line charging susceptance
        let t = if br.tap == 0.0 { 1.0 } else { br.tap }; // default tap ratio = 1
        let tap = Complex64::from_polar(t, br.shift * PI / 180.0); // add phase shifters

        let y_tt = y_s + cmplx!(0.0, b_c / 2.0);
        let y_ff = y_tt / (tap * tap.conj());
        let y_ft = -y_s / tap.conj();
        let y_tf = -y_s / tap;

        let (f, t) = (br.from_bus, br.to_bus);

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

    for (i, b) in bus.iter().enumerate() {
        y_bus.push(i, i, b.y_sh(base_mva));
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
    v: &Arr<Complex64>,
    cartesian: bool,
) -> (CSR<usize, Complex64>, CSR<usize, Complex64>) {
    const J: Complex64 = Complex64 { re: 0.0, im: 1.0 };

    let i_bus = y_bus * v;

    let diag_v = CSR::<usize, Complex64>::with_diag(v.to_vec());
    let diag_i_bus = CSR::<usize, Complex64>::with_diag(i_bus.to_vec());

    if cartesian {
        // dSbus/dVr = conj(diagIbus) + diagV * conj(Ybus)
        // dSbus/dVi = 1j * (conj(diagIbus) - diagV * conj(Ybus))

        let d_sbus_d_vr = conj(&diag_i_bus) + &diag_v * conj(y_bus);
        let d_sbus_d_vi = (conj(&diag_i_bus) - &diag_v * conj(y_bus)) * J;

        (d_sbus_d_vr, d_sbus_d_vi)
    } else {
        let v_norm = v.iter().map(|v| v / cmplx!(v.norm())).collect();
        let diag_v_norm = CSR::<usize, Complex64>::with_diag(v_norm);

        // dSbus/dVa = 1j * diagV * conj(diagIbus - Ybus * diagV)
        // dSbus/dVm = diagV * conj(Ybus * diagVnorm) + conj(diagIbus) * diagVnorm

        let d_sbus_d_va = &diag_v * conj(&(&diag_i_bus - y_bus * &diag_v)) * J;
        let d_sbus_d_vm =
            &diag_v * conj(&(y_bus * &diag_v_norm)) + conj(&diag_i_bus) * &diag_v_norm;

        (d_sbus_d_va, d_sbus_d_vm)
    }
}
