/// Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
/// by Ray Zimmerman, PSERC Cornell
/// Copyright (c) 2022, Richard Lincoln. All rights reserved.
use crate::mpc::*;
use crate::mpopt::MPOpt;
use crate::traits::Conj;
use ndarray::{arr1, Array1};
use num_complex::Complex64;
use sprs::{CsMat, CsMatBase, CsMatView, TriMat};
use std::f64::consts::PI;
use std::vec;

#[macro_export]
macro_rules! cmplx {
    () => {
        num_complex::Complex64::new(0.0, 0.0)
    };
    ($arg1:expr) => {
        num_complex::Complex64::new($arg1 as f64, 0.0)
    };
    ($arg1:expr, $arg2:expr) => {
        num_complex::Complex64::new($arg1 as f64, $arg2 as f64)
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
    gen.iter().for_each(|g| {
        if g.status {
            bus_gen_status[g.bus] = true
        }
    });

    // Form index lists for slack, PV, and PQ buses.
    // ref = find(bus(:, BUS_TYPE) == REF & bus_gen_status);   /// reference bus index
    // pv  = find(bus(:, BUS_TYPE) == PV  & bus_gen_status);   /// PV bus indices
    // pq  = find(bus(:, BUS_TYPE) == PQ | ~bus_gen_status);   /// PQ bus indices
    // let mut refbus = vec![0; bus.len()];
    // let mut pv = vec![0; bus.len()];
    // let mut pq = vec![0; bus.len()];
    // for (i, b) in bus.iter().enumerate() {
    //     match &b.bus_type {
    //         BusType::REF => {
    //             if bus_gen_status[i] {
    //                 refbus.push(b.i);
    //             }
    //         }
    //         BusType::PV => {
    //             if bus_gen_status[i] {
    //                 pv.push(b.i)
    //             }
    //         }
    //         bus_type => {
    //             if bus_type == BusType::PQ || !bus_gen_status[i] {
    //                 pq.push(b.i);
    //             }
    //         }
    //     }
    // }
    let refbus = bus
        .iter()
        .filter(|&b| b.bus_type == BusType::REF && bus_gen_status[b.i])
        .map(|b| b.i)
        .collect::<Vec<usize>>();
    let pv = bus
        .iter()
        .filter(|&b| b.bus_type == BusType::PV && bus_gen_status[b.i])
        .map(|b| b.i)
        .collect::<Vec<usize>>();
    let pq = bus
        .iter()
        .filter(|&b| b.bus_type == BusType::PQ && !bus_gen_status[b.i])
        .map(|b| b.i)
        .collect::<Vec<usize>>();

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
) -> (Array1<Complex64>, Option<CsMat<Complex64>>) {
    let nb = bus.len();

    // get load parameters
    let (sd_z, sd_i, sd_p) = make_sdzip(base_mva, bus, mpopt);

    // let sg_on = if let Some(sg) = sg {
    //     // gen.iter().filter(|&g| g.status == GenStatus::InService)
    //     sg.iter()
    //         .enumerate()
    //         .filter_map(|(&i, &v)| if gen[i].status { Some((i, v)) } else { None })
    //         .collect::<Vec<(usize, Complex64)>>()
    // } else {
    //     gen.iter()
    //         .filter(|&g| g.status)
    //         .map(|g| (g.bus, cmplx![g.pg, g.qg] / cmplx![baseMVA]))
    //         .collect::<Vec<(usize, Complex64)>>()
    // };

    let mut s_bus_g = Array1::<Complex64>::zeros(nb);
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
    let v_mag = match vm {
        Some(vm) => Array1::from(vm.iter().map(|&v| cmplx!(v)).collect::<Vec<Complex64>>()),
        None => Array1::<Complex64>::zeros(nb),
    };
    // let v_mag2 = v_mag.iter().map(|&v| cmplx!(v.re.powi(2)));

    // let s_bus_d = &sd_p + &sd_i * &v_mag + &sd_z * &v_mag2;
    let s_bus_d = &sd_p + &sd_i * &v_mag + &sd_z * &(&v_mag * &v_mag);

    // Form net complex bus power injection vector
    // (power injected by generators + power injected by loads).
    let s_bus = s_bus_g - s_bus_d;

    if !derivative {
        (s_bus, None)
    } else {
        let d_sbus_d_vm = if let Some(_) = vm {
            TriMat::from_triplets(
                (nb, nb),
                Vec::from_iter(0..nb),
                Vec::from_iter(0..nb),
                (&sd_i + cmplx!(2.0) * &v_mag * sd_z).to_vec(),
            )
            .to_csr()
        } else {
            CsMatBase::zero((nb, nb))
        };
        (s_bus, Some(d_sbus_d_vm))
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
) -> (Array1<Complex64>, Array1<Complex64>, Array1<Complex64>) {
    let pw = mpopt.exp.sys_wide_zip_loads.pw.unwrap_or([1.0, 0.0, 0.0]);
    let qw = mpopt.exp.sys_wide_zip_loads.pw.unwrap_or(pw);

    (
        Array1::from(
            bus.iter()
                .map(|b| cmplx!(b.pd * pw[1], b.qd * qw[1]) / cmplx!(base_mva))
                .collect::<Vec<Complex64>>(),
        ),
        Array1::from(
            bus.iter()
                .map(|b| cmplx!(b.pd * pw[1], b.qd * qw[1]) / cmplx!(base_mva))
                .collect::<Vec<Complex64>>(),
        ),
        Array1::from(
            bus.iter()
                .map(|b| cmplx!(b.pd * pw[0], b.qd * qw[0]) / cmplx!(base_mva))
                .collect::<Vec<Complex64>>(),
        ),
    )
}

/// Builds the bus admittance matrix and branch admittance matrices.
pub(crate) fn make_ybus(
    base_mva: f64,
    bus: &[Bus],
    branch: &[Branch],
    yf_yt: bool,
) -> (
    CsMat<Complex64>,
    Option<CsMat<Complex64>>,
    Option<CsMat<Complex64>>,
) {
    let nb = bus.len();
    let nl = branch.len();

    // For each branch, compute the elements of the branch admittance matrix where:
    //
    //      | If |   | Yff  Yft |   | Vf |
    //      |    | = |          | * |    |
    //      | It |   | Ytf  Ytt |   | Vt |
    /*
    {
        let mut y_bus = TriMat::new((nb, nb));
        let (mut y_f, mut yt) = if yf_yt = {
            (
                Some(TriMat::new((nl, nb))),
                Some(TriMat::new((nl, nb))),
            )
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

                y_f.add_triplet(i, f, y_ff);
                y_f.add_triplet(i, t, y_ft);

                y_t.add_triplet(i, f, y_tf);
                y_t.add_triplet(i, t, y_tt);
            }

            y_bus.add_triplet(f, f, y_ff);
            y_bus.add_triplet(f, t, y_ft);
            y_bus.add_triplet(t, f, y_tf);
            y_bus.add_triplet(t, t, y_tt);
        }

        for (i, b) in bus.iter().enumerate() {
            y_bus.add_triplet(i, i, b.y_sh(base_mva));
        }

        (y_bus, y_f, y_t)
    }
    */

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
    y_bus: CsMatView<Complex64>,
    v: &[Complex64],
    cartesian: bool,
) -> (CsMat<Complex64>, CsMat<Complex64>) {
    let i_bus: Array1<Complex64> = &y_bus * &arr1(v);

    let diag_v = spdiag(v);
    let diag_i_bus: CsMat<Complex64> = spdiag(&i_bus.to_vec());

    if cartesian {
        // dSbus/dVr = conj(diagIbus) + diagV * conj(Ybus)
        // dSbus/dVi = 1j * (conj(diagIbus) - diagV * conj(Ybus))
        (
            &diag_i_bus.conj() + &(&diag_v * &y_bus.to_owned().conj()),
            (&diag_i_bus.conj() - &(&diag_v * &y_bus.to_owned().conj())).map(|v| v * cmplx!(0, 1)),
        )
    } else {
        let diag_v_norm = spdiag(
            &v.iter()
                .map(|v| v / cmplx!(v.norm()))
                .collect::<Vec<Complex64>>(),
        );

        // dSbus/dVa = 1j * diagV * conj(diagIbus - Ybus * diagV)
        // dSbus/dVm = diagV * conj(Ybus * diagVnorm) + conj(diagIbus) * diagVnorm
        (
            (&diag_v * &(&diag_i_bus - &(&y_bus * &diag_v)).conj()).map(|v| v * cmplx!(0, 1)),
            &(&diag_v * &(&y_bus * &diag_v_norm).conj()) + &(&diag_i_bus.conj() * &diag_v_norm),
        )
    }
}

fn spdiag(d: &[Complex64]) -> CsMat<Complex64> {
    let n = d.len();
    let tri = TriMat::from_triplets((n, n), (0..n).collect(), (0..n).collect(), d.to_vec());
    tri.to_csr()
}

fn sparse(
    shape: (usize, usize),
    row_inds: &[usize],
    col_inds: &[usize],
    data: &[Complex64],
) -> CsMat<Complex64> {
    let tri = TriMat::from_triplets(shape, row_inds.to_vec(), col_inds.to_vec(), data.to_vec());
    tri.to_csr()
}
