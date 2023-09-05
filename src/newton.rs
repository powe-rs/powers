use crate::mpopt::MPOpt;
use crate::powers::{d_imis_d_v, d_sbus_d_v, norm_inf, SBus, J};
use crate::traits::LinearSolver;

use itertools::{izip, Itertools};
use num_complex::Complex64;
use sparsetools::coo::Coo;
use sparsetools::csr::{CCSR, CSR};

pub trait ProgressMonitor {
    fn update(&self, i: usize, norm_f: f64);
}

/// Solves power flow using full Newton's method (power/polar).
///
/// Solves for bus voltages using a full Newton-Raphson method, using nodal
/// power balance equations and polar coordinate representation of
/// voltages.
///
/// The bus voltage vector contains the set point for generator
/// (including ref bus) buses, and the reference angle of the swing
/// bus, as well as an initial guess for remaining magnitudes and
/// angles.
///
/// Returns the final complex voltages, a flag which indicates whether it
/// converged or not, and the number of iterations performed.
pub(crate) fn newtonpf(
    y_bus: &CSR<usize, Complex64>,
    s_bus_fn: &dyn SBus,
    v0: &[Complex64],
    _ref: &[usize],
    pv: &[usize],
    pq: &[usize],
    lin_solver: &dyn LinearSolver,
    mpopt: &MPOpt,
    progress: Option<&dyn ProgressMonitor>,
) -> Result<(Vec<Complex64>, bool, usize), String> {
    let nb = v0.len();
    let pv_pq = [pv, pq].concat();

    let tol = mpopt.pf.tolerance;
    let max_it = mpopt.pf.max_it_nr;
    // let lin_solver  = mpopt.pf.nr.lin_solver;

    let mut converged = false;
    let mut i = 0;
    let mut v: Vec<Complex64> = v0.to_vec();
    let mut va = v.iter().map(|v| v.arg()).collect_vec();
    let mut vm = v.iter().map(|v| v.norm()).collect_vec();

    // set up indexing for updating V
    let npv = pv.len();
    let npq = pq.len();
    let (j1, j2) = (0, npv); // j1:j2 - V angle of pv buses
    let (j3, j4) = (j2, j2 + npq); // j3:j4 - V angle of pq buses
    let (j5, j6) = (j4, j4 + npq); // j5:j6 - V mag of pq buses

    // let j1_j2 = j1..j2;
    // let j3_j4 = j3..j4;
    // let j5_j6 = j5..j6;

    // let j1_j2 = 0..npv; // V angle of pv buses
    // let j3_j4 = j1_j2.end..(j1_j2.end + npq); // V angle of pq buses
    // let j5_j6 = j3_j4.end..(j3_j4.end + npq); // V mag of pq buses

    // evaluate F(x0)
    let i_bus: Vec<Complex64> = y_bus * &v;
    let s_bus = s_bus_fn.s_bus(&vm);
    // let mis = &v * (y_bus * &v).conj() - s_bus.s_bus(&vm);
    // let mis: Vec<Complex64> = (0..nb).map(|i| v[i] * i_bus[i].conj() - s_bus[i]).collect();
    let mis = izip!(&v, &i_bus, &s_bus)
        .map(|(v, i_bus, s_bus)| v * i_bus.conj() - s_bus)
        .collect_vec();
    let f = [
        pv_pq.iter().map(|&i| mis[i].re).collect_vec(),
        pq.iter().map(|&i| mis[i].im).collect_vec(),
    ]
    .concat();

    // check tolerance
    let norm_f = norm_inf(&f);
    if let Some(pm) = progress {
        pm.update(i, norm_f);
    }
    if norm_f < tol {
        converged = true;
        println!("Converged!");
    }

    // do Newton iterations
    while !converged && i < max_it {
        // update iteration counter
        i = i + 1;

        // evaluate Jacobian
        let (d_sbus_d_va, mut d_sbus_d_vm) = d_sbus_d_v(y_bus, &v, false);
        let neg_d_sd_d_vm = s_bus_fn.d_sbus_d_vm(&vm);
        d_sbus_d_vm = d_sbus_d_vm - neg_d_sd_d_vm;

        let j11 = d_sbus_d_va.select(Some(&pv_pq), Some(&pv_pq))?.real();
        let j12 = d_sbus_d_vm.select(Some(&pv_pq), Some(pq))?.real();
        let j21 = d_sbus_d_va.select(Some(pq), Some(&pv_pq))?.imag();
        let j22 = d_sbus_d_vm.select(Some(pq), Some(pq))?.imag();

        let jac = Coo::compose([
            [&j11.to_coo(), &j12.to_coo()],
            [&j21.to_coo(), &j22.to_coo()],
        ])?;

        // compute update step
        let neg_f: Vec<f64> = f.iter().map(|f| -f).collect();
        let dx = lin_solver.solve(jac.to_csc(), &neg_f)?;

        // update voltage
        // pv.iter().zip(j1..j2).for_each(|(&i, j)| va[i] += dx[j]);
        // pq.iter().zip(j3..j4).for_each(|(&i, j)| va[i] += dx[j]);
        // pq.iter().zip(j5..j6).for_each(|(&i, j)| vm[i] += dx[j]);
        for (i, j) in (j1..j2).enumerate() {
            va[pv[i]] += dx[j];
        }
        for (i, j) in (j3..j4).enumerate() {
            va[pq[i]] += dx[j];
        }
        for (i, j) in (j5..j6).enumerate() {
            vm[pq[i]] += dx[j];
        }

        // update Vm and Va again in case we wrapped around with a negative Vm
        // v = Vec::from_polar(&vm, &va);
        v = izip!(vm, va)
            .map(|(vm, va)| Complex64::from_polar(vm, va))
            .collect();
        va = v.iter().map(|v| v.arg()).collect();
        vm = v.iter().map(|v| v.norm()).collect();

        // evalute F(x)
        let i_bus: Vec<Complex64> = y_bus * &v;
        let s_bus = s_bus_fn.s_bus(&vm);
        // let mis = &v * conj(&(y_bus * &v)) - s_bus.s_bus(&vm);
        // let mis: Vec<Complex64> = (0..nb).map(|i| v[i] * i_bus[i].conj() - s_bus[i]).collect();
        let mis = izip!(&v, &i_bus, &s_bus)
            .map(|(v, i_bus, s_bus)| v * i_bus.conj() - s_bus)
            .collect_vec();
        let f = [
            pv_pq.iter().map(|&i| mis[i].re).collect_vec(),
            pq.iter().map(|&i| mis[i].im).collect_vec(),
        ]
        .concat();

        // check for convergence
        let norm_f = norm_inf(&f);
        if let Some(pm) = progress {
            pm.update(i, norm_f);
        }
        if norm_f < tol {
            converged = true;
            println!(
                "Newton's method power flow (power balance, polar) converged in {} iterations.",
                i
            );
        }
    }

    if !converged {
        println!(
            "Newton's method power flow (power balance, polar) did not converge in {} iterations.",
            i
        );
    }

    Ok((v, converged, i))
}

/// Solves power flow using full Newton's method (current/cartesian).
///
/// Solves for bus voltages using a full Newton-Raphson method, using nodal
/// current balance equations and polar coordinate representation of
/// voltages.
///
/// The bus voltage vector contains the set point for generator
/// (including ref bus) buses, and the reference angle of the swing
/// bus, as well as an initial guess for remaining magnitudes and
/// angles.
///
/// Returns the final complex voltages, a flag which indicates whether it
/// converged or not, and the number of iterations performed.
pub(crate) fn newtonpf_i_polar(
    y_bus: &CSR<usize, Complex64>,
    s_bus_fn: &dyn SBus,
    v0: &[Complex64],
    ref_: &[usize],
    pv: &[usize],
    pq: &[usize],
    lin_solver: &dyn LinearSolver,
    mpopt: &MPOpt,
    progress: Option<&dyn ProgressMonitor>,
) -> Result<(Vec<Complex64>, bool, usize), String> {
    let pv_pq = [pv, pq].concat();

    let tol = mpopt.pf.tolerance;
    let max_it = mpopt.pf.max_it_nr;

    let mut converged = false;
    let mut i = 0;
    let mut v = v0.to_vec();
    let mut va: Vec<f64> = v.iter().map(|v| v.arg()).collect();
    let mut vm: Vec<f64> = v.iter().map(|v| v.norm()).collect();
    // let vm_pv: Vec<f64> = pv.iter().map(|&i| vm[i]).collect();
    let nb = v0.len();

    // set up indexing for updating V
    let npv = pv.len();
    let npq = pq.len();
    let (j1, j2) = (0, npv); // j1:j2 - V angle of pv buses
    let (j3, j4) = (j2, j2 + npq); // j3:j4 - V angle of pq buses
    let (j5, j6) = (j4, j4 + npv); // j5:j6 - Q of pv buses
    let (j7, j8) = (j6, j6 + npv); // j7:j8 - V mag of pq buses
    let j1_j2: Vec<usize> = (j1..j2).collect();
    let j3_j4: Vec<usize> = (j3..j4).collect();
    let j5_j6: Vec<usize> = (j5..j6).collect();
    let j7_j8: Vec<usize> = (j7..j8).collect();

    // evaluate F(x0)
    let mut s_bus = s_bus_fn.s_bus(&vm);
    let mis: Vec<Complex64> = {
        // let ybus_pv = y_bus.select(Some(&pv), None)?;
        // let i_bus_pv: Vec<Complex64> = ybus_pv * &v;
        let i_bus: Vec<Complex64> = y_bus * &v;
        // for i in 0..npv {
        pv.iter().for_each(|&i| {
            s_bus[i].im = (v[i] * i_bus[i]).conj().im;
        });
        // (0..nb)
        //     .map(|i| i_bus[i] - (s_bus[i] / v[i]).conj())
        //     .collect()
        izip!(&i_bus, &s_bus, &v)
            .map(|(i_bus, s_bus, v)| i_bus - (s_bus / v).conj())
            .collect_vec()
    };
    let f = [
        pv_pq.iter().map(|&i| mis[i].re).collect_vec(),
        pv_pq.iter().map(|&i| mis[i].im).collect_vec(),
    ]
    .concat();

    // check tolerance
    let norm_f = norm_inf(&f);
    if let Some(pm) = progress {
        pm.update(i, norm_f); // max Ir & Ii mismatch (p.u.)
    }
    if norm_f < tol {
        converged = true;
        println!("Converged!");
    }

    fn compose(j: [[&Coo<usize, f64>; 3]; 2]) -> Result<Coo<usize, f64>, String> {
        let j1x = Coo::h_stack(j[0][0], j[0][1], Some(j[0][2]))?;
        let j2x = Coo::h_stack(j[1][0], j[1][1], Some(j[1][2]))?;
        Coo::v_stack(&j1x, &j2x, None)
    }

    // do Newton iterations
    while !converged && i < max_it {
        // update iteration counter
        i = i + 1;

        // evaluate Jacobian
        /*
        dImis_dQ = sparse(pv, pv, 1j./conj(V(pv)), n, n);
        [dImis_dVa, dImis_dVm] = dImis_dV(Sb, Ybus, V);
        dImis_dVm(:, pv) = dImis_dQ(:, pv);

        j11 = real(dImis_dVa([pv; pq], [pv; pq]));
        j12 = real(dImis_dVm([pv; pq], [pv; pq]));
        j21 = imag(dImis_dVa([pv; pq], [pv; pq]));
        j22 = imag(dImis_dVm([pv; pq], [pv; pq]));

        J = [   j11 j12;
                j21 j22;    ];
        */
        let dImis_dQ = Coo::new(
            nb,
            nb,
            pv.to_vec(),
            pv.to_vec(),
            pv.iter().map(|&i| J / v[i].conj()).collect(),
        )?
        .to_csr();
        let (dImis_dVa, dImis_dVm) = d_imis_d_v(&s_bus, &y_bus, &v, false)?;

        /*
        j11 = real(dImis_dVa([pv; pq], [pv; pq]));
        j12 = real(dImis_dQ([pv; pq], pv));
        j13 = real(dImis_dVm([pv; pq], pq));
        j21 = imag(dImis_dVa([pv; pq], [pv; pq]));
        j22 = imag(dImis_dQ([pv; pq], pv));
        j23 = imag(dImis_dVm([pv; pq], pq));

        J = [   j11 j12 j13;
                j21 j22 j23;    ];
        */
        let j11 = dImis_dVa.select(Some(&pv_pq), Some(&pv_pq))?.real();
        let j12 = dImis_dQ.select(Some(&pv_pq), Some(pv))?.real();
        let j13 = dImis_dVm.select(Some(&pv_pq), Some(pq))?.real();
        let j21 = dImis_dVa.select(Some(&pv_pq), Some(&pv_pq))?.imag();
        let j22 = dImis_dQ.select(Some(&pv_pq), Some(pv))?.imag();
        let j23 = dImis_dVm.select(Some(&pv_pq), Some(pq))?.imag();

        let jac = compose([
            [&j11.to_coo(), &j12.to_coo(), &j13.to_coo()],
            [&j21.to_coo(), &j22.to_coo(), &j23.to_coo()],
        ])?;

        // compute update step
        // let dx = lin_solver.solve(jac.view(), f.iter().map(|&f_i| -f_i).collect())?;
        let neg_f: Vec<f64> = f.iter().map(|f| -f).collect();
        let dx = lin_solver.solve(jac.to_csc(), &neg_f)?;

        // for (i, j) in (j1..j2).enumerate() {
        //     va[pv[i]] += dx[j];
        // }
        // for (i, j) in (j3..j4).enumerate() {
        //     va[pq[i]] += dx[j];
        // }
        // for (i, j) in (j5..j6).enumerate() {
        //     vm[pq[i]] += dx[j];
        // }

        // update voltage
        if npv != 0 {
            /*
            Va(pv) = Va(pv) + dx(j1:j2);
            Sb(pv) = real(Sb(pv)) + 1j * (imag(Sb(pv)) + dx(j5:j6));
            */
            // let dx_va = dx.select(&(j1..j2).collect::<Vec<usize>>());
            // let dx_sb = dx.select(&(j5..j6).collect::<Vec<usize>>());
            // for (i, &j) in pv.iter().enumerate() {
            //     va[j] += dx_va[i];
            //     s_bus[j] += Complex64::new(0.0, dx_sb[i]);
            // }

            (0..npv).for_each(|i| va[pv[i]] += dx[j1_j2[i]]);
            (0..npv).for_each(|i| s_bus[pv[i]].im += dx[j5_j6[i]]);
            // for i in 0..npv {
            //     va[pv[i]] += dx[j1_j2[i]];
            //     s_bus[pv[i]].im += dx[j5_j6[i]];
            // }
        }
        if npq != 0 {
            /*
            Va(pq) = Va(pq) + dx(j3:j4);
            Vm(pq) = Vm(pq) + dx(j7:j8);
            */
            (0..npq).for_each(|i| va[pq[i]] += dx[j3_j4[i]]);
            (0..npq).for_each(|i| vm[pq[i]] += dx[j7_j8[i]]);
            // for i in 0..npq {
            //     va[pq[i]] += dx[j3_j4[i]];
            //     vm[pq[i]] += dx[j7_j8[i]];
            // }
        }

        // update Vm and Va again in case we wrapped around with a negative Vm
        // v = Vec::from_polar(&vm, &va);
        v = izip!(vm, va)
            .map(|(vm, va)| Complex64::from_polar(vm, va))
            .collect();
        va = v.iter().map(|v| v.arg()).collect();
        vm = v.iter().map(|v| v.norm()).collect();

        // evalute F(x)
        // mis = Ybus * V - conj(Sb ./ V);
        let i_bus = y_bus * &v;
        let mis = izip!(&i_bus, &s_bus, &v)
            .map(|(i_bus, s_bus, v)| i_bus - (s_bus / v).conj())
            .collect_vec();
        // let mis: Vec<Complex64> = (0..nb)
        //     .map(|i| i_bus[i] - (s_bus[i] / v[i]).conj())
        //     .collect();

        let f = [
            pv_pq.iter().map(|&i| mis[i].re).collect::<Vec<f64>>(),
            pv_pq.iter().map(|&i| mis[i].im).collect::<Vec<f64>>(),
        ]
        .concat();

        // check for convergence
        let norm_f = norm_inf(&f);
        if let Some(pm) = progress {
            pm.update(i, norm_f);
        }
        if norm_f < tol {
            converged = true;
            println!(
                "Newton's method power flow (current balance, polar) converged in {} iterations.",
                i
            );
        }
    }

    if !converged {
        println!(
            "Newton's method power flow (current balance, polar) did not converge in {} iterations.",
            i
        );
    }

    Ok((v, converged, i))
}

/// Solves power flow using full Newton's method (current/cartesian).
///
/// Solves for bus voltages using a full Newton-Raphson method, using nodal
/// current balance equations and cartesian coordinate representation of
/// voltages.
///
/// The bus voltage vector contains the set point for generator
/// (including ref bus) buses, and the reference angle of the swing
/// bus, as well as an initial guess for remaining magnitudes and
/// angles.
///
/// Returns the final complex voltages, a flag which indicates whether it
/// converged or not, and the number of iterations performed.
pub(crate) fn newtonpf_i_cart(
    y_bus: &CSR<usize, Complex64>,
    s_bus_fn: &dyn SBus,
    v0: &[Complex64],
    ref_: &[usize],
    pv: &[usize],
    pq: &[usize],
    lin_solver: &dyn LinearSolver,
    mpopt: &MPOpt,
    progress: Option<&dyn ProgressMonitor>,
) -> Result<(Vec<Complex64>, bool, usize), String> {
    let pv_pq = [pv, pq].concat();

    let tol = mpopt.pf.tolerance;
    let max_it = mpopt.pf.max_it_nr;

    let mut converged = false;
    let mut i = 0;
    let mut v = v0.to_vec();
    // let mut va: Vec<f64> = v.iter().map(|v| v.arg()).collect();
    let mut vm: Vec<f64> = v.iter().map(|v| v.norm()).collect();
    // let vm_pv = pv.iter().map(|&i| vm[i]).collect();
    let nb = v0.len();

    // set up indexing for updating V
    let npv = pv.len();
    let npq = pq.len();
    let (j1, j2) = (0, npv); // j1:j2 - Q of pv buses
    let (j3, j4) = (j2, j2 + npq); // j3:j4 - Vr of pq buses
    let (j5, j6) = (j4, j4 + npv); // j5:j6 - Vr of pv buses
    let (j7, j8) = (j6, j6 + npq); // j7:j9 - Vi of pq buses
    let (j9, j10) = (j8, j8 + npv); // j9:j10 - Vi of pv buses
    let j1_j2: Vec<usize> = (j1..j2).collect();
    let j3_j4: Vec<usize> = (j3..j4).collect();
    let j5_j6: Vec<usize> = (j5..j6).collect();
    let j7_j8: Vec<usize> = (j7..j8).collect();
    let j9_j10: Vec<usize> = (j9..j10).collect();

    // evaluate F(x0)
    let mut s_bus = s_bus_fn.s_bus(&vm);
    let mis = {
        // let ybus_pv = y_bus.select(Some(&pv), None)?;
        // let i_bus_pv: Vec<Complex64> = ybus_pv * &v;
        // for i in 0..npv {
        //     s_bus[pv[i]].im = (v[pv[i]] * i_bus_pv[i]).conj().im;
        // }
        let i_bus: Vec<Complex64> = y_bus * &v;
        pv.iter()
            .for_each(|&i| s_bus[i].im = (v[i] * i_bus[i]).conj().im);

        izip!(&i_bus, &s_bus, &v)
            .map(|(i_bus, s_bus, v)| i_bus - (s_bus / v).conj())
            .collect_vec()
    };

    let f = [
        pv_pq.iter().map(|&i| mis[i].re).collect_vec(),
        pv_pq.iter().map(|&i| mis[i].im).collect_vec(),
        pv.iter()
            .map(|&i| (v[i] * v[i].conj()).re - (vm[i] * vm[i]))
            .collect_vec(),
    ]
    .concat();

    // let sb_pv = Vec::from_parts(
    //     &s_bus.select(&pv).real(),
    //     &(v.select(&pv) * conj(&(y_bus * &v))).imag(),
    // );
    // s_bus.set(pv, &sb_pv);
    // let mis = (y_bus * &v) - conj(&(&s_bus / &v));
    //
    //     let mut f = {
    //         let mis_pvpq = mis.select(&pv_pq);
    //         let v_pv = v.select(&pv);
    //         Vec::concat(&[
    //             &mis_pvpq.real(),
    //             &mis_pvpq.imag(),
    //             &(&(&v_pv * conj(&v_pv)).real() - vm_pv.pow(2)),
    //         ])
    //     };

    // check tolerance
    let norm_f = norm_inf(&f);
    if let Some(pm) = progress {
        pm.update(i, norm_f); // max Ir & Ii mismatch (p.u.)
    }
    if norm_f < tol {
        converged = true;
        println!("Converged!"); // TODO: verbose feature
    }

    // do Newton iterations
    while !converged && i < max_it {
        // update iteration counter
        i = i + 1;

        // evaluate Jacobian
        /*
        dImis_dQ = sparse(pv, pv, 1j./conj(V(pv)), n, n);
        dV2_dVr = sparse(1:npv, npq+(1:npv), 2*real(V(pv)), npv, npv+npq);
        dV2_dVi = sparse(1:npv, npq+(1:npv), 2*imag(V(pv)), npv, npv+npq);
        [dImis_dVr, dImis_dVi] = dImis_dV(Sb, Ybus, V, 1);
             */
        // let v_pv = v.select(&pv);
        let dImis_dQ = Coo::new(
            nb,
            nb,
            pv.to_vec(),
            pv.to_vec(),
            pv.iter().map(|&i| J / v[i].conj()).collect(),
        )?
        .to_csr();
        let dV2_dVr = Coo::new(
            npv,
            npv + npq,
            (0..npv).collect(),
            (npq..npq + npv).collect(),
            pv.iter().map(|&i| 2.0 * v[i].re).collect(),
        )?
        .to_csr();
        let dV2_dVi = Coo::new(
            npv,
            npv + npq,
            (0..npv).collect(),
            (npq..npq + npv).collect(),
            pv.iter().map(|&i| 2.0 * v[i].im).collect(),
        )?
        .to_csr();
        let (dImis_dVr, dImis_dVi) = d_imis_d_v(&s_bus, &y_bus, &v, true)?;

        // handling of derivatives for voltage dependent loads
        // (not yet implemented) goes here

        let j11 = dImis_dQ.select(Some(pv), Some(pv))?.real();
        let j12 = dImis_dVr.select(Some(&pv_pq), Some(&pv_pq))?.real();
        let j13 = dImis_dVi.select(Some(&pv_pq), Some(&pv_pq))?.real();

        let j21 = dImis_dQ.select(Some(&pv_pq), Some(&pv))?.imag();
        let j22 = dImis_dVr.select(Some(&pv_pq), Some(&pv_pq))?.imag();
        let j23 = dImis_dVi.select(Some(&pv_pq), Some(&pv_pq))?.imag();

        let j31 = CSR::zeros(npv, npv);
        let j32 = dV2_dVr;
        let j33 = dV2_dVi;
        /*
        j11 = real(dImis_dQ([pv; pq], pv));
        j12 = real(dImis_dVr([pv; pq], [pq; pv]));
        j13 = real(dImis_dVi([pv; pq], [pq; pv]));

        j21 = imag(dImis_dQ([pv; pq], pv));
        j22 = imag(dImis_dVr([pv; pq], [pq; pv]));
        j23 = imag(dImis_dVi([pv; pq], [pq; pv]));

        j31 = sparse(npv, npv);
        j32 = dV2_dVr;
        j33 = dV2_dVi;
        */

        let jac = Coo::compose3([
            [&j11.to_coo(), &j12.to_coo(), &j13.to_coo()],
            [&j21.to_coo(), &j22.to_coo(), &j23.to_coo()],
            [&j31.to_coo(), &j32.to_coo(), &j33.to_coo()],
        ])?;

        // compute update step
        let neg_f = f.iter().map(|f| -f).collect_vec();
        let dx: Vec<f64> = lin_solver.solve(jac.to_csc(), &neg_f)?;

        // update voltage
        if npv != 0 {
            /*
            V(pv) = V(pv) + dx(j5:j6) + 1j * dx(j9:j10);
            Sb(pv) = real(Sb(pv)) + 1j * (imag(Sb(pv)) + dx(j1:j2));
            */
            // let dx_v: Vec<Complex64> = Vec::from_parts(
            //     &dx.select(&(j5..j6).collect::<Vec<usize>>()),
            //     &dx.select(&(j9..j10).collect::<Vec<usize>>()),
            // );
            // let dx_sb = dx.select(&(j1..j2).collect::<Vec<usize>>());
            // for (i, &j) in pv.iter().enumerate() {
            //     v[j] += dx_v[i];
            //     s_bus[j] += Complex64::new(0.0, dx_sb[i]);
            // }

            for i in 0..npv {
                v[pv[i]] += Complex64::new(dx[j5_j6[i]], dx[j9_j10[i]]);
                s_bus[pv[i]].im += dx[j1_j2[i]];
            }
        }
        if npq != 0 {
            // V(pq) = V(pq) + dx(j3:j4) + 1j * dx(j7:j8);
            // let dx_v = Vec::<Complex64>::from_parts(
            //     &dx.select(&(j3..j4).collect::<Vec<usize>>()),
            //     &dx.select(&(j7..j8).collect::<Vec<usize>>()),
            // );
            // for (i, &j) in pq.iter().enumerate() {
            //     v[j] += dx_v[i];
            // }

            for i in 0..npq {
                v[pq[i]] += Complex64::new(dx[j3_j4[i]], dx[j7_j8[i]]);
            }
        }

        // evalute F(x)
        /*
        mis = Ybus * V - conj(Sb ./ V);
        F = [   real(mis([pv; pq]));
                imag(mis([pv; pq]));
                V(pv) .* conj(V(pv)) - Vmpv.^2  ];
        */
        let i_bus = y_bus * &v;
        // let mis = (0..nb)
        //     .map(|i| i_bus[i] - (s_bus[i] / v[i]).conj())
        //     .collect::<Vec<Complex64>>();
        let mis = izip!(&i_bus, &s_bus, &v)
            .map(|(i_bus, s_bus, v)| i_bus - (s_bus / v).conj())
            .collect_vec();
        let f = [
            pv_pq.iter().map(|&i| mis[i].re).collect_vec(),
            pv_pq.iter().map(|&i| mis[i].im).collect_vec(),
            pv.iter()
                .map(|&i| (v[i] * v[i].conj()).re - (vm[i] * vm[i]))
                .collect_vec(),
        ]
        .concat();

        // check for convergence
        let norm_f = norm_inf(&f);
        if let Some(pm) = progress {
            pm.update(i, norm_f);
        }
        if norm_f < tol {
            converged = true;
            // TODO: verbose feature
            println!(
                "Newton's method power flow (current balance, cartesian) converged in {} iterations.",
                i
            );
        }
    }

    if !converged {
        println!(
            "Newton's method power flow (current balance, cartesian) did not converge in {} iterations.",
            i
        );
    }

    Ok((v, converged, i))
}

/// Solves power flow using full Newton's method (current/hybrid).
///
/// Solves for bus voltages using a full Newton-Raphson method, using nodal
/// current balance equations and a hybrid representation of voltages, where
/// a polar update is computed using a cartesian Jacobian.
///
/// The bus voltage vector contains the set point for generator
/// (including ref bus) buses, and the reference angle of the swing
/// bus, as well as an initial guess for remaining magnitudes and
/// angles.
///
/// Returns the final complex voltages, a flag which indicates whether it
/// converged or not, and the number of iterations performed.
pub(crate) fn newtonpf_I_hybrid(
    y_bus: &CSR<usize, Complex64>,
    s_bus_fn: &dyn SBus,
    v0: &[Complex64],
    _ref: &[usize],
    pv: &[usize],
    pq: &[usize],
    lin_solver: &dyn LinearSolver,
    mpopt: &MPOpt,
    progress: Option<&dyn ProgressMonitor>,
) -> Result<(Vec<Complex64>, bool, usize), String> {
    let pv_pq = [pv, pq].concat();

    let tol = mpopt.pf.tolerance;
    let max_it = mpopt.pf.max_it_nr;

    let mut converged = false;
    let mut i = 0;
    let mut v = v0.to_vec();
    let mut va = v.iter().map(|v| v.arg()).collect_vec();
    let mut vm = v.iter().map(|v| v.norm()).collect_vec();
    // let vm_pv = vm.select(pv);
    let nb = v0.len();

    // set up indexing for updating V
    let npv = pv.len();
    let npq = pq.len();
    let (j1, j2) = (0, npv); // j1:j2 - Q of pv buses
    let (j3, j4) = (j2, j2 + npq); // j3:j4 - Vr of pq buses
    let (j5, j6) = (j4, j4 + npv); // j5:j6 - Vi of pv buses
    let (j7, j8) = (j6, j6 + npq); // j7:j8 - Vi of pq buses

    // evaluate F(x0)
    let mut s_bus = s_bus_fn.s_bus(&vm);
    // Sb(pv) = real(Sb(pv)) + 1j * imag(V(pv) .* conj(Ybus(pv, :) * V));
    let mis = {
        // let ybus_pv = y_bus.select(Some(&pv), None)?;
        // let i_bus_pv: Vec<Complex64> = ybus_pv * &v;
        // for i in 0..npv {
        //     s_bus[pv[i]].im = (v[pv[i]] * i_bus_pv[i]).conj().im;
        // }
        let i_bus: Vec<Complex64> = y_bus * &v;
        pv.iter()
            .for_each(|&i| s_bus[i].im = (v[i] * i_bus[i]).conj().im);

        izip!(i_bus, s_bus, v)
            .map(|(i_bus, s_bus, v)| i_bus - (s_bus / v).conj())
            .collect_vec()
        // (0..nb)
        //     .map(|i| i_bus[i] - (s_bus[i] / v[i]).conj())
        //     .collect()
    };

    let f: Vec<f64> = [
        pv_pq.iter().map(|&i| mis[i].re).collect::<Vec<f64>>(),
        pv_pq.iter().map(|&i| mis[i].im).collect::<Vec<f64>>(),
    ]
    .concat();

    // check tolerance
    let norm_f = norm_inf(&f);
    if let Some(pm) = progress {
        pm.update(i, norm_f); // max Ir & Ii mismatch (p.u.)
    }
    if norm_f < tol {
        converged = true;
        println!("Converged!"); // TODO: verbose feature
    }

    // do Newton iterations
    while !converged && i < max_it {
        // update iteration counter
        i = i + 1;
    }

    Err("not implemented".to_string())
}

pub(crate) fn newtonpf_S_cart(
    _y_bus: &CSR<usize, Complex64>,
    _s_bus: &dyn SBus,
    _v0: &[Complex64],
    _ref: &[usize],
    _pv: &[usize],
    _pq: &[usize],
    _lin_solver: &dyn LinearSolver,
    _mpopt: &MPOpt,
    _progress: Option<&dyn ProgressMonitor>,
) -> Result<(Vec<Complex64>, bool, usize), String> {
    Err("not implemented".to_string())
}

pub(crate) fn newtonpf_S_hybrid(
    _y_bus: &CSR<usize, Complex64>,
    _s_bus: &dyn SBus,
    _v0: &[Complex64],
    _ref: &[usize],
    _pv: &[usize],
    _pq: &[usize],
    _lin_solver: &dyn LinearSolver,
    _mpopt: &MPOpt,
    _progress: Option<&dyn ProgressMonitor>,
) -> Result<(Vec<Complex64>, bool, usize), String> {
    Err("not implemented".to_string())
}
