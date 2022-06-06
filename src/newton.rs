use crate::mpopt::MPOpt;
use crate::rpower::d_sbus_d_v;
use crate::traits::LinearSolver;
use densetools::arr::{Arr, CArr};
use densetools::slice::norm_inf;
use num_complex::Complex64;
use sparsetools::coo::Coo;
use sparsetools::csr::{CCSR, CSR};

pub trait SBus {
    fn s_bus(&self, v_m: &[f64]) -> (Arr<Complex64>, Option<CSR<usize, Complex64>>) {
        panic!()
    }
}

pub trait ProgressMonitor {
    fn update(&self, i: usize, norm_f: f64);
}

/// Solves power flow using full Newton's method (power/polar).
pub(crate) fn newtonpf(
    y_bus: &CSR<usize, Complex64>,
    s_bus: &dyn SBus,
    v0: Arr<Complex64>,
    _ref: &[usize],
    pv: &[usize],
    pq: &[usize],
    lin_solver: &dyn LinearSolver,
    mpopt: &MPOpt,
    progress: Option<&dyn ProgressMonitor>,
) -> Result<(Arr<Complex64>, bool, usize), String> {
    let pv_pq = [pv, pq].concat();

    let tol = mpopt.pf.tolerance;
    let max_it = mpopt.pf.max_it_nr;
    // let lin_solver  = mpopt.pf.nr.lin_solver;

    let mut converged = false;
    let mut i = 0;
    let mut v = v0.clone();
    let mut va = v.arg();
    let mut vm = v.norm();

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
    let mis = &v * (y_bus * &v).conj() - s_bus.s_bus(&vm).0;
    let mut f = Arr::concat(&mis.select(&pv_pq).real(), &mis.select(pq).imag());

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
        let (_, neg_d_sd_d_vm) = s_bus.s_bus(&vm);
        d_sbus_d_vm = d_sbus_d_vm - neg_d_sd_d_vm.unwrap();

        let j11 = d_sbus_d_va.select(Some(&pv_pq), Some(&pv_pq))?.real();
        let j12 = d_sbus_d_vm.select(Some(&pv_pq), Some(pq))?.real();
        let j21 = d_sbus_d_va.select(Some(pq), Some(&pv_pq))?.imag();
        let j22 = d_sbus_d_vm.select(Some(pq), Some(pq))?.imag();

        let jac = Coo::compose(&j11.to_coo(), &j12.to_coo(), &j21.to_coo(), &j22.to_coo())?;

        // compute update step
        // let dx = lin_solver.solve(jac.view(), f.iter().map(|&f_i| -f_i).collect())?;
        let dx = lin_solver.solve(jac.to_csc(), &-f)?;

        // update voltage
        pv.iter().zip(j1..j2).for_each(|(&i, j)| va[i] += dx[j]);
        pq.iter().zip(j3..j4).for_each(|(&i, j)| va[i] += dx[j]);
        pq.iter().zip(j5..j6).for_each(|(&i, j)| vm[i] += dx[j]);
        // for (i, j) in (j1..j2).enumerate() {
        //     va[pv[i]] = va[pv[i]] + dx[j];
        // }
        // for (i, j) in (j3..j4).enumerate() {
        //     va[pq[i]] = va[pq[i]] + dx[j];
        // }
        // for (i, j) in (j5..j6).enumerate() {
        //     vm[pq[i]] = vm[pq[i]] + dx[j];
        // }

        // update Vm and Va again in case we wrapped around with a negative Vm
        v = Arr::from_polar(&vm, &va);
        va = v.arg();
        vm = v.norm();

        // evalute F(x)
        let mis = &v * (y_bus * &v).conj() - s_bus.s_bus(&vm).0;
        f = Arr::concat(&mis.select(&pv_pq).real(), &mis.select(pq).imag());

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
