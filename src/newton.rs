use crate::math::{csr_select, norm_inf};
use crate::mpopt::MPOpt;
use crate::rpower::d_sbus_d_v;
use crate::traits::LinearSolver;
use ndarray::Array1;
use num_complex::Complex64;
use sprs::{hstack, vstack, CsMat, CsMatView};

pub trait SBus {
    fn s_bus(&self, v_m: &[f64]) -> (Array1<Complex64>, Option<CsMat<Complex64>>) {
        panic!()
    }
}

pub trait ProgressMonitor {
    fn update(&self, i: usize, norm_f: f64);
}

/// Solves power flow using full Newton's method (power/polar).
pub(crate) fn newtonpf(
    y_bus: CsMatView<Complex64>,
    s_bus: &dyn SBus,
    v0: Array1<Complex64>,
    _ref: &[usize],
    pv: &[usize],
    pq: &[usize],
    lin_solver: &dyn LinearSolver,
    mpopt: &MPOpt,
    progress: Option<&dyn ProgressMonitor>,
) -> Result<(Array1<Complex64>, bool, usize), String> {
    let tol = mpopt.pf.tolerance;
    let max_it = mpopt.pf.max_it_nr;
    // let lin_solver  = mpopt.pf.nr.lin_solver;

    let mut converged = false;
    let mut i = 0;
    let mut v = v0.clone();
    let mut va = v.iter().map(|v_i| v_i.arg()).collect::<Vec<f64>>();
    let mut vm = v.iter().map(|v_i| v_i.norm()).collect::<Vec<f64>>();

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
    let i_bus_conj: Array1<Complex64> = (&y_bus * &v.view()).iter().map(|i| i.conj()).collect();
    let mis: Array1<Complex64> = &v * i_bus_conj - s_bus.s_bus(&vm).0;
    let mut f = Array1::from(
        [
            pv.iter().map(|&pv_i| mis[pv_i].re).collect::<Vec<f64>>(),
            pq.iter().map(|&pq_i| mis[pq_i].re).collect::<Vec<f64>>(),
            pq.iter().map(|&pq_i| mis[pq_i].im).collect::<Vec<f64>>(),
        ]
        .concat(),
    );

    // check tolerance
    let norm_f = norm_inf(f.view());
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
        let (d_sbus_d_va, mut d_sbus_d_vm) = d_sbus_d_v(y_bus, v.as_slice().unwrap(), false);
        let (_, neg_d_sd_d_vm) = s_bus.s_bus(&vm);
        if let Some(neg_d_sd_d_vm) = neg_d_sd_d_vm {
            d_sbus_d_vm = &d_sbus_d_vm - &neg_d_sd_d_vm;
        }

        let j11 = csr_select(
            d_sbus_d_va.view(),
            &[pv, pq].concat(),
            &[pv, pq].concat(),
            true,
        );
        let j12 = csr_select(d_sbus_d_vm.view(), &[pv, pq].concat(), pq, true);
        let j21 = csr_select(d_sbus_d_va.view(), pq, &[pv, pq].concat(), false);
        let j22 = csr_select(d_sbus_d_vm.view(), pq, pq, false);
        // j11 = real(dSbus_dVa([pv; pq], [pv; pq]));
        // j12 = real(dSbus_dVm([pv; pq], pq));
        // j21 = imag(dSbus_dVa(pq, [pv; pq]));
        // j22 = imag(dSbus_dVm(pq, pq));

        let jac = vstack(&[
            hstack(&[j11.view(), j12.view()]).view(),
            hstack(&[j21.view(), j22.view()]).view(),
        ]);

        // compute update step
        // let dx = lin_solver.solve(jac.view(), f.iter().map(|&f_i| -f_i).collect())?;
        let dx = lin_solver.solve(jac.view(), (-1.0 * f).view())?;

        // update voltage
        for (i, j) in (j1..j2).enumerate() {
            va[pv[i]] = va[pv[i]] + dx[j];
        }
        for (i, j) in (j3..j4).enumerate() {
            va[pq[i]] = va[pq[i]] + dx[j];
        }
        for (i, j) in (j5..j6).enumerate() {
            vm[pq[i]] = vm[pq[i]] + dx[j];
        }
        // update Vm and Va again in case we wrapped around with a negative Vm
        v = vm
            .iter()
            .zip(&va)
            .map(|(&vm, &va)| Complex64::from_polar(vm, va))
            .collect();
        va = v.iter().map(|v| v.arg()).collect();
        vm = v.iter().map(|v| v.norm()).collect();

        // evalute F(x)
        let i_bus_conj: Array1<Complex64> = (&y_bus * &v.view()).iter().map(|i| i.conj()).collect();
        let mis: Array1<Complex64> = &v * i_bus_conj - s_bus.s_bus(&vm).0;
        f = Array1::from(
            [
                pv.iter().map(|&pv_i| mis[pv_i].re).collect::<Vec<f64>>(),
                pq.iter().map(|&pq_i| mis[pq_i].re).collect::<Vec<f64>>(),
                pq.iter().map(|&pq_i| mis[pq_i].im).collect::<Vec<f64>>(),
            ]
            .concat(),
        );

        // check for convergence
        let norm_f = norm_inf(f.view());
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
