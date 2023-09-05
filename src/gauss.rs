use crate::mpopt::MPOpt;
use crate::powers::SBus;
use crate::traits::LinearSolver;
use num_complex::Complex64;
use sparsetools::csr::CSR;

pub trait ProgressMonitor {
    fn update(&self, i: usize, norm_f: f64);
}

pub struct PrintProgress {}

impl ProgressMonitor for PrintProgress {
    fn update(&self, i: usize, norm_f: f64) {
        if i == 0 {
            println!(" it    max P & Q mismatch (p.u.)");
            println!("----  ---------------------------");
            // println!("%3d        %10.3e", i, normF);
            println!("{}        {}", i, norm_f);
        } else {
            // println!("%3d        %10.3e", i, normF);
            println!("{}        {}", i, norm_f);
        }
    }
}

pub(crate) fn gausspf(
    y_bus: &CSR<usize, Complex64>,
    s_bus: &dyn SBus,
    v0: &[Complex64],
    _ref: &[usize],
    pv: &[usize],
    pq: &[usize],
    lin_solver: &dyn LinearSolver,
    mpopt: &MPOpt,
    progress: Option<&dyn ProgressMonitor>,
) -> Result<(Vec<Complex64>, bool, usize), String> {
    Err("not implemented".to_string())
}
