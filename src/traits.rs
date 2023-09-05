use sparsetools::csc::CSC;

pub trait LinearSolver {
    fn solve(&self, a_mat: CSC<usize, f64>, b: &[f64]) -> Result<Vec<f64>, String>;
}

pub trait LinearSolverN {
    fn solve(&self, a_mat: CSC<usize, f64>, b: &[f64], nrhs: usize) -> Result<Vec<f64>, String>;
}
