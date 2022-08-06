use densetools::arr::Arr;
use densetools::mat::Mat;
use sparsetools::csc::CSC;

pub trait LinearSolver {
    fn solve(&self, a_mat: CSC<usize, f64>, b: &Arr<f64>) -> Result<Arr<f64>, String>;
}

pub trait LinearSolverN {
    fn solve(&self, a_mat: CSC<usize, f64>, b: &Mat<f64>) -> Result<Mat<f64>, String>;
}
