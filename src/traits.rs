use densetools::arr::Arr;
use sparsetools::csc::CSC;

pub trait LinearSolver {
    fn solve(&self, a_mat: CSC<usize, f64>, b: &Arr<f64>) -> Result<Arr<f64>, String>;
}
