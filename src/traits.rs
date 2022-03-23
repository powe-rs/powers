use ndarray::{Array1, ArrayView1};
use num_complex::Complex64;
use sprs::{CsMat, CsMatView};

pub(crate) trait Conj {
    fn conj(&self) -> Self;
}

impl Conj for CsMat<Complex64> {
    fn conj(&self) -> Self {
        self.map(|a| a.conj())
    }
}

pub trait LinearSolver {
    fn solve(&self, a_mat: CsMatView<f64>, b: ArrayView1<f64>) -> Result<Array1<f64>, String>;
}
