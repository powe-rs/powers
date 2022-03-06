use num_complex::Complex64;
use sprs::CsMat;

pub(crate) trait Conj {
    fn conj(&self) -> Self;
}

impl Conj for CsMat<Complex64> {
    fn conj(&self) -> Self {
        self.map(|a| a.conj())
    }
}
