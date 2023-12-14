mod current;
mod power;

pub(crate) use current::*;
pub(crate) use power::*;

pub trait ProgressMonitor {
    fn update(&self, i: usize, norm_f: f64);
}
