use crate::{mpc::MPC, mpopt::MPOpt};
use casecsv::{Branch, Bus, Gen};

pub(crate) fn radial_pf(
    base_mva: f64,
    bus0: &[Bus],
    gen0: &[Gen],
    branch0: &[Branch],
    mpopt: &MPOpt,
) -> Result<(MPC, bool, usize), String> {
    Err("not implemented".to_string())
}
