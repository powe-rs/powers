use crate::mpc::MPC;
use anyhow::Result;
use casecsv::read::{read_dir, read_zip};
use std::path::PathBuf;

pub fn load_case_file(case_path: &PathBuf) -> Result<MPC> {
    let is_case = match case_path.extension() {
        None => false,
        Some(os_str) => os_str.to_str() == Some("case"),
    };

    let (case, bus, gen, branch, _gencost, _dcline) = if is_case {
        read_zip(case_path)?
    } else {
        read_dir(case_path)?
    };

    let mpc = MPC {
        base_mva: case.base_mva,
        bus,
        gen,
        branch,
        ..Default::default()
    };
    Ok(mpc)
}
