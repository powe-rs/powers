use crate::mpc::MPC;
use anyhow::Result;
use caseformat::{read_dir, read_zip};
use std::fs::File;
use std::path::PathBuf;
use crate::ext_to_int;

pub fn load_case(case_path: &PathBuf) -> Result<MPC> {
    let is_case = match case_path.extension() {
        None => false,
        Some(os_str) => match os_str.to_str() {
            Some("case") | Some("zip") => true,
            _ => false,
        },
    };

    let (case, bus, gen, branch, _gencost, _dcline, _readme, _license) = if is_case {
        read_zip(File::open(case_path)?)?
    } else {
        read_dir(case_path)?
    };

    let mpc = MPC {
        name: case.name,
        base_mva: case.base_mva,
        bus,
        gen,
        branch,
        ..Default::default()
    };

    Ok(ext_to_int(&mpc))
}
