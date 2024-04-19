use caseformat::{Branch, Bus, Gen};
use std::collections::HashMap;

#[derive(Clone, Default, PartialEq)]
pub enum State {
    Internal,
    #[default]
    External,
}

#[derive(Clone)]
pub struct Saved {
    pub(crate) bus: Vec<Bus>,
    pub(crate) gen: Vec<Gen>,
    pub(crate) branch: Vec<Branch>,
}

#[derive(Clone)]
pub struct Order {
    pub state: State,
    pub(crate) internal: Option<Saved>,
    pub(crate) external: Option<Saved>,
    pub bus: BusOrder,
    pub gen: GenOrder,
    pub branch: BranchOrder,
    // pub area: AreaOrder,
}

impl Order {
    pub(crate) fn new(nb: usize, ng: usize, nl: usize) -> Self {
        Self {
            state: State::External,
            internal: None,
            external: None,
            bus: BusOrder::new(nb),
            gen: GenOrder::new(ng),
            branch: BranchOrder::new(nl),
            // area: AreaOrder::default(),
        }
    }
}

#[derive(Clone)]
pub struct Status {
    pub on: Vec<usize>,
    pub off: Vec<usize>,
}

impl Status {
    fn with_capacity(capacity: usize) -> Self {
        Self {
            on: Vec::with_capacity(capacity),
            off: Vec::default(),
        }
    }
}

#[derive(Clone)]
pub struct BusOrder {
    // pub e2i: Vec<usize>,
    pub e2i: HashMap<usize, usize>,
    pub i2e: Vec<usize>,
    pub status: Status,
}

impl BusOrder {
    fn new(nb: usize) -> Self {
        Self {
            e2i: HashMap::with_capacity(nb),
            i2e: vec![0; nb],
            status: Status::with_capacity(nb),
        }
    }
}

#[derive(Clone)]
pub struct GenOrder {
    pub e2i: Vec<usize>,
    pub i2e: Vec<usize>,
    pub status: Status,
}

impl GenOrder {
    fn new(ng: usize) -> Self {
        Self {
            e2i: vec![0; ng],
            i2e: vec![0; ng],
            status: Status::with_capacity(ng),
        }
    }
}

#[derive(Clone)]
pub struct BranchOrder {
    pub status: Status,
}

impl BranchOrder {
    fn new(nl: usize) -> Self {
        Self {
            status: Status::with_capacity(nl),
        }
    }
}

#[derive(Clone)]
pub struct AreaOrder {
    pub status: Status,
}

impl Default for AreaOrder {
    fn default() -> Self {
        Self {
            status: Status::with_capacity(0),
        }
    }
}

#[cfg(test)]
mod tests {
    use anyhow::{format_err, Result};
    use caseformat::NONE;
    use std::env;
    use std::iter::zip;
    use std::path::Path;

    use crate::loadcase::load_case;
    use crate::mpc::MPC;

    use crate::ext_to_int::ext_to_int;
    use crate::int_to_ext::int_to_ext;

    fn setup() -> Result<(MPC, MPC)> {
        let manifest_dir = env::var("CARGO_MANIFEST_DIR")?;
        let casedata_dir = Path::new(&manifest_dir).join("casedata");

        let mpce = load_case(&casedata_dir.join("t_case_ext.case"))?;
        let mut mpci = load_case(&casedata_dir.join("t_case_int.case"))?;

        mpci.bus.iter_mut().for_each(|b| b.bus_i -= 1);
        mpci.gen.iter_mut().for_each(|g| g.gen_bus -= 1);
        mpci.branch.iter_mut().for_each(|br| br.f_bus -= 1);
        mpci.branch.iter_mut().for_each(|br| br.t_bus -= 1);

        Ok((mpce, mpci))
    }

    fn compare_case(expected: &MPC, actual: &MPC) -> Result<()> {
        if let Some((b_exp, b_act)) =
            zip(&expected.bus, &actual.bus).find(|(b_exp, b_act)| b_exp != b_act)
        {
            return Err(format_err!(
                "buses must be equal:\nexpected: {:?}\nactual: {:?}",
                b_exp,
                b_act
            ));
        }
        if let Some((g_exp, g_act)) =
            zip(&expected.gen, &actual.gen).find(|(g_exp, g_act)| g_exp != g_act)
        {
            return Err(format_err!(
                "gens must be equal:\nexpected: {:?}\nactual: {:?}",
                g_exp,
                g_act
            ));
        }
        if let Some((br_exp, br_act)) =
            zip(&expected.branch, &actual.branch).find(|(br_exp, br_act)| br_exp != br_act)
        {
            return Err(format_err!(
                "branches must be equal:\nexpected: {:?}\nactual: {:?}",
                br_exp,
                br_act
            ));
        }
        Ok(())
    }

    #[test]
    fn test_ext_to_int() -> Result<()> {
        let (mpce, mpci) = setup()?;

        let mpc = ext_to_int(&mpce);
        compare_case(&mpci, &mpc)?;

        let mpc = ext_to_int(&mpc);
        compare_case(&mpci, &mpc)?;

        let mpc = int_to_ext(&mpc)?;
        compare_case(&mpce, &mpc)?;

        Ok(())
    }

    #[test]
    fn test_ext_to_int_all_buses_isolated() -> Result<()> {
        let (mut mpc, _) = setup()?;
        mpc.bus.iter_mut().for_each(|b| b.bus_type = NONE);

        let mpci = ext_to_int(&mpc);
        if !mpci.bus.is_empty() {
            return Err(format_err!(
                "internal case must be empty: {}",
                mpci.bus.len()
            ));
        }

        let mpce = int_to_ext(&mpci)?;
        compare_case(&mpc, &mpce)?;

        Ok(())
    }
}
