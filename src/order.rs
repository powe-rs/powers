use crate::mpc::MPC;

use anyhow::{format_err, Ok, Result};
use casecsv::{Branch, Bus, Gen};
use std::collections::{HashMap, HashSet};

#[derive(Clone, Default, PartialEq)]
pub enum State {
    Internal,
    #[default]
    External,
}

#[derive(Clone)]
pub struct Saved {
    bus: Vec<Bus>,
    gen: Vec<Gen>,
    branch: Vec<Branch>,
}

#[derive(Clone)]
pub struct Order {
    pub state: State,
    pub internal: Option<Saved>,
    pub external: Option<Saved>,
    pub bus: BusOrder,
    pub gen: GenOrder,
    pub branch: BranchOrder,
    pub area: AreaOrder,
}

impl Order {
    fn new(nb: usize, ng: usize, nl: usize) -> Self {
        Self {
            state: State::External,
            internal: None,
            external: None,
            bus: BusOrder::new(nb),
            gen: GenOrder::new(ng),
            branch: BranchOrder::new(nl),
            area: AreaOrder::default(),
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

pub fn ext_to_int(mpc: &MPC) -> MPC {
    if let Some(order) = mpc.order.as_ref() {
        if order.state == State::Internal {
            return mpc.clone();
        }
    }

    // sizes
    let nb = mpc.bus.len();
    let ng = mpc.gen.len();
    let nl = mpc.branch.len();
    // let ng0 = ng;
    // let dc = if mpc.a_mat.is_some() && mpc.a_mat.as_ref().unwrap().cols() < 2 * nb + 2 * ng {
    //     true
    // } else if mpc.n_mat.is_some() && mpc.n_mat.as_ref().unwrap().cols() < 2 * nb + 2 * ng {
    //     true
    // } else {
    //     false
    // };

    // initialize order
    // let mut order = mpc.order.take().unwrap_or(Order::new(nb, ng, nl));
    let mut order = Order::new(nb, ng, nl);

    // save data with external ordering
    order.external = Some(Saved {
        bus: mpc.bus.clone(),
        gen: mpc.gen.clone(),
        branch: mpc.branch.clone(),
    });
    // if let Some(areas) = mpc.areas.as_ref() {
    //     if areas.is_empty() {
    //         mpc.areas = None; // if areas field is empty delete it (so it gets ignored)
    //     } else {
    //         o.external.areas = Some(areas.clone()); // otherwise save it
    //     }
    // }

    // determine which buses, branches, gens are connected & in-service
    let bs: HashSet<usize> = mpc
        .bus
        .iter()
        .filter(|b| b.is_pq() || b.is_pv() || b.is_ref())
        .map(|b| b.bus_i)
        .collect();
    for (i, b) in mpc.bus.iter().enumerate() {
        if b.is_pq() || b.is_pv() || b.is_ref() {
            order.bus.status.on.push(i);
        } else {
            order.bus.status.off.push(i);
        }
    }
    for (i, g) in mpc.gen.iter().enumerate() {
        if g.is_on() && bs.contains(&g.bus) {
            order.gen.status.on.push(i);
        } else {
            order.gen.status.off.push(i);
        }
    }
    for (i, br) in mpc.branch.iter().enumerate() {
        if br.is_on() && bs.contains(&br.f_bus) && bs.contains(&br.t_bus) {
            order.branch.status.on.push(i);
        } else {
            order.branch.status.off.push(i);
        }
    }

    let mut mpc = mpc.clone();

    // delete stuff that is "out"
    mpc.bus.retain(|b| b.is_pq() || b.is_pv() || b.is_ref());
    mpc.branch
        .retain(|br| br.is_on() && bs.contains(&br.f_bus) && bs.contains(&br.t_bus));
    mpc.gen.retain(|g| g.is_on() && bs.contains(&g.bus));

    // update sizes
    let nb = mpc.bus.len();
    let ng = mpc.gen.len();

    // apply consecutive bus numbering
    for (i, b) in mpc.bus.iter().enumerate() {
        order.bus.i2e[i] = b.bus_i;
        order.bus.e2i.insert(b.bus_i, i);
    }
    if nb != 0 {
        for b in mpc.bus.iter_mut() {
            b.bus_i = order.bus.e2i[&b.bus_i];
        }
        for g in mpc.gen.iter_mut() {
            g.bus = order.bus.e2i[&g.bus];
        }
        for br in mpc.branch.iter_mut() {
            br.f_bus = order.bus.e2i[&br.f_bus];
            br.t_bus = order.bus.e2i[&br.t_bus];
        }
    }

    // if reorder_gens {
    //     // reorder gens in order of increasing bus number
    //     order.gen.i2e = argsort(&mut mpc.gen.iter().map(|g| g.bus).collect_vec(), false);
    //     order.gen.e2i = argsort(&mut order.gen.i2e.clone(), false);
    //     todo!("reorder_gens");
    // } else {
    // don't reorder gens in order of increasing bus number, but
    // keep the mappings in place for backward compatibility
    order.gen.i2e = (0..ng).collect();
    order.gen.e2i = order.gen.i2e.clone();
    // }

    order.internal = None;
    order.state = State::Internal;
    mpc.order = Some(order);

    mpc
}

pub(crate) fn int_to_ext(mpc: &MPC) -> Result<MPC> {
    match &mpc.order {
        Some(order) => {
            if order.state == State::External {
                return Ok(mpc.clone());
            }
        }
        None => {
            return Err(format_err!(
                "mpc must have \"order\" field for conversion back to external numbering"
            ));
        }
    }
    let mut mpc = mpc.clone();

    // TODO: execute userfcn callbacks for 'int2ext' stage

    // TODO: convert back "extra" fields

    let order = mpc.order.as_mut().unwrap();

    // Save data matrices with internal ordering & restore originals.
    let internal = Saved {
        bus: mpc.bus.drain(..).collect(),
        branch: mpc.branch.drain(..).collect(),
        gen: mpc.gen.drain(..).collect(),
    };
    if let Some(external) = order.external.as_mut() {
        mpc.bus = external.bus.drain(..).collect();
        mpc.branch = external.branch.drain(..).collect();
        mpc.gen = external.gen.drain(..).collect();
    }

    // TODO: zero pad data matrices on right if necessary

    // Update data (in bus, branch, and gen only).
    for (i, &j) in order.bus.status.on.iter().enumerate() {
        mpc.bus[j] = internal.bus[i].clone();
    }
    for (i, &j) in order.branch.status.on.iter().enumerate() {
        mpc.branch[j] = internal.branch[i].clone();
    }
    for (i, &j) in order.gen.status.on.iter().enumerate() {
        mpc.gen[j] = internal.gen[order.gen.e2i[i]].clone();
    }
    order.internal = Some(internal);

    // revert to original bus numbers
    for &i in &order.bus.status.on {
        let j = mpc.bus[i].bus_i;
        mpc.bus[i].bus_i = order.bus.i2e[j];
    }
    for &i in &order.branch.status.on {
        let j = mpc.branch[i].f_bus;
        mpc.branch[i].f_bus = order.bus.i2e[j];
    }
    for &i in &order.branch.status.on {
        let j = mpc.branch[i].t_bus;
        mpc.branch[i].t_bus = order.bus.i2e[j];
    }
    for &i in &order.gen.status.on {
        let j = mpc.gen[i].bus;
        mpc.gen[i].bus = order.bus.i2e[j];
    }
    order.external = None;

    Ok(mpc)
}

#[cfg(test)]
mod tests {
    use anyhow::{format_err, Result};
    use casecsv::NONE;
    use itertools::izip;
    use std::env;
    use std::path::Path;

    use crate::loadcase::load_case;
    use crate::mpc::MPC;

    use super::{ext_to_int, int_to_ext};

    fn setup() -> Result<(MPC, MPC)> {
        let manifest_dir = env::var("CARGO_MANIFEST_DIR")?;
        let casedata_dir = Path::new(&manifest_dir).join("casedata");

        let mpce = load_case(&casedata_dir.join("t_case_ext.case"))?;
        let mut mpci = load_case(&casedata_dir.join("t_case_int.case"))?;

        mpci.bus.iter_mut().for_each(|b| b.bus_i -= 1);
        mpci.gen.iter_mut().for_each(|g| g.bus -= 1);
        mpci.branch.iter_mut().for_each(|br| br.f_bus -= 1);
        mpci.branch.iter_mut().for_each(|br| br.t_bus -= 1);

        Ok((mpce, mpci))
    }

    fn compare_case(expected: &MPC, actual: &MPC) -> Result<()> {
        if let Some((b_exp, b_act)) =
            izip!(&expected.bus, &actual.bus).find(|(b_exp, b_act)| b_exp != b_act)
        {
            return Err(format_err!(
                "buses must be equal:\nexpected: {:?}\nactual: {:?}",
                b_exp,
                b_act
            ));
        }
        if let Some((g_exp, g_act)) =
            izip!(&expected.gen, &actual.gen).find(|(g_exp, g_act)| g_exp != g_act)
        {
            return Err(format_err!(
                "gens must be equal:\nexpected: {:?}\nactual: {:?}",
                g_exp,
                g_act
            ));
        }
        if let Some((br_exp, br_act)) =
            izip!(&expected.branch, &actual.branch).find(|(br_exp, br_act)| br_exp != br_act)
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
