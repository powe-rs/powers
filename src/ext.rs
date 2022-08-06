use crate::mpc::{BusType, MPC};
use crate::mpopt::MPOpt;
use densetools::arr::Arr;
use densetools::slice::argsort;
use std::collections::HashMap;

#[derive(Clone)]
pub struct Order {
    pub state: String, // 'i' or 'e'
    pub int: Option<MPC>,
    pub external: MPC,
    pub bus: BusOrder,
    pub gen: GenOrder,
    pub branch: BranchAreas,
    pub area: BranchAreas,
}

impl Order {
    pub(crate) fn new(nb: usize) -> Self {
        Self {
            state: "e".to_string(),
            int: None,
            external: MPC::new(),
            bus: BusOrder::new(),
            gen: GenOrder::new(),
            branch: BranchAreas::new(),
            area: BranchAreas::new(),
        }
    }

    pub(crate) fn is_external(&self) -> bool {
        self.state == "e"
    }

    pub(crate) fn is_internal(&self) -> bool {
        self.state == "i"
    }
}

#[derive(Clone)]
pub struct Status {
    pub on: Vec<usize>,
    pub off: Vec<usize>,
}

impl Status {
    fn new() -> Self {
        Self {
            on: vec![],
            off: vec![],
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
    fn new() -> Self {
        Self {
            e2i: HashMap::new(),
            i2e: vec![],
            status: Status::new(),
        }
    }
}

#[derive(Clone)]
pub struct GenOrder {
    pub e2i: Vec<usize>,
    // pub e2i: HashMap<usize, usize>,
    pub i2e: Vec<usize>,
    pub status: Status,
}

impl GenOrder {
    fn new() -> Self {
        Self {
            e2i: vec![],
            i2e: vec![],
            status: Status::new(),
        }
    }
}

#[derive(Clone)]
pub struct BranchAreas {
    pub status: Status,
}

impl BranchAreas {
    fn new() -> Self {
        Self {
            status: Status::new(),
        }
    }
}

pub fn ext2int(mpc: &MPC, mpopt: &MPOpt, reorder_gens: bool) -> MPC {
    let mut mpc: MPC = mpc.clone();
    if let Some(order) = mpc.order.take().as_ref() {
        if order.is_external() {
            return mpc.clone();
        }
    }

    // sizes
    let nb = mpc.bus.len();
    let ng = mpc.gen.len();
    let ng0 = ng;
    let dc = if mpc.a_mat.is_some() && mpc.a_mat.as_ref().unwrap().cols() < 2 * nb + 2 * ng {
        true
    } else if mpc.n_mat.is_some() && mpc.n_mat.as_ref().unwrap().cols() < 2 * nb + 2 * ng {
        true
    } else {
        false
    };

    // initialize order
    // let o = if first { Order::new() } else { mpc.order.unwrap() };
    // let mut o = mpc.order.unwrap_or(Order::new(nb));
    let mut o = match mpc.order.take().as_ref() {
        Some(order) => order.clone(),
        None => Order::new(nb),
    };

    // save data with external ordering
    o.external.bus = mpc.bus.clone();
    o.external.branch = mpc.branch.clone();
    o.external.gen = mpc.gen.clone();
    if let Some(areas) = mpc.areas.as_ref() {
        if areas.is_empty() {
            mpc.areas = None; // if areas field is empty delete it (so it gets ignored)
        } else {
            o.external.areas = Some(areas.clone()); // otherwise save it
        }
    }

    // determine which buses, branches, gens are connected & in-service
    let mut bs = HashMap::new();
    for (i, b) in mpc.bus.iter().enumerate() {
        if b.is_on() {
            o.bus.status.on[i] = 1;
        } else {
            o.bus.status.off[i] = 1;
        }
        bs.insert(b.i, b.is_on());
    }
    for (i, g) in mpc.gen.iter().enumerate() {
        if g.status && bs[&g.bus] {
            o.gen.status.on[i] = 1;
        } else {
            o.gen.status.off[i] = 1;
        }
    }
    for (i, br) in mpc.branch.iter().enumerate() {
        if br.status && bs[&br.from_bus] && bs[&br.to_bus] {
            o.branch.status.on[i] = 1;
        } else {
            o.branch.status.off[i] = 1;
        }
    }
    let on_bus = mpc
        .bus
        .iter()
        .enumerate()
        .filter(|(i, _)| o.bus.status.on[*i] != 0)
        .map(|(_, b)| b.clone())
        .collect();
    mpc.bus = on_bus;

    // update sizes
    let nb = mpc.bus.len();
    let ng = mpc.gen.len();

    // apply consecutive bus numbering
    // o.bus.i2e = mpc.bus.iter().map(|b| b.i).collect();
    for (i, b) in mpc.bus.iter().enumerate() {
        o.bus.i2e[i] = b.i;
        o.bus.e2i.insert(b.i, i);
    }
    if nb != 0 {
        for b in mpc.bus.iter_mut() {
            b.i = o.bus.e2i[&b.i];
        }
        for g in mpc.gen.iter_mut() {
            g.bus = o.bus.e2i[&g.bus];
        }
        for br in mpc.branch.iter_mut() {
            br.from_bus = o.bus.e2i[&br.from_bus];
            br.to_bus = o.bus.e2i[&br.to_bus];
        }
    }

    if reorder_gens {
        // reorder gens in order of increasing bus number
        o.gen.i2e = argsort(
            &mut mpc.gen.iter().map(|g| g.bus).collect::<Vec<usize>>(),
            false,
        );
        o.gen.e2i = argsort(&mut o.gen.i2e.clone(), false);
    } else {
        // don't reorder gens in order of increasing bus number, but
        // keep the mappings in place for backward compatibility
        o.gen.i2e = (0..ng).collect();
        o.gen.e2i = o.gen.i2e.clone();
    }

    if o.int.is_some() {
        o.int = None;
    }
    o.state = "i".to_string();
    mpc.order = Box::new(Some(o));

    mpc
}

pub fn int2ext(mpc: &MPC, mpopt: Option<&MPOpt>) -> MPC {
    let mpc: MPC = mpc.clone();

    mpc
}
