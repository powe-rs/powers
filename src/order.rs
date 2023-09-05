use crate::{mpc::MPC, mpopt::MPOpt, powers::argsort};

use anyhow::{format_err, Result};
use casecsv::{Branch, Bus, Gen};
use itertools::Itertools;
use std::collections::HashMap;

#[derive(Clone, Default, PartialEq)]
pub enum State {
    Internal,
    #[default]
    External,
}

#[derive(Clone, Default)]
pub struct Order {
    pub state: State,
    pub internal: Option<(Vec<Bus>, Vec<Gen>, Vec<Branch>)>,
    pub external: Option<(Vec<Bus>, Vec<Gen>, Vec<Branch>)>,
    pub bus: BusOrder,
    pub gen: GenOrder,
    pub branch: BranchOrder,
    pub area: AreaOrder,
}

#[derive(Clone, Default)]
pub struct Status {
    pub on: Vec<usize>,
    pub off: Vec<usize>,
}

#[derive(Clone, Default)]
pub struct BusOrder {
    // pub e2i: Vec<usize>,
    pub e2i: HashMap<usize, usize>,
    pub i2e: Vec<usize>,
    pub status: Status,
}

#[derive(Clone, Default)]
pub struct GenOrder {
    pub e2i: Vec<usize>,
    pub i2e: Vec<usize>,
    pub status: Status,
}

#[derive(Clone, Default)]
pub struct BranchOrder {
    pub status: Status,
}

#[derive(Clone, Default)]
pub struct AreaOrder {
    pub status: Status,
}

pub fn ext2int(mpc: &MPC, mpopt: &MPOpt, reorder_gens: bool) -> MPC {
    let mut mpc: MPC = mpc.clone();
    if let Some(order) = mpc.order.take().as_ref() {
        if order.state == State::External {
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
        None => Order::default(),
    };

    // save data with external ordering
    o.external = Some((mpc.bus.clone(), mpc.gen.clone(), mpc.branch.clone()));
    // o.external.bus = mpc.bus.clone();
    // o.external.branch = mpc.branch.clone();
    // o.external.gen = mpc.gen.clone();
    // if let Some(areas) = mpc.areas.as_ref() {
    //     if areas.is_empty() {
    //         mpc.areas = None; // if areas field is empty delete it (so it gets ignored)
    //     } else {
    //         o.external.areas = Some(areas.clone()); // otherwise save it
    //     }
    // }

    // determine which buses, branches, gens are connected & in-service
    let mut bs = HashMap::new();
    for (i, b) in mpc.bus.iter().enumerate() {
        if b.bus_type != 4 {
            o.bus.status.on[i] = 1;
        } else {
            o.bus.status.off[i] = 1;
        }
        bs.insert(b.bus_i, b.bus_type != 4);
    }
    for (i, g) in mpc.gen.iter().enumerate() {
        if g.is_on() && bs[&g.bus] {
            o.gen.status.on[i] = 1;
        } else {
            o.gen.status.off[i] = 1;
        }
    }
    for (i, br) in mpc.branch.iter().enumerate() {
        if br.is_on() && bs[&br.f_bus] && bs[&br.t_bus] {
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
        o.bus.i2e[i] = b.bus_i;
        o.bus.e2i.insert(b.bus_i, i);
    }
    if nb != 0 {
        for b in mpc.bus.iter_mut() {
            b.bus_i = o.bus.e2i[&b.bus_i];
        }
        for g in mpc.gen.iter_mut() {
            g.bus = o.bus.e2i[&g.bus];
        }
        for br in mpc.branch.iter_mut() {
            br.f_bus = o.bus.e2i[&br.f_bus];
            br.t_bus = o.bus.e2i[&br.t_bus];
        }
    }

    if reorder_gens {
        // reorder gens in order of increasing bus number
        o.gen.i2e = argsort(&mut mpc.gen.iter().map(|g| g.bus).collect_vec(), false);
        o.gen.e2i = argsort(&mut o.gen.i2e.clone(), false);
        todo!("reorder_gens");
    } else {
        // don't reorder gens in order of increasing bus number, but
        // keep the mappings in place for backward compatibility
        o.gen.i2e = (0..ng).collect();
        o.gen.e2i = o.gen.i2e.clone();
    }

    if o.internal.is_some() {
        o.internal = None;
    }
    o.state = State::Internal;
    mpc.order = Some(o);

    mpc
}

pub(crate) fn int2ext(mpc: &MPC, mpopt: Option<&MPOpt>) -> Result<MPC> {
    let mut mpc = mpc.clone();

    if mpc.order.is_none() {
        return Err(format_err!(
            "mpc must have \"order\" field for conversion back to external numbering"
        ));
    }
    let mut o = mpc.order.as_mut().unwrap();

    if o.state == State::Internal {
        // TODO: execute userfcn callbacks for 'int2ext' stage

        // TODO: convert back "extra" fields

        // Save data matrices with internal ordering & restore originals.
        if let Some(internal) = o.internal.as_mut() {
            internal.0 = mpc.bus.clone();
            internal.2 = mpc.branch.clone();
            internal.1 = mpc.gen.clone();
        } else {
        }
        if let Some(external) = &o.external {
            mpc.bus = external.0.clone();
            mpc.branch = external.2.clone();
            mpc.gen = external.1.clone();
        } else {
        }

        // TODO: zero pad data matrices on right if necessary

        // Update data (in bus, branch, and gen only).
        for (i, &j) in o.bus.status.on.iter().enumerate() {
            mpc.bus[j] = o.internal.as_ref().unwrap().0[i].clone();
        }
        for (i, &j) in o.branch.status.on.iter().enumerate() {
            mpc.branch[j] = o.internal.as_ref().unwrap().2[i].clone();
        }
        for (i, &j) in o.gen.status.on.iter().enumerate() {
            mpc.gen[j] = o.internal.as_ref().unwrap().1[o.gen.e2i[i]].clone();
        }

        // revert to original bus numbers
        for &i in &o.bus.status.on {
            let j = mpc.bus[i].bus_i;
            mpc.bus[i].bus_i = o.bus.i2e[j];
        }
        for &i in &o.branch.status.on {
            let j = mpc.branch[i].f_bus;
            mpc.branch[i].f_bus = o.bus.i2e[j];
        }
        for &i in &o.branch.status.on {
            let j = mpc.branch[i].t_bus;
            mpc.branch[i].t_bus = o.bus.i2e[j];
        }
        for &i in &o.gen.status.on {
            let j = mpc.gen[i].bus;
            mpc.gen[i].bus = o.bus.i2e[j];
        }
        o.external = None;
    }

    Ok(mpc)
}
