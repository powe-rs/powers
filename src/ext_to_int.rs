use std::collections::HashSet;
use crate::{MPC};
use crate::order::{Order, Saved, State};

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
        if g.is_on() && bs.contains(&g.gen_bus) {
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
    mpc.gen.retain(|g| g.is_on() && bs.contains(&g.gen_bus));

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
            g.gen_bus = order.bus.e2i[&g.gen_bus];
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
