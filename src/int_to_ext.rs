use crate::order::{Saved, State};
use crate::MPC;
use anyhow::format_err;

pub fn int_to_ext(mpc: &MPC) -> anyhow::Result<MPC> {
    match &mpc.order {
        Some(order) => {
            if order.state == State::External {
                return anyhow::Ok(mpc.clone());
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
        let j = mpc.gen[i].gen_bus;
        mpc.gen[i].gen_bus = order.bus.i2e[j];
    }
    order.external = None;

    anyhow::Ok(mpc)
}
