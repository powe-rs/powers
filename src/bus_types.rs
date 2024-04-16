use std::collections::HashSet;
use caseformat::{Bus, Gen};

/// Builds index lists for each type of bus (REF, PV, PQ).
///
/// Generators with "out-of-service" status are treated as PQ buses with
/// zero generation (regardless of Pg/Qg values in gen). Expects `BUS` and
/// `GEN` have been converted to use internal consecutive bus numbering.
pub fn bus_types(bus: &[Bus], gen: &[Gen]) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
    // Buses with generators that are ON.
    // let mut bus_gen_status = vec![false; bus.len()];
    // for g in gen {
    //     if g.is_on() {
    //         bus_gen_status[g.bus] = true
    //     }
    // }
    let bus_gen_status = gen
        .iter()
        .filter(|g| g.is_on())
        .map(|g| g.gen_bus)
        .collect::<HashSet<usize>>();

    // Form index lists for slack, PV, and PQ buses.
    // ref = find(bus(:, BUS_TYPE) == REF & bus_gen_status);   /// reference bus index
    // pv  = find(bus(:, BUS_TYPE) == PV  & bus_gen_status);   /// PV bus indices
    // pq  = find(bus(:, BUS_TYPE) == PQ | ~bus_gen_status);   /// PQ bus indices
    let refbus = bus
        .iter()
        .filter(|b| b.is_ref() && bus_gen_status.contains(&b.bus_i))
        .map(|b| b.bus_i)
        .collect::<Vec<usize>>();
    let pv = bus
        .iter()
        .filter(|b| b.is_pv() && bus_gen_status.contains(&b.bus_i))
        .map(|b| b.bus_i)
        .collect::<Vec<usize>>();
    let pq = bus
        .iter()
        .filter(|b| b.is_pq() || !bus_gen_status.contains(&b.bus_i))
        .map(|b| b.bus_i)
        .collect::<Vec<usize>>();

    (refbus, pv, pq)
}
