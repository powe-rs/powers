use crate::mpc::{Branch, Bus, BusType, Gen};
use crate::tests::idx;
use densetools::arr::Arr;
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
struct MPCase {
    #[serde(rename = "baseMVA")]
    base_mva: f64,
    bus: Vec<Vec<f64>>,
    gen: Vec<Vec<f64>>,
    branch: Vec<Vec<f64>>,
}

fn array_to_bus(a: &Arr<f64>) -> Bus {
    Bus {
        i: a[idx::BUS_I] as usize,

        bus_type: match a[idx::BUS_TYPE] as usize {
            idx::PQ => BusType::PQ,
            idx::PV => BusType::PV,
            idx::REF => BusType::REF,
            idx::NONE => BusType::NONE,
            _ => BusType::PQ,
        },

        pd: a[idx::PD],
        qd: a[idx::QD],
        gs: a[idx::GS],
        bs: a[idx::BS],

        area: a[idx::BUS_AREA] as usize,

        vm: a[idx::VM],
        va: a[idx::VA],

        base_kv: a[idx::BASE_KV],

        vmax: a[idx::VMAX],
        vmin: a[idx::VMIN],
    }
}

fn array_to_gen(a: &Arr<f64>) -> Gen {
    Gen {
        bus: a[idx::GEN_BUS] as usize,

        pg: a[idx::PG],
        qg: a[idx::QG],

        qmax: a[idx::QMAX],
        qmin: a[idx::QMIN],

        vg: a[idx::VG],

        mbase: a[idx::MBASE],
        status: a[idx::GEN_STATUS] != 0.0,

        pmax: a[idx::PMAX],
        pmin: a[idx::PMIN],
    }
}

fn array_to_branch(a: &Arr<f64>) -> Branch {
    Branch {
        from_bus: a[idx::F_BUS] as usize,
        to_bus: a[idx::T_BUS] as usize,

        r: a[idx::BR_R],
        x: a[idx::BR_X],
        b: a[idx::BR_B],

        rate_a: a[idx::RATE_A],

        tap: a[idx::TAP],
        shift: a[idx::SHIFT],

        status: a[idx::BR_STATUS] != 0.0,

        pf: a[idx::PF],
        qf: a[idx::QF],
        pt: a[idx::PT],
        qt: a[idx::QT],
    }
}
