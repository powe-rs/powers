use crate::cmplx;
use crate::mpc::{Bus, BusType, Gen};
use crate::mpopt::MPOpt;
use crate::rpower::make_sdzip;
use densetools::arr::Arr;

#[derive(PartialEq)]
pub enum LoadZone {
    /// Load zones defined by the areas specified in the bus data.
    Area = 0,
    /// Use a single zone for the entire system.
    All = 1,
    /// Use a different zone for each bus.
    Bus = 2,
}

#[derive(PartialEq)]
pub enum LoadType {
    /// Sum both fixed and dispatchable loads.
    Both = 0,
    /// Sum only fixed loads.
    Fixed = 1,
    /// Sum only dispatchable loads.
    Dispatchable = 2,
}

/// Returns vector of total load in each load zone.
///
/// `load_zone` - (optional) nb element vector where the value of
/// each element is either zero or the index of the load zone
/// to which the corresponding bus belongs. If LOAD_ZONE(b) = k
/// then the loads at bus b will added to the values of PD(k) and
/// QD(k). If LOAD_ZONE is empty, the default is defined as the areas
/// specified in the BUS matrix, i.e. LOAD_ZONE = BUS(:, BUS_AREA)
/// and load will have dimension = MAX(BUS(:, BUS_AREA)). LOAD_ZONE
/// can also take the following string values:
///     'all'  - use a single zone for the entire system (return scalar)
///     'area' - use LOAD_ZONE = BUS(:, BUS_AREA), same as default
///     'bus'  - use a different zone for each bus (i.e. to compute
///         final values of bus-wise loads, including voltage dependent
///         fixed loads and or dispatchable loads)
///
/// `load_type`  -  enum specifying types of loads to include, default
/// is `BOTH` if GEN is provided, otherwise `FIXED`
///
///         'FIXED'        : sum only fixed loads
///         'DISPATCHABLE' : sum only dispatchable loads
///         'BOTH'         : sum both fixed and dispatchable loads
///
/// `nominal` -  use nominal load for dispatchable loads. Otherwise, use
/// actual realized load for dispatchable loads
pub(crate) fn total_load(
    bus: &[Bus],
    gen: Option<&[Gen]>,
    load_zone: LoadZone,
    load_zones: Option<&[Option<usize>]>,
    load_type: LoadType,
    nominal: bool,
    mpopt: &MPOpt,
    want_q: bool,
) -> (Vec<f64>, Option<Vec<f64>>) {
    let nb = bus.len();

    // let want_fixed = load_type == LoadType::Both || load_type == LoadType::Fixed;
    let want_fixed = match load_type {
        LoadType::Both | LoadType::Fixed => true,
        _ => false,
    };
    // let want_disp = load_type == LoadType::Both || load_type == LoadType::Dispatchable;
    let want_disp = match load_type {
        LoadType::Both | LoadType::Dispatchable => true,
        _ => false,
    };

    let lz = match load_zone {
        LoadZone::Bus => Vec::from_iter(0..nb), // each bus is its own zone
        LoadZone::All => vec![1; nb],           // make a single zone of all buses
        LoadZone::Area => bus.iter().map(|b| b.area).collect(),
    };
    let nz = *lz.iter().max().unwrap_or(&0); // number of load zones

    // fixed load at each bus, & initialize dispatchable
    let (p_df, q_df) = if want_fixed {
        let (sd_z, sd_i, sd_p) = make_sdzip(1.0, bus, mpopt);

        let vm = Arr::with_vec(bus.iter().map(|b| cmplx!(b.vm)).collect());

        let s_bus_d = &sd_p + &sd_i * &vm + &sd_z * &(&vm * &vm);
        let p_df = s_bus_d.iter().map(|s| s.re).collect(); // real power
        let q_df = if want_q {
            Some(s_bus_d.iter().map(|s| s.im).collect()) // reactive power
        } else {
            None
        };
        (p_df, q_df)
    } else {
        let p_df = vec![0.0; nb];
        let q_df = if want_q { Some(vec![0.0; nb]) } else { None };
        (p_df, q_df)
    };

    // dispatchable load at each bus
    let (p_dd, q_dd) = if want_disp {
        let mut p_dd = vec![0.0; nb];
        let mut q_dd = if want_q { Some(vec![0.0; nb]) } else { None };

        if let Some(gen) = gen {
            for g in gen {
                if g.is_load() && g.status {
                    if nominal {
                        p_dd[g.bus] += -g.pmin;
                        if want_q {
                            let q_dd = q_dd.as_mut().unwrap();
                            // TODO: (gen(ld, QMIN) == 0) .* gen(ld, QMAX) + (gen(ld, QMAX) == 0) .* gen(ld, QMIN)
                            if g.qmin == 0.0 {
                                q_dd[g.bus] += g.qmax;
                            } else if g.qmax == 0.0 {
                                q_dd[g.bus] += g.qmin;
                            }
                        }
                    } else {
                        p_dd[g.bus] += -g.pg;
                        if want_q {
                            let q_dd = q_dd.as_mut().unwrap();
                            q_dd[g.bus] += -g.qg;
                        }
                    }
                }
            }
        }

        (p_dd, q_dd)
    } else {
        let p_dd = vec![0.0; nb];
        let q_dd = if want_q { Some(vec![0.0; nb]) } else { None };
        (p_dd, q_dd)
    };

    // compute load sums
    if load_zone == LoadZone::Bus {
        // Pd = (Pdf + Pdd) .* (bus(:, BUS_TYPE) ~= NONE)
        let p_d = bus
            .iter()
            .enumerate()
            .map(|(i, b)| {
                if b.bus_type != BusType::NONE {
                    p_df[i] + p_dd[i]
                } else {
                    0.0
                }
            })
            .collect();

        let q_d = if want_q {
            Some(
                bus.iter()
                    .enumerate()
                    .map(|(i, b)| {
                        if b.bus_type != BusType::NONE {
                            let q_df = q_df.as_ref().unwrap();
                            let q_dd = q_dd.as_ref().unwrap();
                            q_df[i] + q_dd[i]
                        } else {
                            0.0
                        }
                    })
                    .collect(),
            )
        } else {
            None
        };

        (p_d, q_d)
    } else {
        let mut p_d = vec![0.0; nz];
        let mut q_d = if want_q { Some(vec![0.0; nz]) } else { None };

        for k in 0..nz {
            for (i, b) in bus.iter().enumerate() {
                if lz[i] == k && b.bus_type != BusType::NONE {
                    p_d[k] += p_df[i] + p_dd[i]; // TODO: sum(Pdf(idx)) + sum(Pdd(idx))

                    if want_q {
                        let qd = q_d.as_mut().unwrap();
                        let q_df = q_df.as_ref().unwrap();
                        let q_dd = q_dd.as_ref().unwrap();

                        qd[k] += q_df[i] + q_dd[i]; // TODO: sum(Qdf(idx)) + sum(Qdd(idx))
                    }
                }
            }
        }

        (p_d, q_d)
    }
}
