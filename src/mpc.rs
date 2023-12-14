use crate::order::Order;
use casecsv::{Branch, Bus, Gen};

#[derive(Clone, Default)]
pub struct MPC {
    pub name: String,

    /// System MVA base used for converting power into per-unit quantities.
    /// Default value is 100.
    pub base_mva: f64,

    /// Power system nodes, including static loads and shunts.
    pub bus: Vec<Bus>,

    /// Generators and dispatchable loads.
    pub gen: Vec<Gen>,

    /// Transmission lines/cables and transformers.
    pub branch: Vec<Branch>,

    // pub areas: Option<Vec<Area>>,
    pub(crate) order: Option<Order>,

    // pub(crate) a_mat: Option<CSR<usize, f64>>,
    // pub(crate) n_mat: Option<CSR<usize, f64>>,
    pub(crate) success: Option<bool>,
    pub(crate) iterations: Option<usize>,
}
