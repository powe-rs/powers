use num_complex::Complex64;

/// MPC is a MATPOWER case that models a power system as a directed graph
/// structure.
pub struct MPC {
    /// System MVA base used for converting power into per-unit quantities.
    /// Default value is 100.
    pub base_mva: f64,

    /// Power system nodes, including static loads and shunts.
    pub bus: Vec<Bus>,

    /// Generators and dispatchable loads.
    pub gen: Vec<Gen>,

    /// Transmission lines/cables and transformers.
    pub branch: Vec<Branch>,
}

#[derive(Debug, PartialEq, Copy, Clone)]
pub enum BusType {
    /// Fixed active and reactive power.
    PQ = 0,
    /// Fixed voltage magnitude and active power.
    PV = 1,
    /// Reference voltage angle. Slack active and reactive power.
    REF = 2,
    /// Isolated bus.
    NONE = 3,
}

/// Bus is a node in the power system graph structure.
/// Static loads and shunts are included in the Bus definition.
#[derive(Clone)]
pub struct Bus {
    /// Bus number.
    pub i: usize,

    pub bus_type: BusType,

    /// Real power demand (MW).
    pub pd: f64,

    /// Reactive power demand (MVAr).
    pub qd: f64,

    /// Shunt conductance (MW at V = 1.0 p.u.).
    pub gs: f64,

    /// Shunt susceptance (MVAr at V = 1.0 p.u.).
    pub bs: f64,

    /// Area number, 1-100.
    pub area: usize,

    /// Voltage magnitude (p.u.).
    pub vm: f64,

    /// Voltage angle (degrees).
    pub va: f64,

    /// Base voltage (kV).
    pub base_kv: f64,

    /// Maximum voltage magnitude (p.u.).
    pub vmax: f64,

    /// Minimum voltage magnitude (p.u.).
    pub vmin: f64,
}

impl Bus {
    pub(crate) fn y_sh(&self, base_mva: f64) -> Complex64 {
        Complex64::new(self.gs, self.bs) / Complex64::new(base_mva, 0.0)
    }
}

/// Gen is a generator or dispatchable load.
#[derive(Clone)]
pub struct Gen {
    /// Bus number.
    pub bus: usize,

    /// Real power output (MW).
    pub pg: f64,

    /// Reactive power output (MVAr).
    pub qg: f64,

    /// Maximum reactive power output (MVAr).
    pub qmax: f64,

    /// Minimum reactive power output (MVAr).
    pub qmin: f64,

    /// Voltage magnitude setpoint (p.u.).
    pub vg: f64,

    /// Total MVA base of this machine, defaults to base_mva.
    pub mbase: f64,

    pub status: bool,

    /// Maximum real power output (MW).
    pub pmax: f64,

    /// Minimum real power output (MW).
    pub pmin: f64,
}

impl Gen {
    /// Checks for dispatchable loads.
    pub(crate) fn is_load(&self) -> bool {
        self.pmin < 0.0 && self.pmax == 0.0
    }
}

/// Branch represents either a transmission line/cable or a two winding
/// transformer.
#[derive(Clone)]
pub struct Branch {
    /// From bus number.
    pub from_bus: usize,

    /// To bus number.
    pub to_bus: usize,

    /// Resistance (p.u.).
    pub r: f64,

    /// Reactance (p.u.).
    pub x: f64,

    /// Total line charging susceptance (p.u.).
    pub b: f64,

    /// MVA rating A (long term rating).
    pub rate_a: f64,

    /// Transformer off nominal tap ratio.
    pub tap: f64,

    /// Transformer phase shift angle (degrees).
    pub shift: f64,

    /// Initial branch status.
    pub status: bool,

    /// Real power injected at "from" bus end (MW).
    pub pf: f64,

    /// Reactive power injected at "from" bus end (MVAr).
    pub qf: f64,

    /// Real power injected at "to" bus end (MW).
    pub pt: f64,

    /// Reactive power injected at "to" bus end (MVAr).
    pub qt: f64,
}

impl Branch {
    pub(crate) fn y_s(&self) -> Complex64 {
        if !self.status {
            Complex64::new(0.0, 0.0)
        } else {
            Complex64::new(1.0, 0.0) / Complex64::new(self.r, self.x)
        }
    }
}
