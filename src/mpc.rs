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
    // HVDC interconnections.
    // dc_line: Vec<DCLine>,
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

    /// Loss zone (1-999).
    pub zone: usize,

    /// Maximum voltage magnitude (p.u.).
    pub vmax: f64,

    /// Minimum voltage magnitude (p.u.).
    pub vmin: f64,

    /// Lagrange multiplier on real power mismatch (u/MW).
    pub lam_p: f64,

    /// Lagrange multiplier on reactive power mismatch (u/MVAr).
    pub lam_q: f64,

    /// Kuhn-Tucker multiplier on upper voltage limit (u/p.u.).
    pub mu_vmax: f64,

    /// Kuhn-Tucker multiplier on lower voltage limit (u/p.u.).
    pub mu_vmin: f64,
}

impl Bus {
    pub(crate) fn y_sh(&self, base_mva: f64) -> Complex64 {
        Complex64::new(self.gs, self.bs) / Complex64::new(base_mva, 0.0)
    }
}

#[derive(Debug, PartialEq, Copy, Clone)]
pub enum GenStatus {
    InService = 0,
    OutOfService = 1,
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

    // pub status: GenStatus,
    pub status: bool,

    /// Maximum real power output (MW).
    pub pmax: f64,

    /// Minimum real power output (MW).
    pub pmin: f64,

    /// Lower real power output of PQ capability curve (MW).
    pub pc1: f64,

    /// Upper real power output of PQ capability curve (MW).
    pub pc2: f64,

    /// Minimum reactive power output at Pc1 (MVAr).
    pub qc1min: f64,

    /// Maximum reactive power output at Pc1 (MVAr).
    pub qc1max: f64,

    /// Minimum reactive power output at Pc2 (MVAr).
    pub qc2min: f64,

    /// Maximum reactive power output at Pc2 (MVAr).
    pub qc2max: f64,

    /// Ramp rate for load following/AGC (MW/min).
    pub ramp_agc: f64,

    /// Ramp rate for 10 minute reserves (MW).
    pub ramp_10: f64,

    /// Ramp rate for 30 minute reserves (MW).
    pub ramp_30: f64,

    /// Ramp rate for reactive power (2 sec timescale) (MVAr/min).
    pub ramp_q: f64,

    /// Area participation factor.
    pub apf: f64,

    /// Kuhn-Tucker multiplier on upper Pg limit (u/MW).
    pub mu_pmax: f64,

    /// Kuhn-Tucker multiplier on lower Pg limit (u/MW).
    pub mu_pmin: f64,

    /// Kuhn-Tucker multiplier on upper Qg limit (u/MVAr).
    pub mu_qmax: f64,

    /// Kuhn-Tucker multiplier on lower Qg limit (u/MVAr).
    pub mu_qmin: f64,

    /// Real power cost function.
    pub pcost: Option<GenCost>,

    /// Reactive power cost function;
    pub qcost: Option<GenCost>,
}

impl Gen {
    /// Checks for dispatchable loads.
    pub(crate) fn is_load(&self) -> bool {
        self.pmin < 0.0 && self.pmax == 0.0
    }
}

#[derive(Debug, PartialEq, Copy, Clone)]
pub enum CostModel {
    /// Piecewise linear cost model defined by breakpoints.
    PwLinear = 0,
    Polynomial = 1,
}

/// GenCost defines a generator cost function.
#[derive(Clone)]
pub struct GenCost {
    /// Cost model.
    pub model: CostModel,

    /// Startup cost in US dollars.
    pub startup: f64,

    /// Shutdown cost in US dollars.
    pub shutdown: f64,

    /// Number of cost coefficients to follow for polynomial
    /// cost function, or number of data points for piecewise linear.
    pub n_cost: usize,

    /// parameters defining total cost function f(p),
    /// units of f and p are $/hr and MW (or MVAr), respectively.
    ///
    /// PW_LINEAR: p0, f0, p1, f1, ..., pn, fn
    /// where p0 < p1 < ... < pn and the cost f(p) is defined by
    /// the coordinates (p0,f0), (p1,f1), ..., (pn,fn) of the
    /// end/break-points of the piecewise linear cost function
    ///
    /// POLYNOMIAL: cn, ..., c1, c0
    /// n+1 coefficients of an n-th order polynomial cost function,
    /// starting with highest order, where cost is
    /// f(p) = cn*p^n + ... + c1*p + c0
    pub cost: f64,
}

#[derive(Debug, PartialEq, Copy, Clone)]
pub enum BranchStatus {
    InService = 0,
    OutOfService = 1,
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

    /// MVA rating B (short term rating).
    pub rate_b: f64,

    /// MVA rating C (emergency rating).
    pub rate_c: f64,

    /// Transformer off nominal tap ratio.
    pub tap: f64,

    /// Transformer phase shift angle (degrees).
    pub shift: f64,

    /// Initial branch status.
    // pub status: BranchStatus,
    pub status: bool,

    /// Minimum angle difference, angle(Vf) - angle(Vt) (degrees).
    pub ang_min: f64,

    /// Maximum angle difference, angle(Vf) - angle(Vt) (degrees).
    pub ang_max: f64,

    /// Real power injected at "from" bus end (MW).
    pub pf: f64,

    /// Reactive power injected at "from" bus end (MVAr).
    pub qf: f64,

    /// Real power injected at "to" bus end (MW).
    pub pt: f64,

    /// Reactive power injected at "to" bus end (MVAr).
    pub qt: f64,

    /// Kuhn-Tucker multiplier on MVA limit at "from" bus (u/MVA).
    pub mu_sf: f64,

    /// Kuhn-Tucker multiplier on MVA limit at "to" bus (u/MVA).
    pub mu_st: f64,

    /// Kuhn-Tucker multiplier lower angle difference limit (u/degree).
    pub mu_ang_min: f64,

    /// Kuhn-Tucker multiplier upper angle difference limit (u/degree).
    pub mu_ang_max: f64,
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
