// Bus //

// define bus types
pub(crate) const PQ: usize = 1;
pub(crate) const PV: usize = 2;
pub(crate) const REF: usize = 3;
pub(crate) const NONE: usize = 4;

// define the indices
pub(crate) const BUS_I: usize = 0; // bus number (1 to 29997)
pub(crate) const BUS_TYPE: usize = 1; // bus type
pub(crate) const PD: usize = 2; // Pd, real power demand (MW)
pub(crate) const QD: usize = 3; // Qd, reactive power demand (MVAr)
pub(crate) const GS: usize = 4; // Gs, shunt conductance (MW at V = 1.0 p.u.)
pub(crate) const BS: usize = 5; // Bs, shunt susceptance (MVAr at V = 1.0 p.u.)
pub(crate) const BUS_AREA: usize = 6; // area number, 1-100
pub(crate) const VM: usize = 7; // Vm, voltage magnitude (p.u.)
pub(crate) const VA: usize = 8; // Va, voltage angle (degrees)
pub(crate) const BASE_KV: usize = 9; // baseKV, base voltage (kV)
pub(crate) const ZONE: usize = 10; // zone, loss zone (1-999)
pub(crate) const VMAX: usize = 11; // maxVm, maximum voltage magnitude (p.u.)
pub(crate) const VMIN: usize = 12; // minVm, minimum voltage magnitude (p.u.)

// included in opf solution, not necessarily in input
// assume objective function has units, u
pub(crate) const LAM_P: usize = 13; // Lagrange multiplier on real power mismatch (u/MW)
pub(crate) const LAM_Q: usize = 14; // Lagrange multiplier on reactive power mismatch (u/MVAr)
pub(crate) const MU_VMAX: usize = 15; // Kuhn-Tucker multiplier on upper voltage limit (u/p.u.)
pub(crate) const MU_VMIN: usize = 16; // Kuhn-Tucker multiplier on lower voltage limit (u/p.u.)

// Gen //

pub(crate) const GEN_BUS: usize = 0; // bus number
pub(crate) const PG: usize = 1; // Pg, real power output (MW)
pub(crate) const QG: usize = 2; // Qg, reactive power output (MVAr)
pub(crate) const QMAX: usize = 3; // Qmax, maximum reactive power output at Pmin (MVAr)
pub(crate) const QMIN: usize = 4; // Qmin, minimum reactive power output at Pmin (MVAr)
pub(crate) const VG: usize = 5; // Vg, voltage magnitude setpoint (p.u.)
pub(crate) const MBASE: usize = 6; // mBase, total MVA base of this machine, defaults to baseMVA
pub(crate) const GEN_STATUS: usize = 7; // status, 1 - machine in service, 0 - machine out of service
pub(crate) const PMAX: usize = 8; // Pmax, maximum real power output (MW)
pub(crate) const PMIN: usize = 9; // Pmin, minimum real power output (MW)
pub(crate) const PC1: usize = 10; // Pc1, lower real power output of PQ capability curve (MW)
pub(crate) const PC2: usize = 11; // Pc2, upper real power output of PQ capability curve (MW)
pub(crate) const QC1MIN: usize = 12; // Qc1min, minimum reactive power output at Pc1 (MVAr)
pub(crate) const QC1MAX: usize = 13; // Qc1max, maximum reactive power output at Pc1 (MVAr)
pub(crate) const QC2MIN: usize = 14; // Qc2min, minimum reactive power output at Pc2 (MVAr)
pub(crate) const QC2MAX: usize = 15; // Qc2max, maximum reactive power output at Pc2 (MVAr)
pub(crate) const RAMP_AGC: usize = 16; // ramp rate for load following/AGC (MW/min)
pub(crate) const RAMP_10: usize = 17; // ramp rate for 10 minute reserves (MW)
pub(crate) const RAMP_30: usize = 18; // ramp rate for 30 minute reserves (MW)
pub(crate) const RAMP_Q: usize = 19; // ramp rate for reactive power (2 sec timescale) (MVAr/min)
pub(crate) const APF: usize = 20; // area participation factor

// included in opf solution, not necessarily in input
// assume objective function has units, u
pub(crate) const MU_PMAX: usize = 21; // Kuhn-Tucker multiplier on upper Pg limit (u/MW)
pub(crate) const MU_PMIN: usize = 22; // Kuhn-Tucker multiplier on lower Pg limit (u/MW)
pub(crate) const MU_QMAX: usize = 23; // Kuhn-Tucker multiplier on upper Qg limit (u/MVAr)
pub(crate) const MU_QMIN: usize = 24; // Kuhn-Tucker multiplier on lower Qg limit (u/MVAr)

// Branch //

// define the indices
pub(crate) const F_BUS: usize = 0; // f, from bus number
pub(crate) const T_BUS: usize = 1; // t, to bus number
pub(crate) const BR_R: usize = 2; // r, resistance (p.u.)
pub(crate) const BR_X: usize = 3; // x, reactance (p.u.)
pub(crate) const BR_B: usize = 4; // b, total line charging susceptance (p.u.)
pub(crate) const RATE_A: usize = 5; // rateA, MVA rating A (long term rating)
pub(crate) const RATE_B: usize = 6; // rateB, MVA rating B (short term rating)
pub(crate) const RATE_C: usize = 7; // rateC, MVA rating C (emergency rating)
pub(crate) const TAP: usize = 8; // ratio, transformer off nominal turns ratio
pub(crate) const SHIFT: usize = 9; // angle, transformer phase shift angle (degrees)
pub(crate) const BR_STATUS: usize = 10; // initial branch status, 1 - in service, 0 - out of service
pub(crate) const ANGMIN: usize = 11; // minimum angle difference, angle(Vf) - angle(Vt) (degrees)
pub(crate) const ANGMAX: usize = 12; // maximum angle difference, angle(Vf) - angle(Vt) (degrees)

// included in power flow solution, not necessarily in input
pub(crate) const PF: usize = 13; // real power injected at "from" bus end (MW)
pub(crate) const QF: usize = 14; // reactive power injected at "from" bus end (MVAr)
pub(crate) const PT: usize = 15; // real power injected at "to" bus end (MW)
pub(crate) const QT: usize = 16; // reactive power injected at "to" bus end (MVAr)

// included in opf solution, not necessarily in input
// assume objective function has units, u
pub(crate) const MU_SF: usize = 17; // Kuhn-Tucker multiplier on MVA limit at "from" bus (u/MVA)
pub(crate) const MU_ST: usize = 18; // Kuhn-Tucker multiplier on MVA limit at "to" bus (u/MVA)
pub(crate) const MU_ANGMIN: usize = 19; // Kuhn-Tucker multiplier lower angle difference limit
pub(crate) const MU_ANGMAX: usize = 20; // Kuhn-Tucker multiplier upper angle difference limit
