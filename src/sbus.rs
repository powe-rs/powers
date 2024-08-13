use caseformat::{Bus, Gen};
use num_complex::Complex64;
use sparsetools::csr::{CCSR, CSR};

// pub type SBus = fn(v_m: &[f64]) -> (Arr<Complex64>, Option<CSR<usize, Complex64>>);

pub trait SBus {
    fn s_bus(&self, v_m: &[f64]) -> Vec<Complex64>;
    fn d_sbus_d_vm(&self, v_m: &[f64]) -> CSR<usize, Complex64>;
}

pub struct MakeSBus<'a> {
    pub base_mva: f64,
    pub bus: &'a [Bus],
    pub gen: &'a [Gen],
    pub pw: Option<[f64; 3]>,
    pub qw: Option<[f64; 3]>,
}

impl<'a> SBus for MakeSBus<'a> {
    fn s_bus(&self, v_m: &[f64]) -> Vec<Complex64> {
        make_sbus(
            self.base_mva,
            self.bus,
            self.gen,
            self.pw,
            self.qw,
            Some(v_m),
            None,
        )
    }

    fn d_sbus_d_vm(&self, v_m: &[f64]) -> CSR<usize, Complex64> {
        make_d_sbus_d_vm(
            self.base_mva,
            self.bus,
            self.gen,
            self.pw,
            self.qw,
            Some(v_m),
            None,
        )
    }
}

/// Builds the vector of complex bus power injections.
///
/// Returns the vector of complex bus power injections, that is, generation
/// minus load. Power is expressed in per unit. If the MPOPT and VM arguments
/// are present it evaluates any ZIP loads based on the provided voltage
/// magnitude vector. If VM is empty, it assumes nominal voltage. If SG is
/// provided, it is a complex ng x 1 vector of generator power injections in
/// p.u., and overrides the PG and QG columns in GEN, using GEN only for
/// connectivity information.
///
///   [SBUS, DSBUS_DVM] = MAKESBUS(BASEMVA, BUS, GEN, MPOPT, VM)
///
/// See also MAKEYBUS.
pub fn make_sbus(
    base_mva: f64,
    bus: &[Bus],
    gen: &[Gen],
    pw: Option<[f64; 3]>,
    qw: Option<[f64; 3]>,
    vm: Option<&[f64]>,
    sg: Option<&[Complex64]>,
) -> Vec<Complex64> {
    let nb = bus.len();
    let base_mva = Complex64::new(base_mva, 0.0);

    // Form net complex bus power injection vector
    // (power injected by generators + power injected by loads).
    let mut s_bus = vec![Complex64::default(); nb];

    if let Some(sg) = sg {
        gen.iter()
            .zip(sg)
            .filter(|(g, _)| g.is_on())
            .for_each(|(g, s_pu)| {
                s_bus[g.gen_bus] += s_pu;
            });
    } else {
        gen.iter().filter(|g| g.is_on()).for_each(|g| {
            s_bus[g.gen_bus] += Complex64::new(g.pg, g.qg) / base_mva;
        });
    }

    // get load parameters
    let pw = pw.unwrap_or([1.0, 0.0, 0.0]);
    let qw = qw.unwrap_or(pw);

    bus.iter()
        .filter(|b| b.pd != 0.0 || b.qd != 0.0)
        .for_each(|b| {
            // Compute per-bus loads in p.u.
            let sd_z = Complex64::new(b.pd * pw[2], b.qd * qw[2]) / base_mva;
            let sd_i = Complex64::new(b.pd * pw[1], b.qd * qw[1]) / base_mva;
            let sd_p = Complex64::new(b.pd * pw[0], b.qd * qw[0]) / base_mva;

            let vm_i = match vm {
                Some(vm) => Complex64::new(vm[b.bus_i], 0.0),
                None => Complex64::new(1.0, 0.0),
            };

            let sd = sd_p + (sd_i * vm_i) + (sd_z * (vm_i * vm_i));

            s_bus[b.bus_i] -= sd;
        });

    s_bus
}

/// Computes partial derivatives of power injection w.r.t. voltage.
///
/// The derivatives can be take with respect to polar or cartesian coordinates
/// of voltage, depending on the 3rd argument.
pub fn d_sbus_d_v(
    y_bus: &CSR<usize, Complex64>,
    v: &[Complex64],
    cartesian: bool,
) -> (CSR<usize, Complex64>, CSR<usize, Complex64>) {
    let i_bus = y_bus * v;

    let diag_v = CSR::<usize, Complex64>::with_diagonal(v.to_vec());
    let diag_i_bus = CSR::<usize, Complex64>::with_diagonal(i_bus);

    if cartesian {
        // dSbus/dVr = conj(diagIbus) + diagV * conj(Ybus)
        // dSbus/dVi = 1j * (conj(diagIbus) - diagV * conj(Ybus))

        let d_sbus_d_vr = diag_i_bus.conj() + &diag_v * y_bus.conj();
        let d_sbus_d_vi = (diag_i_bus.conj() - &diag_v * y_bus.conj()) * Complex64::i();

        (d_sbus_d_vr, d_sbus_d_vi)
    } else {
        let v_norm = v
            .iter()
            .map(|v| v / Complex64::new(v.norm(), 0.0))
            .collect();
        let diag_v_norm = CSR::<usize, Complex64>::with_diagonal(v_norm);

        // dSbus/dVa = 1j * diagV * conj(diagIbus - Ybus * diagV)
        // dSbus/dVm = diagV * conj(Ybus * diagVnorm) + conj(diagIbus) * diagVnorm

        let mut d_sbus_d_va = &diag_v * (&diag_i_bus - y_bus * &diag_v).conj() * Complex64::i();
        let d_sbus_d_vm =
            &diag_v * (y_bus * &diag_v_norm).conj() + diag_i_bus.conj() * &diag_v_norm;

        d_sbus_d_va.sort_indexes();

        (d_sbus_d_va, d_sbus_d_vm)
    }
}

/// Computes the partial derivative of the bus injections with respect to
/// voltage magnitude. If `vm` is empty, it assumes no voltage dependence
/// and returns a sparse zero matrix.
pub fn make_d_sbus_d_vm(
    base_mva: f64,
    bus: &[Bus],
    _gen: &[Gen],
    pw: Option<[f64; 3]>,
    qw: Option<[f64; 3]>,
    vm: Option<&[f64]>,
    _sg: Option<&[Complex64]>,
) -> CSR<usize, Complex64> {
    let nb = bus.len();
    let base_mva = Complex64::new(base_mva, 0.0);

    // get load parameters
    let pw = pw.unwrap_or([1.0, 0.0, 0.0]);
    let qw = qw.unwrap_or(pw);

    if let Some(vm) = vm {
        const TWO: Complex64 = Complex64 { re: 2.0, im: 0.0 };
        let diag: Vec<Complex64> = bus
            .iter()
            .map(|b| {
                if b.pd != 0.0 || b.qd != 0.0 {
                    let sd_z = Complex64::new(b.pd * pw[2], b.qd * qw[2]) / base_mva;
                    let sd_i = Complex64::new(b.pd * pw[1], b.qd * qw[1]) / base_mva;

                    let vm_i = Complex64::new(vm[b.bus_i], 0.0);

                    -(sd_i + TWO * vm_i * sd_z)
                } else {
                    Complex64::default()
                }
            })
            .collect();
        CSR::with_diagonal(diag)
    } else {
        CSR::with_size(nb, nb)
    }
}
