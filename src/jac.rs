use anyhow::Result;
use num_complex::Complex64;
use sparsetools::coo::Coo;
use sparsetools::csc::CSC;
use sparsetools::csr::{CCSR, CSR};

use crate::debug::format_polar_vec;
use crate::{bus_types, d_sbus_d_v, MPC};

/// Forms the power flow Jacobian.
///
/// Input is a MATPOWER case struct and the bus admittance matrix.
/// Bus numbers must be consecutive beginning at 0 (i.e. internal ordering).
/// If the `full_jac` argument is true, it returns the full Jacobian
/// (sensitivities of all bus injections w.r.t all voltage angles/magnitudes)
/// as opposed to the reduced version used in the Newton power flow updates.
/// The units for all quantities are in per unit with radians for voltage
/// angles.
#[allow(non_snake_case)]
pub fn make_jac(
    mpc: &MPC,
    Ybus: &CSR<usize, Complex64>,
    full_jac: bool,
) -> Result<CSC<usize, f64>> {
    // extract voltage
    let mut V: Vec<Complex64> = mpc
        .bus
        .iter()
        .map(|bus| Complex64::from_polar(bus.vm, bus.va.to_radians()))
        .collect();

    // make sure we use generator setpoint voltage for PV and slack buses
    for g in mpc.gen.iter().filter(|g| g.is_on()) {
        let gbus = g.gen_bus;
        if mpc.bus[gbus].is_pv() || mpc.bus[gbus].is_ref() {
            V[gbus] = Complex64::new(g.vg / (V[gbus] * V[gbus]).norm(), 0.0);
        }
    }
    log::debug!("V0: {}", format_polar_vec(&V));

    // build Jacobian
    let (dSbus_dVa, dSbus_dVm) = d_sbus_d_v(Ybus, &V, false);

    let (j11, j12, j21, j22) = if full_jac {
        let j11 = dSbus_dVa.real();
        let j12 = dSbus_dVm.real();
        let j21 = dSbus_dVa.imag();
        let j22 = dSbus_dVm.imag();

        (j11, j12, j21, j22)
    } else {
        // get bus index lists of each type of bus
        let (_, pv, pq) = bus_types(&mpc.bus, &mpc.gen);
        let pv_pq = [pv, pq.clone()].concat();

        let j11 = dSbus_dVa.select(Some(&pv_pq), Some(&pv_pq))?.real();
        let j12 = dSbus_dVm.select(Some(&pv_pq), Some(&pq))?.real();
        let j21 = dSbus_dVa.select(Some(&pq), Some(&pv_pq))?.imag();
        let j22 = dSbus_dVm.select(Some(&pq), Some(&pq))?.imag();

        (j11, j12, j21, j22)
    };

    let J = Coo::compose([
        [&j11.to_coo(), &j12.to_coo()],
        [&j21.to_coo(), &j22.to_coo()],
    ])?
    .to_csc();
    log::trace!("J:\n{}", J.to_csr().to_table());

    Ok(J)
}
