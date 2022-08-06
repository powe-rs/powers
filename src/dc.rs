use crate::mpc::{Branch, Bus};
use crate::traits::LinearSolver;
use densetools::arr::Arr;
use densetools::slice::set_slice;
use sparsetools::coo::Coo;
use sparsetools::csr::CSR;
use std::f64::consts::PI;

/// Solves a DC power flow.
///
/// Solves for the bus voltage angles at all but the reference bus,
/// given the full system B matrix and the vector of bus real power injections,
/// the initial vector of bus voltage angles (in radians), and column vectors
/// with the lists of bus indices for the swing bus, PV buses, and PQ buses,
/// respectively. Returns a vector of bus voltage angles in radians.
pub(crate) fn dc_pf(
    b_mat: &CSR<usize, f64>,
    p_bus: &Arr<f64>,
    v_a0: &Arr<f64>,
    ref_: &[usize],
    pv: &[usize],
    pq: &[usize],
    lin_solver: &dyn LinearSolver,
) -> Result<(Arr<f64>, bool), String> {
    let v_a_threshold = 1e5; // arbitrary threshold on |Va| for declaring failure

    // initialize result vector
    let mut v_a = v_a0.clone();
    let mut success = true; // successful by default

    // update angles for non-reference buses
    let pvpq = [pv, pq].concat();

    // Va([pv; pq]) = B([pv; pq], [pv; pq]) \ ...
    //                     (Pbus([pv; pq]) - B([pv; pq], ref) * Va0(ref));

    let a_mat = b_mat.select(Some(&pvpq), Some(&pvpq))?;
    let b_ref = b_mat.select(Some(&pvpq), Some(ref_))?;
    let p_bus_pvpq = p_bus.select(&pvpq);
    let v_a_ref = v_a0.select(&ref_);

    let rhs = p_bus_pvpq - (b_ref * v_a_ref);

    let v_a_pvpq = lin_solver.solve(a_mat.to_csc(), &rhs)?;

    set_slice(&mut v_a, &pvpq, &v_a_pvpq);

    if v_a.abs().max() > v_a_threshold {
        success = false;
    }

    Ok((v_a, success))
}

/// Builds the B matrices and phase shift injections for DC power flow.
///
/// Returns the B matrices and phase shift injection vectors needed for
/// a DC power flow. The bus real power injections are related to bus
/// voltage angles by
///     P = BBUS * Va + PBUSINJ
/// The real power flows at the from end the lines are related to the bus
/// voltage angles by
///     Pf = BF * Va + PFINJ
/// Does appropriate conversions to p.u.
/// Bus numbers must be consecutive beginning at 1 (i.e. internal ordering).
pub(crate) fn make_b_dc(
    _base_mva: f64,
    bus: &[Bus],
    branch: &[Branch],
) -> (CSR<usize, f64>, CSR<usize, f64>, Arr<f64>, Arr<f64>) {
    let (rows, cols) = (branch.len(), bus.len());
    let nnz = 2 * branch.len();

    // Build Bf such that Bf * Va is the vector of real branch powers injected
    // at each branch's "from" bus.
    let mut b_f = Coo::empty(rows, cols, nnz);
    // Build connection matrix Cft = Cf - Ct for line and from - to buses.
    let mut c_ft = Coo::empty(rows, cols, nnz); // connection matrix

    fn br_b(br: &Branch) -> f64 {
        let b = if br.status { 1.0 / br.x } else { 0.0 }; // series susceptance
        let tap = if br.tap == 0.0 { 1.0 } else { br.tap }; // default tap ratio = 1
        b / tap
    }
    for (i, br) in branch.iter().enumerate() {
        let b = br_b(br);

        let (f, t) = (br.from_bus, br.to_bus);

        b_f.push(i, f, b);
        b_f.push(i, t, -b);

        c_ft.push(i, f, 1.0);
        c_ft.push(i, t, -1.0);
    }
    let b_f = b_f.to_csr();
    let c_ft = c_ft.to_csr();

    let b_bus = &c_ft.t() * &b_f;

    // Build phase shift injection vectors.
    let pfinj = Arr::with_vec(
        branch
            .iter()
            .map(|br| br_b(br) * -br.shift * PI / 180.0)
            .collect(),
    ); // injected at the from bus ...

    let pbusinj = &c_ft.t() * &pfinj; // ... and extracted at the to bus

    (b_bus, b_f, pbusinj, pfinj)
}
