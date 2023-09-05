use crate::powers::norm_inf;
use crate::traits::LinearSolver;
use casecsv::{Branch, Bus};
use itertools::{izip, Itertools};
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
    p_bus: &[f64],
    va0: &[f64],
    ref_: &[usize],
    pv: &[usize],
    pq: &[usize],
    lin_solver: &dyn LinearSolver,
) -> Result<(Vec<f64>, bool), String> {
    let va_threshold = 1e5; // arbitrary threshold on |Va| for declaring failure

    // initialize result vector
    let mut va = va0.to_vec();
    let mut success = true; // successful by default

    // update angles for non-reference buses
    let pvpq = [pv, pq].concat();

    // Va([pv; pq]) = B([pv; pq], [pv; pq]) \ ...
    //                     (Pbus([pv; pq]) - B([pv; pq], ref) * Va0(ref));

    let b_pvpq = b_mat.select(Some(&pvpq), Some(&pvpq))?;
    let b_ref = b_mat.select(Some(&pvpq), Some(ref_))?;
    let p_bus_pvpq = pvpq.iter().map(|&i| p_bus[i]).collect_vec();
    let va_ref = ref_.iter().map(|&i| va0[i]).collect_vec();

    let rhs = izip!(p_bus_pvpq, b_ref * &va_ref)
        .map(|(p_bus, p_ref)| p_bus - p_ref)
        .collect_vec();

    let va_pvpq = lin_solver.solve(b_pvpq.to_csc(), &rhs)?;

    pvpq.iter()
        .enumerate()
        .for_each(|(i, &j)| va[j] = va_pvpq[i]);
    // set_slice(&mut va, &pvpq, &va_pvpq);

    // if va.abs().max() > va_threshold {
    if norm_inf(&va) > va_threshold {
        success = false;
    }

    Ok((va, success))
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
) -> (CSR<usize, f64>, CSR<usize, f64>, Vec<f64>, Vec<f64>) {
    let (rows, cols) = (branch.len(), bus.len());
    let nnz = 2 * branch.len();

    // Build Bf such that Bf * Va is the vector of real branch powers injected
    // at each branch's "from" bus.
    let mut b_f = Coo::empty(rows, cols, nnz);
    // Build connection matrix Cft = Cf - Ct for line and from - to buses.
    let mut c_ft = Coo::empty(rows, cols, nnz); // connection matrix

    fn br_b(br: &Branch) -> f64 {
        let b = if br.is_on() { 1.0 / br.x } else { 0.0 }; // series susceptance
        let tap = if br.tap == 0.0 { 1.0 } else { br.tap }; // default tap ratio = 1
        b / tap
    }
    for (i, br) in branch.iter().enumerate() {
        let b = br_b(br);

        let (f, t) = (br.f_bus, br.t_bus);

        b_f.push(i, f, b);
        b_f.push(i, t, -b);

        c_ft.push(i, f, 1.0);
        c_ft.push(i, t, -1.0);
    }
    let b_f = b_f.to_csr();
    let c_ft = c_ft.to_csr();

    let b_bus = &c_ft.t() * &b_f;

    // Build phase shift injection vectors.
    let pfinj = branch
        .iter()
        .map(|br| br_b(br) * -br.shift * PI / 180.0)
        .collect::<Vec<f64>>(); // injected at the from bus ...

    let pbusinj = &c_ft.t() * &pfinj; // ... and extracted at the to bus

    (b_bus, b_f, pbusinj, pfinj)
}
