use caseformat::{Branch, Bus};
use sparsetools::coo::Coo;
use sparsetools::csr::CSR;
use std::f64::consts::PI;

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
pub fn make_b_dc(
    _base_mva: f64,
    bus: &[Bus],
    branch: &[Branch],
) -> (CSR<usize, f64>, CSR<usize, f64>, Vec<f64>, Vec<f64>) {
    let (rows, cols) = (branch.len(), bus.len());
    let nnz = 2 * branch.len();

    // Build Bf such that Bf * Va is the vector of real branch powers injected
    // at each branch's "from" bus.
    let mut b_f = Coo::with_capacity(rows, cols, nnz);
    // Build connection matrix Cft = Cf - Ct for line and from - to buses.
    let mut c_ft = Coo::with_capacity(rows, cols, nnz); // connection matrix

    fn br_b(br: &Branch) -> f64 {
        let b = if br.is_on() { 1.0 / br.br_x } else { 0.0 }; // series susceptance
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
    let c_ft = c_ft.to_csc();

    let mut b_bus = &c_ft.t() * &b_f;
    b_bus.sort_indexes();

    // Build phase shift injection vectors.
    let pfinj: Vec<f64> = branch
        .iter()
        .map(|br| br_b(br) * -br.shift * PI / 180.0)
        .collect(); // injected at the from bus ...

    let pbusinj = c_ft.transpose() * &pfinj; // ... and extracted at the to bus

    (b_bus, b_f, pbusinj, pfinj)
}
