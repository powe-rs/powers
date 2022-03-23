use ndarray::ArrayView1;
use num_complex::Complex64;
use sprs::{CsMat, CsMatBase, CsMatView};
use std::collections::HashMap;

/// Select specific rows and columns from A.
///
/// Input: A is assumed to have canonical format.
///
/// A translation of `get_csr_submatrix` from the sparsetools C++ library in SciPy.
///
/// https://github.com/scipy/scipy/tree/main/scipy/sparse/sparsetools
pub(crate) fn csr_select(
    // a_p: &[usize],
    // a_j: &[usize],
    // a_x: &[Complex64],
    a_csr: CsMatView<Complex64>,
    row_idx: &[usize],
    col_idx: &[usize],
    real: bool,
) -> CsMat<f64> {
    assert!(a_csr.is_csr());

    // let a_p = a_csr.indptr();
    // let a_j = a_csr.indices();
    // let a_x = a_csr.data();

    let (a_p, a_j, a_x) = a_csr.into_raw_storage();

    let new_n_row = row_idx.len();
    let mut new_nnz = 0;

    for &row in row_idx {
        let start = a_p[row];
        let end = a_p[row + 1];

        let mut col_ptrs = HashMap::new();
        for jj in start..end {
            col_ptrs.insert(a_j[jj], jj);
        }

        for col in col_idx {
            if col_ptrs.contains_key(col) {
                new_nnz += 1;
            }
        }
    }

    let mut b_p = vec![0; new_n_row + 1];
    let mut b_j = vec![0; new_nnz];
    // let mut b_x = vec![Complex64::new(0.0, 0.0); new_nnz];
    let mut b_x = vec![0.0; new_nnz];

    b_p[0] = 0;
    let mut kk = 0;
    for (r, &row) in row_idx.iter().enumerate() {
        let start = a_p[row];
        let end = a_p[row + 1];

        let mut col_ptrs = HashMap::new();
        for jj in start..end {
            col_ptrs.insert(a_j[jj], jj);
        }

        for (c, col) in col_idx.iter().enumerate() {
            if let Some(&jj) = col_ptrs.get(col) {
                b_j[kk] = c;
                b_x[kk] = if real { a_x[jj].re } else { a_x[jj].im };
                kk += 1;
            }
        }
        b_p[r + 1] = kk;
    }

    CsMatBase::new(a_csr.shape(), b_p, b_j, b_x)
}

pub(crate) fn norm_inf(a: ArrayView1<f64>) -> f64 {
    a.iter()
        .max_by(|&a, &b| a.partial_cmp(b).unwrap())
        .unwrap()
        .abs()
}
