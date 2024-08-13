use num_complex::Complex64;
use sparsetools::coo::Coo;
use sparsetools::csr::CSR;
use std::iter::zip;

/// Computes partial derivatives of current balance w.r.t. voltage.
///
/// The derivatives can be taken with respect to polar or cartesian coordinates
/// of voltage, depending on the `vcart` argument.
/// Returns two matrices containing partial derivatives of the complex bus
/// current balance w.r.t voltage magnitude and voltage angle, respectively
/// (for all buses).
pub fn d_imis_dv(
    s_bus: &Vec<Complex64>,
    y_bus: &CSR<usize, Complex64>,
    v: &[Complex64],
    vcart: bool,
) -> anyhow::Result<(CSR<usize, Complex64>, CSR<usize, Complex64>)> {
    let n = v.len();

    let (d_imis_dv1, d_imis_dv2) = if vcart {
        /*
        diagSV2c = sparse(1:n, 1:n, conj(Sbus./(V.^2)), n, n);
        */
        let diag_sv2c = Coo::new(
            n,
            n,
            (0..n).collect(),
            (0..n).collect(),
            // (0..n).map(|i| (s_bus[i] / (v[i] * v[i]).conj())).collect(),
            zip(s_bus, v)
                .map(|(s_bus, v)| s_bus / (v * v).conj())
                .collect(),
        )?
        .to_csr();
        let d_imis_dvr = y_bus + &diag_sv2c;
        let d_imis_dvi = Complex64::i() * (y_bus - diag_sv2c);

        (d_imis_dvr, d_imis_dvi)
    } else {
        // let v_m = Vec::from_real(&v.norm());
        // let v_norm: Vec<Complex64> = (0..n).map(|i| v[i] / cmplx!(v[i].norm())).collect();
        let v_norm = v
            .iter()
            .map(|v| v / Complex64::new(v.norm(), 0.0))
            .collect();
        // let i_bus: Vec<Complex64> = (0..n).map(|i| s_bus[i] / v[i]).collect();
        let i_bus = zip(s_bus, v).map(|(s_bus, v)| s_bus / v).collect();
        let i_bus_vm: Vec<Complex64> = zip(&i_bus, v)
            .map(|(i_bus, v)| i_bus / Complex64::new(v.norm(), 0.0))
            .collect();
        /*
        diagV       = sparse(1:n, 1:n, V, n, n);
        diagIbus    = sparse(1:n, 1:n, Ibus, n, n);
        diagIbusVm  = sparse(1:n, 1:n, Ibus./Vm, n, n);
        diagVnorm   = sparse(1:n, 1:n, V./abs(V), n, n);
        */
        let diag_v = CSR::with_diagonal(v.to_vec());
        let diag_ibus = CSR::with_diagonal(i_bus);
        let diag_ibus_vm = CSR::with_diagonal(i_bus_vm);
        // let diag_v_norm = CSR::with_diag((v / v.norm()).to_vec());
        let diag_v_norm = CSR::with_diagonal(v_norm);

        let d_imis_dva = Complex64::i() * (y_bus * diag_v - diag_ibus);
        let d_imis_dvm = y_bus * diag_v_norm + diag_ibus_vm; // dImis/dVm

        (d_imis_dva, d_imis_dvm)
    };

    Ok((d_imis_dv1, d_imis_dv2))
}
