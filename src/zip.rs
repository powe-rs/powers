use caseformat::Bus;
use num_complex::Complex64;

/// Builds vectors of nominal complex bus power demands for ZIP loads.
///
/// Returns a struct with three fields, each an nb x 1 vectors. The fields
/// 'z', 'i' and 'p' correspond to the nominal p.u. complex power
/// (at 1 p.u. voltage magnitude) of the constant impedance, constant current,
/// and constant power portions, respectively of the ZIP load model.
pub fn make_sdzip(
    base_mva: f64,
    bus: &[Bus],
    pw: Option<[f64; 3]>,
    qw: Option<[f64; 3]>,
) -> (Vec<Complex64>, Vec<Complex64>, Vec<Complex64>) {
    let pw = pw.unwrap_or([1.0, 0.0, 0.0]);
    let qw = qw.unwrap_or(pw);

    let base_mva = Complex64::new(base_mva, 0.0);

    let mut sd_z = Vec::with_capacity(bus.len());
    let mut sd_i = Vec::with_capacity(bus.len());
    let mut sd_p = Vec::with_capacity(bus.len());
    for b in bus {
        sd_z.push(Complex64::new(b.pd * pw[2], b.qd * qw[2]) / base_mva);
        sd_i.push(Complex64::new(b.pd * pw[1], b.qd * qw[1]) / base_mva);
        sd_p.push(Complex64::new(b.pd * pw[0], b.qd * qw[0]) / base_mva);
    }
    (sd_z, sd_i, sd_p)
}
