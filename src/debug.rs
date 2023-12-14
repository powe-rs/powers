use num_complex::Complex64;
use pretty_dtoa::{dtoa, FmtFloatConfig};
use std::f64::consts::PI;

const FLOAT_CONFIG: FmtFloatConfig = FmtFloatConfig::default()
    .add_point_zero(false)
    .max_significant_digits(9);

pub fn format_f64_vec(v: &[f64]) -> String {
    let a: Vec<String> = v.iter().map(|f| dtoa(*f, FLOAT_CONFIG)).collect();
    format!("[{}]", a.join(", "))
}

fn format_complex(z: &Complex64) -> String {
    format!(
        "{}{}j{}",
        dtoa(z.re, FLOAT_CONFIG),
        if z.im.signum() < 0.0 { "-" } else { "+" },
        dtoa(z.im.abs(), FLOAT_CONFIG)
    )
    .to_string()
}

pub fn format_rect_vec(v: &[Complex64]) -> String {
    let a: Vec<String> = v.iter().map(|z| format_complex(z)).collect();
    format!("[{}]", a.join(", "))
}

fn format_polar(z: &Complex64) -> String {
    format!(
        "{}\u{2220}{}\u{00B0}",
        dtoa(z.norm(), FLOAT_CONFIG),
        dtoa(z.arg() * 180.0 / PI, FLOAT_CONFIG)
    )
    .to_string()
}

pub fn format_polar_vec(v: &[Complex64]) -> String {
    let a: Vec<String> = v.iter().map(|z| format_polar(z)).collect();
    format!("[{}]", a.join(", "))
}
