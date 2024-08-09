// Copyright (c) 2022-2024, Richard Lincoln. All rights reserved.

use num_complex::Complex64;

pub const J: Complex64 = Complex64 { re: 0.0, im: 1.0 };

#[macro_export]
macro_rules! cmplx {
    () => {
        num_complex::Complex64::new(0.0, 0.0)
    };
    ($arg1:expr) => {
        num_complex::Complex64::new($arg1, 0.0)
    };
    ($arg1:expr, $arg2:expr) => {
        num_complex::Complex64::new($arg1, $arg2)
    };
}

//
// /// Computes the infinity norm: `max(abs(a))`
// pub fn norm_inf(a: &[f64]) -> Option<f64> {
//     if a.is_empty() {
//         return None;
//     }
//     let mut max = f64::NEG_INFINITY;
//     a.iter().for_each(|v| {
//         let absvi = v.abs();
//         if absvi > max {
//             max = absvi
//         }
//     });
//     Some(max)
// }
//
// /// Returns the 2-norm (Euclidean) of `a`.
// pub fn norm(a: &[f64]) -> f64 {
//     let mut sqsum = 0.0;
//     for i in 0..a.len() {
//         sqsum += a[i] * a[i];
//     }
//     f64::sqrt(sqsum)
// }
