use crate::mpc::MPC;
use crate::mpopt::Alg;
use crate::rpower::make_ybus;
use sprs::CsMat;

/// Builds the two matrices B prime and B double prime used in the fast
/// decoupled power flow.
fn make_b(mpc: &MPC, alg: Alg, double_prime: bool) -> (CsMat<f64>, Option<CsMat<f64>>) {
    // Form Bp (B prime).
    let mut bus = mpc.bus.clone(); // modify a copy of bus
    for b in bus.iter_mut() {
        b.bs = 0.0; // zero out shunts at buses
    }

    let b_p = {
        let mut branch = mpc.branch.clone(); // modify a copy of branch
        for br in branch.iter_mut() {
            br.b = 0.0; // zero out line charging shunts
            br.tap = 1.0; // cancel out taps
            if alg == Alg::FDXB {
                br.r = 0.0; // zero out line resistance
            }
        }
        let (y_p, _, _) = make_ybus(mpc.base_mva, &bus, &branch, false);
        y_p.map(|y| -y.im)
    };

    let b_pp = if double_prime {
        // Form Bpp (B double prime).
        let mut branch = mpc.branch.clone(); // modify a copy of branch
        for br in branch.iter_mut() {
            br.shift = 0.0; // zero out phase shifters
            if alg == Alg::FDBX {
                br.r = 0.0; // zero out line resistance
            }
        }
        let (y_pp, _, _) = make_ybus(mpc.base_mva, &bus, &branch, false);
        Some(y_pp.map(|y| -y.im))
    } else {
        None
    };

    (b_p, b_pp)
}
