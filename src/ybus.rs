use std::f64::consts::PI;
use caseformat::{Branch, Bus};
use num_complex::Complex64;
use sparsetools::coo::Coo;

/// Builds the bus admittance matrix and branch admittance matrices.
pub fn make_ybus(
    base_mva: f64,
    bus: &[Bus],
    branch: &[Branch],
    yf_yt: bool,
) -> (
    Coo<usize, Complex64>,
    Option<Coo<usize, Complex64>>,
    Option<Coo<usize, Complex64>>,
) {
    let nb = bus.len();
    let nl = branch.len();

    // For each branch, compute the elements of the branch admittance matrix where:
    //
    //      | If |   | Yff  Yft |   | Vf |
    //      |    | = |          | * |    |
    //      | It |   | Ytf  Ytt |   | Vt |
    let mut y_bus = Coo::with_size(nb, nb);
    let mut y_f = if yf_yt {
        Some(Coo::with_size(nl, nb))
    } else {
        None
    };
    let mut y_t = if yf_yt {
        Some(Coo::with_size(nl, nb))
    } else {
        None
    };

    for (i, br) in branch.iter().enumerate() {
        let y_s = if br.is_on() {
            Complex64::new(1.0, 0.0) / Complex64::new(br.br_r, br.br_x)
        } else {
            Complex64::default()
        }; // series admittance
        let b_c = if br.is_on() { br.br_b } else { 0.0 }; // line charging susceptance
        let t = if br.tap == 0.0 { 1.0 } else { br.tap }; // default tap ratio = 1
        let tap = Complex64::from_polar(t, br.shift * PI / 180.0); // add phase shifters

        let y_tt = y_s + Complex64::new(0.0, b_c / 2.0);
        let y_ff = y_tt / (tap * tap.conj());
        let y_ft = -y_s / tap.conj();
        let y_tf = -y_s / tap;

        let (f, t) = (br.f_bus, br.t_bus);

        if yf_yt {
            let y_f = y_f.as_mut().unwrap();
            let y_t = y_t.as_mut().unwrap();

            y_f.push(i, f, y_ff);
            y_f.push(i, t, y_ft);

            y_t.push(i, f, y_tf);
            y_t.push(i, t, y_tt);
        }

        y_bus.push(f, f, y_ff);
        y_bus.push(f, t, y_ft);
        y_bus.push(t, f, y_tf);
        y_bus.push(t, t, y_tt);
    }

    let y_sh = bus
        .iter()
        .map(|b| Complex64::new(b.gs, b.bs) / Complex64::new(base_mva, 0.0))
        .collect::<Vec<Complex64>>();

    for (i, _) in bus.iter().enumerate() {
        y_bus.push(i, i, y_sh[i]);
    }
    (y_bus, y_f, y_t)

    /*
    // series admittance
    let y_s = branch.iter().map(|br| br.y_s()).collect::<Vec<Complex64>>();
    // line charging susceptance
    let b_c = branch
        .iter()
        .map(|br| if br.status { br.b } else { 0.0 })
        .collect::<Vec<f64>>();
    let tap = branch
        .iter()
        .map(|br| {
            let t = if br.tap == 0.0 { 1.0 } else { br.tap }; // default tap ratio = 1
            Complex64::from_polar(t, br.shift * PI / 180.0) // add phase shifters
        })
        .collect::<Vec<Complex64>>();

    // Ytt = Ys + 1j*Bc/2;
    let y_tt = y_s
        .iter()
        .zip(&b_c)
        .map(|(ys, bc)| ys + cmplx!(0.0, bc / 2.0))
        .collect::<Vec<Complex64>>();
    // Yff = Ytt ./ (tap .* conj(tap));
    let y_ff = y_tt
        .iter()
        .zip(&tap)
        .map(|(ytt, t)| ytt / (t * t.conj()))
        .collect::<Vec<Complex64>>();
    // Yft = - Ys ./ conj(tap);
    let y_ft = y_s
        .iter()
        .zip(&tap)
        .map(|(ys, t)| -ys / t.conj())
        .collect::<Vec<Complex64>>();
    // Ytf = - Ys ./ tap;
    let y_tf = y_s
        .iter()
        .zip(tap)
        .map(|(ys, t)| -ys / t)
        .collect::<Vec<Complex64>>();

    // Compute shunt admittance.
    //
    // If Psh is the real power consumed by the shunt at V = 1.0 p.u.
    // and Qsh is the reactive power injected by the shunt at V = 1.0 p.u.
    // then Psh - j Qsh = V * conj(Ysh * V) = conj(Ysh) = Gs - j Bs,
    // i.e. Ysh = Psh + j Qsh, so ...
    //Ysh = (bus(:, GS) + 1j * bus(:, BS)) / baseMVA; %% vector of shunt admittances
    let y_sh = bus
        .iter()
        .map(|b| b.y_sh(base_mva))
        .collect::<Vec<Complex64>>();

    let f = branch.iter().map(|br| br.from_bus).collect::<Vec<usize>>();
    let t = branch.iter().map(|br| br.to_bus).collect::<Vec<usize>>();

    {
        let lr = Vec::from_iter(0..nl);
        let l1 = vec![cmplx!(1.0); nl];

        // Build connection matrices.
        // let c_f = TriMat::from_triplets((nl, nb), lr.clone(), f, l1.clone()).to_csr();
        let c_f = sparse((nl, nb), &lr, &f, &l1);
        // let c_t = TriMat::from_triplets((nl, nb), lr.clone(), t, l1.clone()).to_csr();
        let c_t = sparse((nl, nb), &lr, &t, &l1);

        // Build Yf and Yt such that Yf * V is the vector of complex branch currents injected
        // at each branch's "from" bus, and Yt is the same for the "to" bus end
        // let y_f = TriMat::from_triplets((nl, nl), lr.clone(), lr.clone(), y_ff).to_csr() * &c_f
        //     + TriMat::from_triplets((nl, nl), lr.clone(), lr.clone(), y_ft).to_csr() * &c_t;
        let y_f = &(&spdiag(&y_ff) * &c_f) + &(&spdiag(&y_ft) * &c_t);
        // let y_t = TriMat::from_triplets((nl, nl), lr.clone(), lr.clone(), y_tf).to_csr() * c_f
        //     + TriMat::from_triplets((nl, nl), lr.clone(), lr.clone(), y_tt).to_csr() * c_t;
        let y_t = &(&spdiag(&y_tf) * &c_f) + &(&spdiag(&y_tt) * &c_t);

        // Branch admittances.
        let y_branch = &(&c_f.transpose_view() * &y_f) + &(&c_t.transpose_view() * &y_t);
        // Shunt admittance.
        let y_shunt = spdiag(&y_sh);

        let y_bus = &y_branch + &y_shunt;

        if !yf_yt {
            (y_bus, None, None)
        } else {
            (y_bus, Some(y_f), Some(y_t))
        }
    }
    */
    /*{
        let y_bus = TriMat::from_triplets(
            (nb, nb),
            [&f, &f, &t, &t].concat(),
            [&f, &t, &f, &t].concat(),
            [&y_ff, &y_ft, &y_tf, &y_tt].concat(),
        )
        .to_csr();

        if !yf_yt {
            (y_bus, None, None)
        } else {
            let i = [0..nl, 0..nl].concat();
            let y_f =
                TriMat::from_triplets((nl, nl), i, [&f, &t].concat(), [&y_ff, &y_ft].concat()).to_csr();
            let y_t =
                TriMat::from_triplets((nl, nl), i, [&f, &t].concat(), [&y_tf, &y_tt].concat()).to_csr();

            (y_bus, Some(y_f), Some(y_t))
        }

        (y_bus, y_f, y_t)
    }*/
}
