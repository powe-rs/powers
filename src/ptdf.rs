use crate::dc::make_b_dc;
use crate::dense::{find, Mat};
use crate::mpc::MPC;
use anyhow::{format_err, Result};
use casecsv::Branch;
use sparsetools::coo::Coo;
use spsolve::Solver;

// #[derive(Clone)]
enum PTDFSlack {
    /// Single slack bus. The reference bus is used by default.
    Single(usize),

    /// Weights specifying the proportion of the slack taken up at each bus.
    Weights(Vec<f64>),

    /// An `nb x nb` matrix where each column specifies how the slack should
    /// be handled for injections at that bus.
    ///
    /// This option only applies when computing the full PTDF matrix.
    Matrix(Mat<f64>),
}

// #[derive(Clone)]
enum PTDFBus {
    /// A matrix with `nb` rows whose columns each sum to zero, where each column
    /// defines a specific (slack independent) transfer.
    ///
    /// The returned PTDF matrix will have the same number of columns.
    Txfr(Mat<f64>),

    /// A column vector of bus indices. The columns of the returned PTDF matrix
    /// correspond to these indices, but they are computed individually rather
    /// than computing the full PTDF matrix and selecting the desired columns.
    BusIdx(Vec<usize>),
}

impl PTDFBus {
    fn bus_idx(self) -> Vec<usize> {
        match self {
            PTDFBus::BusIdx(v) => v,
            _ => {
                panic!("called `PTDFBus::bus_idx()` on a `Txfr` value")
            }
        }
    }
}

/// Builds the DC PTDF matrix for a given choice of slack.
pub fn make_ptdf(
    mpc: &MPC,
    mut slack: Option<PTDFSlack>,
    bus_opt: Option<PTDFBus>,
    solver: &dyn Solver<usize, f64>,
) -> Result<Mat<f64>> {
    if slack.is_none() {
        for (i, bus) in mpc.bus.iter().enumerate() {
            if bus.is_ref() {
                slack = Some(PTDFSlack::Single(i));
                break;
            }
        }
        if slack.is_none() {
            return Err(format_err!("slack must not be none"));
        }
    }

    // compute full PTDF?
    let mut compute_full_h = true;
    let nb = mpc.bus.len();
    // let nbr = mpc.branch.len();
    let mut txfr = false; // default assumes not just for specific transfers
    let mut d_p: Option<Mat<f64>> = None;
    if let Some(ptdf_bus) = bus_opt {
        compute_full_h = false;
        if let PTDFBus::Txfr(txfr_mat) = ptdf_bus {
            if txfr_mat.rows() != nb {
                return Err(format_err!(
                    "txfr: matrix must have nb rows; nb = {} rows = {}",
                    txfr_mat.rows(),
                    nb
                ));
            }
            for r in 0..txfr_mat.rows() {
                let r_sum: f64 = txfr_mat.row(r).sum();
                if r_sum != 0.0 {
                    return Err(format_err!(
                        "txfr: columns of each row must sum to zero; row = {} sum = {}",
                        r,
                        r_sum
                    ));
                }
            }
            txfr = true; // cols of H for specific (slack independent) transfers
            d_p = Some(txfr_mat);
        }
    }

    // set the slack bus to be used to compute initial PTDF
    let mut slack_bus = 0; // use bus 1 for temp slack bus
    if let Some(ptdf_slack) = slack {
        if let PTDFSlack::Single(slack_bus_index) = ptdf_slack {
            slack_bus = slack_bus_index
        }
    }

    let noref: Vec<usize> = (1..nb).collect(); // use bus 1 for voltage angle reference
    let noslack: Vec<usize> = (0..nb).filter(|&i| i != slack_bus).collect();
    // let noslack = find(
    //     (0..nb)
    //         .map(|i| if i != slack_bus { 1 } else { 0 })
    //         .collect(),
    // );
    // let noslack = find(&ne(&Arr::range(nb), slack_bus));

    // check that bus numbers are equal to indices to bus (one set of bus numbers)
    for (i, bus) in mpc.bus.iter().enumerate() {
        if bus.bus_i != i {
            return Err(format_err!(
                "make_ptdf: buses must be numbered consecutively in bus matrix"
            ));
        }
    }

    //  compute PTDF for single slack_bus  //
    let (b_bus, b_f, _, _) = make_b_dc(mpc.base_mva, &mpc.bus, &mpc.branch);

    // set up injections/transfers
    let (nbi0, nbi, bidx) = if compute_full_h {
        // full H for all columns
        let nbi = nb; // number of buses of interest (i.e. all of them)
        d_p = Some(Mat::identity(nb).col_major().build()?); // TODO: sparse RHS
        (0, nbi, None)
    } else {
        if txfr {
            // for specific transfers
            let nbi = d_p.unwrap().cols(); // number of transfers
            (0, nbi, None)
        } else {
            let bus_idx = bus_opt.unwrap().bus_idx();

            // for subset of columns of H for specific buses
            // if necessary add missing slacks to bus index list
            let nbi0 = bus_idx.len(); // number of original buses of interest
            let mut bidx: Vec<usize> = bus_idx.clone();
            if let Some(ptdf_slack) = &slack {
                if let PTDFSlack::Weights(weights) = ptdf_slack {
                    let slacks = find(&weights); // get all slack buses
                    for sw in slacks {
                        if !bus_idx.contains(&sw) {
                            bidx.push(sw);
                        }
                    }
                }
            }
            let nbi = bidx.len(); // number of buses of interest

            // define the equivalent transfer, each column is a single transfer
            //d_p = accumarray([bidx (1:nbi)'], 1, [nb, nbi]);
            // let mut d_p = Coo::empty(nb, nbi, 2 * nbi);
            let mut dp = Mat::new(nb, nbi).col_major().build()?;
            for (_i, _j) in bidx.iter().enumerate() {
                // TODO
                // dp.push(i, j, 1.0);
            }
            d_p = Some(dp);
            (nbi0, nbi, Some(bidx))
        }
    };

    // solve for change in voltage angles
    let d_theta = {
        let bbus_noslack = b_bus.select(Some(&noslack), Some(&noref))?;
        let mut dp_noslack = d_p.unwrap().select_rows(&noslack);

        solver.solve(
            bbus_noslack.cols(),
            bbus_noslack.colidx(),
            bbus_noslack.rowptr(),
            bbus_noslack.values(),
            dp_noslack.values_mut(),
            true,
        )?;
        let d_theta_noref = dp_noslack;

        // let d_theta = Mat::zeros(nb, nbi);
        // let d_theta = d_theta_noref;
        let mut d_theta: Mat<f64> = Mat::new(nb, nbi).col_major().build()?;
        d_theta.values_mut().extend(d_theta_noref.values());
        d_theta
    };

    // compute corresponding change in branch flows
    let mut h_values: Vec<f64> = b_f * d_theta.values();
    let mut h_mat = Mat::new(nb, nbi).values(h_values).col_major().build()?;

    // distribute slack, if requested //
    if !txfr {
        if let Some(ptdf_slack) = &slack {
            match ptdf_slack {
                PTDFSlack::Weights(weights) => {
                    // slack is a vector of weights

                    // normalize weights
                    let slack_sum: f64 = weights.iter().sum();
                    let slack: Vec<f64> = weights.iter().map(|w| w / slack_sum).collect();

                    // conceptually, we want to do ...
                    //    H = H * (eye(nb,nb) - slack * ones(1, nb));
                    // ... we just do it more efficiently

                    if compute_full_h {
                        let v = h_mat.mat_vec(&slack);
                        for k in 0..nb {
                            h_mat.col_mut(k).enumerate().for_each(|(i, hi)| *hi -= v[i]);
                        }
                    } else {
                        // let v = H * select(&slack, bidx.unwrap());
                        let slack_bidx: Vec<f64> =
                            bidx.unwrap().iter().map(|&i| slack[i]).collect();
                        let v = h_mat.mat_vec(&slack_bidx);
                        for k in 0..nbi {
                            h_mat.col_mut(k).enumerate().for_each(|(i, hi)| *hi -= v[i]);
                        }
                        h_mat = h_mat.select_cols(&(0..nbi0).collect::<Vec<usize>>());
                        // remove temp cols added for missing slacks
                    }
                }
                PTDFSlack::Matrix(slack) => {
                    // normalize weights
                    // let slack_sum: f64 = slack.iter().sum();
                    // let slack: Vec<f64> = slack.iter().map(|w| w / slack_sum).collect();

                    let v = h_mat.mat_vec(slack.values());
                    for c in 0..h_mat.cols() {
                        h_mat.col_mut(c).enumerate().for_each(|(i, hi)| *hi -= v[i]);
                    }
                }
                _ => {
                    unreachable!("txfr")
                }
            }
        }
    }

    Ok(h_mat)
}

/// Builds the line outage distribution factor matrix.
///
/// Returns the DC line outage distribution factor matrix for a given PTDF.
/// The matrix is `nbr x nbr`, where `nbr` is the number of branches.
pub fn make_lodf(branch: &Vec<Branch>, ptdf: Mat<f64>) -> Result<Mat<f64>, String> {
    let (nl, nb) = ptdf.shape();

    let mut c_ft = Coo::with_capacity(nb, nl, 2 * nl);
    for (i, br) in branch.iter().enumerate() {
        c_ft.push(br.f_bus, i, 1.0);
        c_ft.push(br.t_bus, i, -1.0);
    }

    /*
    H = PTDF * Cft;
    h = diag(H, 0);
    LODF = H ./ (ones(nl, nl) - ones(nl, 1) * h');
    LODF = LODF - diag(diag(LODF)) - eye(nl, nl);
    */

    let h_values = c_ft * ptdf.values();
    let h_mat = Mat::new(nl, nb).values(h_values).build()?;
    let h: Vec<f64> = h_mat.diagonal().collect();
    let lodf: Mat<f64> = h_mat
        / (Mat::new(nl, nl).ones().build()?
            - Mat::new(nl, 1)
                .ones()
                .build()?
                .mat_mat(&Mat::new(1, nl).values(h).build()?));
    let lodf: Mat<f64> = lodf
        - Mat::with_diagonal(lodf.diagonal().collect()).build()?
        - Mat::identity(nl).build()?;

    Ok(lodf)
}
