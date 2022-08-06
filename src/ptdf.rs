use crate::dc::make_b_dc;
use crate::mpc::{Branch, BusType, MPC};
use crate::traits::LinearSolverN;
use densetools::arr::Arr;
use densetools::mat::Mat;
use densetools::slice::{find, ne, not, select, sum};
use sparsetools::coo::Coo;
use sparsetools::csr::CSR;

#[derive(Clone)]
enum PTDFSlack {
    /// Single slack bus. The reference bus is used by default.
    Single(usize),
    /// Weights specifying the proportion of the slack taken up at each bus.
    Weights(Vec<f64>),
    /// Each column specifies how the slack should be handled for injections
    /// at that bus. This option only applies when computing the full PTDF matrix.
    Matrix(Mat<f64>),
}

#[derive(Clone)]
enum PTDFBus {
    Txfr(Mat<f64>),
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
    linsol: &dyn LinearSolverN,
) -> Result<Mat<f64>, String> {
    if slack.is_none() {
        for (i, bus) in mpc.bus.iter().enumerate() {
            if bus.bus_type == BusType::REF {
                slack = Some(PTDFSlack::Single(i));
                break;
            }
        }
        if slack.is_none() {
            return Err("slack must not be none".to_string());
        }
    }

    // compute full PTDF?
    let mut compute_full_H = true;
    let nb = mpc.bus.len();
    let nbr = mpc.branch.len();
    let mut txfr = false; // default assumes not just for specific transfers
    let mut dP: Option<Mat<f64>> = None;
    if let Some(ptdf_bus) = bus_opt {
        compute_full_H = false;
        if let PTDFBus::Txfr(txfr_mat) = ptdf_bus {
            if txfr_mat.rows() != nb {
                return Err(format!(
                    "txfr: matrix must have nb rows; nb = {} rows = {}",
                    txfr_mat.rows(),
                    nb
                ));
            }
            for r in 0..txfr_mat.rows() {
                let r_sum = txfr_mat.row(r).sum();
                if r_sum != 0.0 {
                    return Err(format!(
                        "txfr: columns of each row must sum to zero; row = {} sum = {}",
                        r, r_sum
                    )
                    .to_string());
                }
            }
            txfr = true; // cols of H for specific (slack independent) transfers
            dP = Some(txfr_mat);
        }
    }

    // set the slack bus to be used to compute initial PTDF
    let mut slack_bus = 0; // use bus 1 for temp slack bus
    if let Some(ptdf_slack) = slack {
        if let PTDFSlack::Single(slack_bus_index) = ptdf_slack {
            slack_bus = slack_bus_index
        }
    }

    let noref = Arr::arange(1, nb, 1); // use bus 1 for voltage angle reference
    let noslack = find(&ne(&Arr::range(nb), slack_bus));

    // check that bus numbers are equal to indices to bus (one set of bus numbers)
    for (i, bus) in mpc.bus.iter().enumerate() {
        if bus.i != i {
            return Err(
                "make_ptdf: buses must be numbered consecutively in bus matrix".to_string(),
            );
        }
    }

    //  compute PTDF for single slack_bus  //
    let (Bbus, Bf, _, _) = make_b_dc(mpc.base_mva, &mpc.bus, &mpc.branch);

    // set up injections/transfers
    let (nbi0, nbi, bidx) = if compute_full_H {
        // full H for all columns
        let nbi = nb; // number of buses of interest (i.e. all of them)
        dP = Some(Mat::identity(nb)); // TODO: sparse RHS
        (0, nbi, None)
    } else {
        if txfr {
            // for specific transfers
            let nbi = dP.unwrap().cols(); // number of transfers
            (0, nbi, None)
        } else {
            let bus_idx = bus_opt.unwrap().bus_idx();

            // for subset of columns of H for specific buses
            // if necessary add missing slacks to bus index list
            let nbi0 = bus_idx.len(); // number of original buses of interest
            let mut bidx: Vec<usize> = bus_idx.clone();
            if let Some(ptdf_slack) = slack {
                if let PTDFSlack::Weights(weights) = ptdf_slack {
                    let slacks = find::<f64, usize>(&weights); // get all slack buses
                    for sw in slacks {
                        if !bus_idx.contains(&sw) {
                            bidx.push(sw);
                        }
                    }
                }
            }
            let nbi = bidx.len(); // number of buses of interest

            // define the equivalent transfer, each column is a single transfer
            //dP = accumarray([bidx (1:nbi)'], 1, [nb, nbi]);
            // let mut dP = Coo::empty(nb, nbi, 2 * nbi);
            let mut dp = Mat::zeros(nb, nbi);
            for (i, j) in bidx.iter().enumerate() {
                // TODO
                // dp.push(i, j, 1.0);
            }
            dP = Some(dp);
            (nbi0, nbi, Some(bidx))
        }
    };

    // solve for change in voltage angles
    let bbus_noslack = Bbus.select(Some(&noslack), Some(&noref))?;
    let dp_noslack = dP.unwrap().select_rows(&noslack);

    let d_theta_noref = linsol.solve(bbus_noslack.to_csc(), &dp_noslack)?;

    // let d_theta = Mat::zeros(nb, nbi);
    // let d_theta = d_theta_noref;
    let mut d_theta = Mat::zeros(1, nbi);
    d_theta.data_mut().extend(d_theta_noref.data());

    // compute corresponding change in branch flows
    let mut H = Bf * d_theta;

    // distribute slack, if requested //
    if !txfr {
        if let Some(ptdf_slack) = slack {
            // if let PTDFSlack::Weights(weights) = ptdf_slack {
            // }
            match ptdf_slack {
                PTDFSlack::Weights(weights) => {
                    // slack is a vector of weights

                    let slack = slack / slack.sum(); // normalize weights

                    // conceptually, we want to do ...
                    //    H = H * (eye(nb,nb) - slack * ones(1, nb));
                    // ... we just do it more efficiently

                    if compute_full_H {
                        let v = H * slack;
                        for k in 0..nb {
                            H.col_mut(k).sub(v);
                        }
                    } else {
                        let v = H * select(&slack, bidx.unwrap());
                        for k in 0..nbi {
                            H.col_mut(k).sub(v);
                        }
                        H = H.select(None, Some(0..nbi0)); // remove temp cols added for missing slacks
                    }
                }
                PTDFSlack::Matrix(slack) => {
                    let slack = slack / sum(&slack); // normalize weights
                    H = H - H * slack;
                }
                _ => {}
            }
        }
    }

    Ok(H)
}

/// Builds the line outage distribution factor matrix.
///
/// Returns the DC line outage distribution factor matrix for a given PTDF.
/// The matrix is `nbr x nbr`, where `nbr` is the number of branches.
pub fn make_lodf(branch: &Vec<Branch>, ptdf: Mat<f64>) -> Result<Mat<f64>, String> {
    let (nl, nb) = ptdf.shape();

    let mut c_ft = Coo::empty(nb, nl, 2 * nl);
    for (i, br) in branch.iter().enumerate() {
        c_ft.push(br.from_bus, i, 1.0);
        c_ft.push(br.to_bus, i, -1.0);
    }

    /*
    H = PTDF * Cft;
    h = diag(H, 0);
    LODF = H ./ (ones(nl, nl) - ones(nl, 1) * h');
    LODF = LODF - diag(diag(LODF)) - eye(nl, nl);
    */

    let h_mat = ptdf * c_ft;
    let h: Arr<f64> = h_mat.diagonal();
    let lodf = h_mat / (Mat::ones(nl, nl) - Mat::ones(nl, 1).mat_mat(&Mat::new(1, nl, h.vec())));
    let lodf: Mat<f64> = lodf - Mat::with_diag(lodf.diag()) - Mat::eye(nl, nl);

    Ok(lodf)
}
