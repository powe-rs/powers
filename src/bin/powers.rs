use anyhow::Result;
use casecsv::Case;
use clap::{Args, Parser, Subcommand};
use powers7::{load_case, runpf, Alg, GenQLimits, MPOpt};
use spsolve::rlu::RLU;
use std::path::PathBuf;

/// Power flow simulation and optimization.
#[derive(Parser)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Power Flow
    #[clap(name = "pf")]
    PowerFlow(PfArgs),

    /// Optimal Power Flow
    OPF(OpfArgs),

    /// Continuation Power Flow
    CPF(CpfArgs),

    /// Power Transfer Distribution Factors
    PTDF(PtdfArgs),

    /// Line Outage Distribution Factors
    LODF(LodfArgs),
}

#[derive(Args)]
struct PfArgs {
    /// The input file
    #[arg(required = true)]
    input: PathBuf,

    /// Output file
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Linearized DC power flow.
    #[arg(long, default_value_t = false)]
    pub dc: bool,

    /// AC power flow algorithm.
    #[arg(long)]
    pub alg: Option<Alg>,

    /// Termination tolerance on per unit P & Q mismatch.
    #[arg(long)]
    pub tol: Option<f64>,

    /// Maximum number of iterations.
    #[arg(long)]
    pub max_it: Option<usize>,

    /// Enforce gen reactive power limits at expense of |V|.
    #[arg(long, default_value_t = false)]
    pub qlim: bool,
}

type OpfArgs = PfArgs;
type CpfArgs = PfArgs;
type PtdfArgs = PfArgs;
type LodfArgs = PfArgs;

fn main() {
    env_logger::Builder::from_default_env()
        .format_level(false)
        .format_target(false)
        .format_timestamp(None)
        .init();

    let cli = Cli::parse();

    match execute(&cli) {
        Ok(_) => {
            std::process::exit(0);
        }
        Err(err) => {
            eprintln!("error: {}", err);
            std::process::exit(2);
        }
    }
}

fn execute(cli: &Cli) -> Result<()> {
    let input = match &cli.command {
        Commands::PowerFlow(args) => &args.input,
        Commands::OPF(args) => &args.input,
        Commands::CPF(args) => &args.input,
        Commands::PTDF(args) => &args.input,
        Commands::LODF(args) => &args.input,
    };

    let mpc = load_case(input)?;

    let mut mpopt = MPOpt::default();
    if let Commands::PowerFlow(args) = &cli.command {
        mpopt.dc = args.dc;
        if let Some(alg) = args.alg {
            mpopt.pf.algorithm = alg;
        }
        if let Some(tol) = args.tol {
            mpopt.pf.tolerance = tol;
        }
        if let Some(max_it) = args.max_it {
            match mpopt.pf.algorithm {
                Alg::NR => {
                    mpopt.pf.max_it_nr = max_it;
                }
                Alg::FDBX | Alg::FDXB => {
                    mpopt.pf.max_it_fd = max_it;
                }
                Alg::GS => {
                    mpopt.pf.max_it_gs = max_it;
                }
                _ => {}
            }
        }
        if args.qlim {
            mpopt.pf.enforce_q_limits = GenQLimits::OneAtATime;
        }
    }

    let solver = RLU::default();

    let (mpc, success) = runpf(&mpc, &mpopt, &solver)?;
    if !success {
        return Err(anyhow::anyhow!("power flow did not succeed"));
    }

    let output = match &cli.command {
        Commands::PowerFlow(args) => &args.output,
        Commands::OPF(args) => &args.output,
        Commands::CPF(args) => &args.output,
        Commands::PTDF(args) => &args.output,
        Commands::LODF(args) => &args.output,
    };

    if let Some(out_path) = output {
        let case = Case::new(mpc.name).base_mva(mpc.base_mva).build()?;
        casecsv::write::write_zip(
            &out_path,
            &case,
            &mpc.bus,
            &mpc.gen,
            &mpc.branch,
            &Vec::default(),
            &Vec::default(),
            None,
        )?
    }

    Ok(())
}
