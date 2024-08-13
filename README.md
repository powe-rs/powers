# powers

Power flow simulation and optimization.

## About

This crate contains common structures and auxiliary functions.
The case structure is defined in the [caseformat](https://crates.io/crates/powers-pf) crate.
The main functions of `powers` are provided in dependent crates:

- [`powers-pf`](https://crates.io/crates/powers-pf) Power Flow
    - Newton's method
    - Fast-Decoupled method
    - Gauss-Seidel method
    - DC power flow

Additional features are available from the author on request:

- `powers-opf` Optimal Power Flow
- `powers-cpf` Continuation Power Flow
- `powers-ptdf` Power Transfer Distribution Factors/Line Outage Distribution Factors

## License

Translated from [MATPOWER](https://matpower.org/) into [Rust](https://www.rust-lang.org/) by Richard Lincoln.

The source code is distributed under the same BSD 3-clause license ([LICENSE](LICENSE) or
https://opensource.org/license/bsd-3-clause/) as MATPOWER.