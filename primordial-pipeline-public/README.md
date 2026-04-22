# InflationPrimordialPipeline

Julia tools for solving canonical single-field inflationary backgrounds,
evolving exact scalar and tensor perturbations, calibrating the primordial
spectrum to a physical `k` axis, and exporting CLASS-compatible external
primordial tables.

## Scope

This package computes primordial scalar and tensor spectra for a specified
inflaton potential, homogeneous initial data, and explicit `N_star` calibration.
It does not determine reheating physics, infer `N_star`, or run CLASS internally.

## Setup

The required Julia packages are declared in `Project.toml`. A `Manifest.toml`
is not required for normal use; Julia can resolve compatible package versions
from the project file.

From the repository root, install the declared dependencies with:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Then check that the package loads:

```bash
julia --project=. -e 'using InflationPrimordialPipeline'
```

For a copy-paste command manual, see [`QUICKSTART.md`](QUICKSTART.md).

## Example

Generate a CLASS `external_Pk` table from the representative quadratic example:

```bash
julia --project=. examples/write_class_external_pk.jl --output artifacts/primordial_external_pk.dat
```

The output table is whitespace-delimited with columns:

```text
k_mpc P_s P_t
```

where `k_mpc` is in `Mpc^-1`, and `P_s` and `P_t` are dimensionless primordial
power spectra.

## Optional CLASS Use

CLASS is not a dependency of this package. To use the generated table with
CLASS, configure CLASS separately with its `external_Pk` primordial-input mode.
