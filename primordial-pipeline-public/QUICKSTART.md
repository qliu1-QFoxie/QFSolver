# Quickstart Command Manual

These commands assume Julia 1.12 is installed and available as `julia`.
Run them from the repository root, the directory containing `Project.toml`.

## 1. Install Julia dependencies

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

This command: 

- `--project=.` tells Julia to use this repository's `Project.toml`.
- `Pkg.instantiate()` installs the packages declared in `Project.toml`.
- No `Manifest.toml` is required. 

## 2. Check that the package loads

```bash
julia --project=. -e 'using InflationPrimordialPipeline; println("InflationPrimordialPipeline loaded")'
```

This command: 

- imports the package from the local repository;
- confirms that the dependency environment is usable.

## 3. Generate a CLASS external_Pk table

```bash
mkdir -p artifacts
julia --project=. examples/write_class_external_pk.jl \
  --output artifacts/primordial_external_pk.dat \
  --header \
  --k-min 0.03 \
  --k-max 0.08 \
  --ppd 16
```

This command: 

- creates an `artifacts/` output directory;
- solves a representative calibrated quadratic-potential case;
- writes a CLASS-compatible primordial table to
  `artifacts/primordial_external_pk.dat`;
- includes a commented header because `--header` is passed;
- requests physical wavenumbers from `0.03` to `0.08 Mpc^-1`;
- uses `16` log-grid points per decade.

The output table columns are:

```text
k_mpc P_s P_t
```

where `k_mpc` is in `Mpc^-1`, and `P_s` and `P_t` are dimensionless primordial
power spectra.

## 4. Inspect the generated table

```bash
head artifacts/primordial_external_pk.dat
```

This command: 

- prints the first rows of the generated primordial table;
- verifies that the file exists and contains numeric spectrum rows.

## 5. One-block minimal run

After entering the repository root, this block installs dependencies, checks the
package import, generates a small CLASS `external_Pk` table, and prints its first
rows:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using InflationPrimordialPipeline; println("InflationPrimordialPipeline loaded")'
mkdir -p artifacts
julia --project=. examples/write_class_external_pk.jl \
  --output artifacts/primordial_external_pk.dat \
  --header \
  --k-min 0.03 \
  --k-max 0.08 \
  --ppd 16
head artifacts/primordial_external_pk.dat
```

## 6. Solve a spectrum from the Julia API

```bash
julia --project=. -e '
using InflationPrimordialPipeline

potential = QuadraticPotential(6.0e-6)
request = build_exact_physical_spectrum_request(
    potential,
    16.0,
    -2.0 / 16.0;
    N_star = 25.0,
    k_request = ExplicitPhysicalKIntervalRequest(0.04, 0.06, 4.0),
    background_config = BackgroundConfig{Float64}(N_max = 80.0),
    mode_config = ModeConfig{Float64}(end_padding = 0.0),
    endpoint_config = ExactEndpointConfig{Float64}(),
)

result = solve_primordial_pipeline(request)
class_export = build_class_primordial_export(result)

println("points = ", length(class_export.k_mpc))
println("k_min = ", first(class_export.k_mpc))
println("k_max = ", last(class_export.k_mpc))
println("P_s_first = ", first(class_export.P_s))
println("P_t_first = ", first(class_export.P_t))
'
```

This command: 

- constructs the built-in quadratic potential `V(phi) = 0.5 m^2 phi^2`;
- converts the solved spectrum into the CLASS export payload;
- prints a small finite summary.

Meaning of the request entries:

- `potential`: the inflationary potential object. Here it is
  `QuadraticPotential(6.0e-6)`, representing `V(phi) = 0.5 m^2 phi^2` with
  `m = 6.0e-6` in reduced Planck units.
- `16.0`: the initial field value `phi0` at the start of the background solve.
- `-2.0 / 16.0`: the initial e-fold derivative `dphi/dN` at the start of the
  background solve. This is part of the model setup and is not inferred by the
  package.
- `N_star = 25.0`: the number of e-folds between the end of inflation and the
  epoch assigned to the physical pivot scale. This fixes the calibration from
  the internal comoving mode label to physical `k` in `Mpc^-1`.
- `k_request = ExplicitPhysicalKIntervalRequest(0.04, 0.06, 4.0)`: request a
  logarithmic physical-`k` grid from `0.04` to `0.06 Mpc^-1`, with `4.0` points
  per decade.
- `background_config = BackgroundConfig{Float64}(N_max = 80.0)`: solve the
  background using `Float64` arithmetic and allow the integration to run up to
  `N = 80` if inflation has not ended earlier.
- `mode_config = ModeConfig{Float64}(end_padding = 0.0)`: solve perturbation
  modes using `Float64` arithmetic and add no extra e-fold padding beyond the
  configured endpoint rule.
- `endpoint_config = ExactEndpointConfig{Float64}()`: use the default exact
  endpoint/freeze-out diagnostics for evaluating final scalar and tensor power.

The pipeline evolves the background, constructs conformal time, solves scalar
and tensor perturbation modes, applies the `N_star` physical-`k` calibration,
and returns dimensionless spectra `P_s(k)` and `P_t(k)`.

## 7. Use a manually defined potential

To solve a potential not built into the package, define a subtype of
`InflatonPotential` and provide three methods:

- `V(p, phi)` for the potential;
- `V_ϕ(p, phi)` for the first derivative;
- `V_ϕϕ(p, phi)` for the second derivative.

The example below defines a quartic potential,

```text
V(phi) = lambda * phi^4 / 4
```

and solves it on a small physical-`k` interval:

```bash
julia --project=. -e '
using InflationPrimordialPipeline

struct QuarticPotential{T<:AbstractFloat} <: InflatonPotential
    lambda::T
end

function InflationPrimordialPipeline.V(p::QuarticPotential, phi)
    p.lambda * phi^4 / 4
end

function InflationPrimordialPipeline.V_ϕ(p::QuarticPotential, phi)
    p.lambda * phi^3
end

function InflationPrimordialPipeline.V_ϕϕ(p::QuarticPotential, phi)
    3 * p.lambda * phi^2
end

potential = QuarticPotential(1.0e-13)
request = build_exact_physical_spectrum_request(
    potential,
    20.0,
    -0.10;
    N_star = 25.0,
    k_request = ExplicitPhysicalKIntervalRequest(0.04, 0.06, 4.0),
    background_config = BackgroundConfig{Float64}(N_max = 100.0),
    mode_config = ModeConfig{Float64}(end_padding = 0.0),
    endpoint_config = ExactEndpointConfig{Float64}(),
)

result = solve_primordial_pipeline(request)
class_export = build_class_primordial_export(result)

println("points = ", length(class_export.k_mpc))
println("first k [Mpc^-1] = ", first(class_export.k_mpc))
println("first P_s = ", first(class_export.P_s))
println("first P_t = ", first(class_export.P_t))
'
```

This command: 

- defines a new potential type (without editing package source code);
- gives the pipeline the potential and its first two field derivatives;
- uses the same public physical-`k` workflow as the built-in example;
- returns calibrated scalar and tensor primordial spectra.

The custom-potential request entries have the same meaning as in the built-in
quadratic example:

- `potential = QuarticPotential(1.0e-13)`: use the manually defined quartic
  potential with coupling `lambda = 1.0e-13`.
- `20.0`: initial field value `phi0`.
- `-0.10`: initial e-fold derivative `dphi/dN`.
- `N_star = 25.0`: physical-`k` calibration choice.
- `ExplicitPhysicalKIntervalRequest(0.04, 0.06, 4.0)`: physical-`k` range and
  grid density.
- `BackgroundConfig{Float64}(N_max = 100.0)`: background integration settings.
- `ModeConfig{Float64}(end_padding = 0.0)`: perturbation-mode integration
  settings.
- `ExactEndpointConfig{Float64}()`: endpoint spectrum-evaluation settings.

The chosen initial field value, field velocity, and `N_star` are part of the
physical model setup. The package does not infer them automatically.

## 8. Another custom potential: quadratic with a step

This example shows a manually entered potential with several parameters:

```text
V(phi) = 0.5 m^2 phi^2 [1 + c tanh((phi - phi_s) / d)]
```

The derivatives are entered explicitly, just as for the quartic example:

```bash
julia --project=. -e '
using InflationPrimordialPipeline

struct StepQuadraticPotential{T<:AbstractFloat} <: InflatonPotential
    m::T
    c::T
    d::T
    phi_s::T
end

function InflationPrimordialPipeline.V(p::StepQuadraticPotential, phi)
    step = tanh((phi - p.phi_s) / p.d)
    0.5 * p.m^2 * phi^2 * (1 + p.c * step)
end

function InflationPrimordialPipeline.V_ϕ(p::StepQuadraticPotential, phi)
    x = (phi - p.phi_s) / p.d
    step = tanh(x)
    sech2 = 1 / cosh(x)^2
    0.5 * p.m^2 * (2 * phi * (1 + p.c * step) + phi^2 * p.c * sech2 / p.d)
end

function InflationPrimordialPipeline.V_ϕϕ(p::StepQuadraticPotential, phi)
    x = (phi - p.phi_s) / p.d
    step = tanh(x)
    sech2 = 1 / cosh(x)^2
    0.5 * p.m^2 * (
        2 * (1 + p.c * step) +
        4 * phi * p.c * sech2 / p.d -
        2 * phi^2 * p.c * sech2 * step / p.d^2
    )
end

potential = StepQuadraticPotential(6.0e-6, 0.0015, 0.03, 14.3)
request = build_exact_physical_spectrum_request(
    potential,
    18.5,
    -0.10810810810810811;
    N_star = 55.0,
    k_request = ExplicitPhysicalKIntervalRequest(0.04, 0.06, 4.0),
    background_config = BackgroundConfig{Float64}(N_max = 140.0),
    mode_config = ModeConfig{Float64}(end_padding = 0.0),
    endpoint_config = ExactEndpointConfig{Float64}(),
)

result = solve_primordial_pipeline(request)
class_export = build_class_primordial_export(result)

println("points = ", length(class_export.k_mpc))
println("first k [Mpc^-1] = ", first(class_export.k_mpc))
println("first P_s = ", first(class_export.P_s))
println("first P_t = ", first(class_export.P_t))
'
```

This adds beyond the quartic example:

- the potential type stores four parameters: `m`, `c`, `d`, and `phi_s`;
- the first and second derivatives include the derivative of the `tanh` feature;
- the initial field value and `N_max` are chosen for this featured model rather
  than reused from the simple quadratic example.

## 9. Optional CLASS use

CLASS is not installed by this package. To use the generated primordial table in
CLASS, configure CLASS separately with its `external_Pk` primordial mode and
point its command to the generated table, for example:

```text
P_k_ini type = external_Pk
command = cat artifacts/primordial_external_pk.dat
```
