module InflationPrimordialPipeline

using OrdinaryDiffEq: ODEProblem, solve, Vern9, ContinuousCallback, terminate!, ReturnCode
using QuadGK: quadgk
using SciMLBase: RightRootFind

include("numeric_validation.jl")
include("potentials.jl")

include("background/types.jl")
include("background/dynamics.jl")
include("background/observables.jl")
include("background/events.jl")
include("background/tau.jl")
include("background/solver.jl")

include("perturbations/types.jl")
include("perturbations/coefficients.jl")
include("perturbations/windows.jl")
include("perturbations/initial_conditions.jl")
include("perturbations/equations.jl")
include("perturbations/solver.jl")
include("perturbations/diagnostics.jl")

include("spectra/types.jl")
include("spectra/helpers.jl")
include("spectra/endpoint_formulas.jl")
include("spectra/diagnostics.jl")
include("spectra/solver.jl")
include("spectra/k_calibration.jl")
include("spectra/class_export.jl")

include("workflow/types.jl")
include("workflow/k_grid.jl")
include("workflow/solver.jl")
include("workflow/convenience.jl")

export AbstractPhysicalKRequest
export ExplicitPhysicalKIntervalRequest
export PivotCenteredPhysicalKRequest
export PrimordialPipelineRequest
export PrimordialPipelineResult
export build_exact_physical_spectrum_request
export solve_exact_physical_spectrum
export InflatonPotential, QuadraticPotential, V, V_ϕ, V_ϕϕ
export BackgroundConfig, BackgroundSolution, BackgroundPoint, τProfile
export solve_background, N_span, inflation_ended, termination_reason
export background_point, ϕ, ϕ_N, ϕ_NN, H², H, ε_H, η_H, a, aH, μ_Q
export N_hc, N_at_k_over_aH
export build_τ, τ, τ_N
export PerturbationKind, scalar, tensor
export ModeConfig, ModeCoefficients, AdiabaticMetrics, ModeInitializationDiagnostics, ModeWindow
export ExactSolveResult, ModeSolveFailure, ScalarModeSolution, TensorModeSolution
export mode_coefficients
export scalar_mode_friction, tensor_mode_friction
export scalar_mode_omega2, tensor_mode_omega2
export scalar_canonical_omega2, tensor_canonical_omega2
export scalar_canonical_frequency_sq, tensor_canonical_frequency_sq
export scalar_adiabatic_metrics, tensor_adiabatic_metrics
export mode_window
export scalar_bd_initial_conditions, tensor_bd_initial_conditions
export scalar_mode_rhs, tensor_mode_rhs
export solve_exact, solve_scalar_mode, solve_tensor_mode
export scalar_mode_state, scalar_mode_value, scalar_mode_derivative
export scalar_mode_canonical_state, scalar_mode_canonical_value, scalar_mode_canonical_derivative
export tensor_mode_state, tensor_mode_value, tensor_mode_derivative
export scalar_mode_wronskian, tensor_mode_wronskian
export scalar_mode_canonical_wronskian, tensor_mode_canonical_wronskian
export scalar_mode_wronskian_residual, tensor_mode_wronskian_residual
export ExactEndpointConfig, ExactEndpointDiagnostics
export ExactEndpointModeResult, ExactEndpointModeFailure, ExactEndpointCaseResult
export build_exact_logk_grid
export exact_scalar_endpoint_observable, exact_tensor_endpoint_observable
export exact_scalar_endpoint_log_slope, exact_tensor_endpoint_log_slope
export build_exact_guard_window, frozen_digits_from_rel_drift
export evaluate_exact_endpoint_diagnostics
export solve_exact_endpoint_mode, solve_exact_endpoint_case
export ExactKCalibrationConfig, ExactKCalibration, ExactKAxisMetadata
export ExactEndpointCaseOutput
export compute_k_calibration, internal_to_physical_k, physical_to_internal_k
export export_exact_endpoint_case
export ClassPrimordialExport
export build_class_primordial_export
export write_class_primordial_table, write_class_combined_table
export write_class_scalar_table, write_class_tensor_table, write_class_primordial_tables
export solve_primordial_pipeline

end
