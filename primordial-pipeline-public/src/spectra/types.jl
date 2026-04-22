Base.@kwdef struct ExactEndpointConfig{T<:AbstractFloat}
    endpoint_offset::T = zero(T)
    requested_guard_width::T = one(T)
    guard_samples::Int = 257
end

Base.@kwdef struct ExactEndpointDiagnostics{T<:AbstractFloat}
    requested_guard_width::T
    achieved_guard_width::T
    guard_samples::Int
    diagnostics_complete::Bool
    diagnostics_failure_reason::Union{Nothing,String}
    guard_start_N::T
    guard_end_N::T
    endpoint_log_slope::T
    endpoint_rel_drift::T
    max_guard_rel_drift::T
    max_guard_abs_log_slope::T
    frozen_digits::T
end

Base.@kwdef struct ExactEndpointModeResult{K,N,V,W,D,I}
    kind::PerturbationKind
    k::K
    window::W
    N_eval::N
    observable_name::Symbol
    observable_value::V
    diagnostics::D
    initialization::I
end

Base.@kwdef struct ExactEndpointModeFailure{K,I}
    kind::PerturbationKind
    k::K
    initialization::I
    failure_stage::Symbol
    failure_reason::String
end

Base.@kwdef struct ExactEndpointCaseResult{K,S,T}
    k_grid::Vector{K}
    scalar_modes::Vector{S}
    tensor_modes::Vector{T}
end

_validate_guard_samples(guard_samples::Integer) =
    guard_samples > 1 ? Int(guard_samples) : throw(ArgumentError("guard_samples must be > 1"))
