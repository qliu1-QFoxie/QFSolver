@inline _perturbation_machine_eps(::Type{T}) where {T<:AbstractFloat} = eps(one(T))

@inline function _default_mode_reltol(::Type{T}) where {T<:AbstractFloat}
    epsilon = _perturbation_machine_eps(T)
    max(T(16) * epsilon, epsilon^(T(8) / T(9)))
end

@inline _default_mode_abstol(::Type{T}) where {T<:AbstractFloat} = _default_mode_reltol(T)

@inline function _default_end_padding(::Type{T}) where {T<:AbstractFloat}
    max(T(0.3), T(32) * sqrt(_perturbation_machine_eps(T)))
end

@inline _default_adiabatic_alpha1_tol(::Type{T}) where {T<:AbstractFloat} = T(1e-3)
@inline _default_adiabatic_alpha2_tol(::Type{T}) where {T<:AbstractFloat} = T(1e-2)
@inline _default_adiabatic_derivative_step(::Type{T}) where {T<:AbstractFloat} = T(1e-3)
@inline _default_adiabatic_search_step(::Type{T}) where {T<:AbstractFloat} = T(0.25)

@enum PerturbationKind scalar tensor

"""
    ModeConfig{T}

When `init_requires_adiabatic=true`, plane-wave BD initial conditions are
accepted only after an exact adiabaticity check based on the canonical
scalar/tensor frequencies.
"""
Base.@kwdef struct ModeConfig{T<:AbstractFloat}
    subhorizon_ratio::T = T(80)
    end_padding::T = _default_end_padding(T)
    abs_tol::T = _default_mode_abstol(T)
    rel_tol::T = _default_mode_reltol(T)
    max_iters::Int = 20_000_000
    solver::Symbol = :vern9
    dt_min::T = zero(T)
    dt_max::T = zero(T)
    init_requires_adiabatic::Bool = true
    adiabatic_alpha1_tol::T = _default_adiabatic_alpha1_tol(T)
    adiabatic_alpha2_tol::T = _default_adiabatic_alpha2_tol(T)
    adiabatic_derivative_step::T = _default_adiabatic_derivative_step(T)
    adiabatic_search_step::T = _default_adiabatic_search_step(T)
    adiabatic_refine::Bool = true
end

struct ModeCoefficients{T<:AbstractFloat}
    epsilon_H::T
    aH::T
    mu_Q::T
end

Base.@kwdef struct AdiabaticMetrics{T<:AbstractFloat}
    omega_sq::T
    alpha1::T
    alpha2::T
    passed::Bool
end

Base.@kwdef struct ModeInitializationDiagnostics{T<:AbstractFloat}
    N_in_ratio::T
    N_in_used::T
    q_in::T
    init_omega_sq::T
    init_alpha1::T
    init_alpha2::T
    init_adiabatic_passed::Bool
    init_kind::String
    init_failure_reason::Union{Nothing,String}
end

struct ModeWindow{T<:AbstractFloat}
    N_in::T
    N_hc::T
    N_stop::T
end

Base.@kwdef struct ExactSolveResult{T<:AbstractFloat,S,I}
    kind::PerturbationKind
    precision_type::Type{T} = T
    config::ModeConfig{T}
    window::ModeWindow{T}
    solution::S
    initialization::I
end

function ExactSolveResult(
    config::ModeConfig{T},
    solution,
    initialization,
) where {T<:AbstractFloat}
    kind = solution isa ScalarModeSolution ? scalar :
        solution isa TensorModeSolution ? tensor :
        throw(ArgumentError(
            "ExactSolveResult requires a ScalarModeSolution or TensorModeSolution"))
    ExactSolveResult{T,typeof(solution),typeof(initialization)}(
        kind=kind,
        config=config,
        window=solution.window,
        solution=solution,
        initialization=initialization,
    )
end

struct ModeSolveFailure{K,W,R} <: Exception
    kind::PerturbationKind
    k::K
    window::W
    retcode::R
end

function Base.showerror(io::IO, err::ModeSolveFailure)
    print(
        io,
        "Mode solve failed with retcode $(err.retcode) for kind=$(err.kind) at k=$(err.k) " *
        "over N ∈ [$(err.window.N_in), $(err.window.N_stop)]",
    )
end
