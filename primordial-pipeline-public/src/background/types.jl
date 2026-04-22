@inline _machine_eps(::Type{T}) where {T<:AbstractFloat} = eps(one(T))

@inline function _default_reltol(::Type{T}) where {T<:AbstractFloat}
    ϵ = _machine_eps(T)
    max(T(16) * ϵ, ϵ^(T(8) / T(9)))
end

@inline _default_abstol(::Type{T}) where {T<:AbstractFloat} = _default_reltol(T)

@inline function _N_tolerance(::Type{T}, N_lo, N_hi) where {T<:AbstractFloat}
    max(T(100) * eps(T) * max(abs(N_lo), abs(N_hi), one(T)), sqrt(eps(T)))
end

@inline function _root_N_tolerance(::Type{T}, N_lo, N_hi, rel_tol::T) where {T<:AbstractFloat}
    scale = max(abs(N_lo), abs(N_hi), one(T))
    max(T(64) * _machine_eps(T) * scale, rel_tol * scale)
end

@inline function _root_residual_tolerance(::Type{T}, rel_tol::T) where {T<:AbstractFloat}
    max(T(64) * _machine_eps(T), rel_tol)
end

Base.@kwdef struct BackgroundConfig{T<:AbstractFloat}
    M_pl::T = one(T)
    N_start::T = zero(T)
    N_max::T = T(200)
    abs_tol::T = _default_abstol(T)
    rel_tol::T = _default_reltol(T)
    max_iters::Int = 20_000_000
    solver::Symbol = :vern9
    dt_min::T = zero(T)
    dt_max::T = zero(T)
end

struct BackgroundSolution{T<:AbstractFloat,S,P<:InflatonPotential}
    ode_solution::S
    potential::P
    config::BackgroundConfig{T}
    N_end::T
    termination::Symbol
end

struct BackgroundPoint{T<:AbstractFloat}
    ϕ::T
    ϕ_N::T
    ϕ_NN::T
    V::T
    V_ϕ::T
    V_ϕϕ::T
    H²::T
    H::T
    ε_H::T
    η_H::T
    a::T
    aH::T
    μ_Q::T
end

struct τProfile{T<:AbstractFloat,B}
    background::B
    N_nodes::Vector{T}
    τ_nodes::Vector{T}
    τ_end::T
    abs_tol::T
    rel_tol::T
end
