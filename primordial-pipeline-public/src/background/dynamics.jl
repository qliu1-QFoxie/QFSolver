struct _BackgroundParams{T<:AbstractFloat,P<:InflatonPotential}
    potential::P
    M_pl::T
    M_pl²::T
end

@inline _solver_algorithm(name::Symbol) =
    name === :vern9 ? Vern9() : throw(ArgumentError("BackgroundConfig.solver must be :vern9"))

@inline function _solver_kwargs(dt_min::T, dt_max::T) where {T<:AbstractFloat}
    has_dt_min = dt_min > zero(T)
    has_dt_max = dt_max > zero(T)
    has_dt_min && has_dt_max && return (; dtmin=dt_min, dtmax=dt_max)
    has_dt_min && return (; dtmin=dt_min)
    has_dt_max && return (; dtmax=dt_max)
    (;)
end

@inline function _state_scales(
    ϕ₀::T,
    ϕ_N₀::T,
    config::BackgroundConfig{T},
) where {T<:AbstractFloat}
    (
        ϕ=max(abs(ϕ₀), config.M_pl, one(T)),
        ϕ_N=max(abs(ϕ_N₀), config.M_pl, one(T)),
    )
end

@inline function _component_abstol(
    scales,
    config::BackgroundConfig{T},
) where {T<:AbstractFloat}
    T[
        config.abs_tol * scales.ϕ,
        config.abs_tol * scales.ϕ_N,
    ]
end

@inline _state_scale_max(scales) = max(scales.ϕ, scales.ϕ_N)

@inline function _scaled_maxnorm(u::Number, scales)
    abs(u) / _state_scale_max(scales)
end

@inline function _scaled_maxnorm(u, scales)
    max(abs(u[1]) / scales.ϕ, abs(u[2]) / scales.ϕ_N)
end

@inline _internalnorm(scales) = (u, N) -> _scaled_maxnorm(u, scales)

@inline _ε_H(ϕ_N, M_pl²) = ϕ_N^2 / (2 * M_pl²)

@inline function _ϕ_NN(p::_BackgroundParams{T}, ϕ, ϕ_N) where {T<:AbstractFloat}
    ε_H = _ε_H(ϕ_N, p.M_pl²)
    H² = _H²(p.potential, ϕ, ϕ_N, p.M_pl²)
    -(T(3) - ε_H) * ϕ_N - V_ϕ(p.potential, ϕ) / H²
end

@inline function _H²(potential::InflatonPotential, ϕ, ϕ_N, M_pl²)
    V(potential, ϕ) / (3 * M_pl² - ϕ_N^2 / 2)
end

@inline function _domain_violation(u, p::_BackgroundParams)
    ϕ = u[1]
    ϕ_N = u[2]
    !(isfinite(ϕ) && isfinite(ϕ_N)) && return true
    denom = 3 * p.M_pl² - ϕ_N^2 / 2
    denom <= zero(denom) && return true
    V(p.potential, ϕ) <= zero(denom) && return true
    false
end

function _validate_config(config::BackgroundConfig{T}) where {T<:AbstractFloat}
    isfinite(config.M_pl) && config.M_pl > zero(T) || throw(ArgumentError("M_pl must be finite and positive"))
    isfinite(config.N_start) || throw(ArgumentError("N_start must be finite"))
    isfinite(config.N_max) && config.N_max > zero(T) || throw(ArgumentError("N_max must be finite and positive"))
    isfinite(config.abs_tol) && config.abs_tol > zero(T) || throw(ArgumentError("abs_tol must be finite and positive"))
    isfinite(config.rel_tol) && config.rel_tol > zero(T) || throw(ArgumentError("rel_tol must be finite and positive"))
    config.max_iters > 0 || throw(ArgumentError("max_iters must be positive"))
    config.solver === :vern9 || throw(ArgumentError("solver must be :vern9"))
    isfinite(config.dt_min) && config.dt_min ≥ zero(T) || throw(ArgumentError("dt_min must be finite and non-negative"))
    isfinite(config.dt_max) && config.dt_max ≥ zero(T) || throw(ArgumentError("dt_max must be finite and non-negative"))
    (config.dt_max == zero(T) || config.dt_max ≥ config.dt_min) || throw(ArgumentError(
        "dt_max must be zero or greater than or equal to dt_min"))
    N_stop = config.N_start + config.N_max
    isfinite(N_stop) || throw(ArgumentError("N_start + N_max must be finite"))
    N_stop > config.N_start || throw(ArgumentError("N_start + N_max must exceed N_start"))
    config
end

function _validate_initial_state(
    potential::InflatonPotential,
    ϕ₀::T,
    ϕ_N₀::T,
    config::BackgroundConfig{T},
) where {T<:AbstractFloat}
    isfinite(ϕ₀) || throw(ArgumentError("ϕ₀ must be finite"))
    isfinite(ϕ_N₀) || throw(ArgumentError("ϕ_N₀ must be finite"))

    M_pl² = config.M_pl^2
    ε₀ = _ε_H(ϕ_N₀, M_pl²)
    ε₀ < one(T) || throw(ArgumentError("Initial ε_H must satisfy ε_H < 1"))

    denom = 3 * M_pl² - ϕ_N₀^2 / 2
    denom > zero(T) || throw(ArgumentError("Initial H² denominator must be positive"))

    V₀ = V(potential, ϕ₀)
    V₀ > zero(T) || throw(ArgumentError("Initial V(ϕ₀) must be positive"))
    nothing
end

function _background_rhs!(du, u, p::_BackgroundParams, N)
    ϕ = u[1]
    ϕ_N = u[2]

    du[1] = ϕ_N
    du[2] = _ϕ_NN(p, ϕ, ϕ_N)
    nothing
end

@inline function _inflation_end_residual(sol, p::_BackgroundParams{T}, N::T) where {T<:AbstractFloat}
    N_final = T(sol.t[end])
    tol = _N_tolerance(T, T(sol.t[1]), N_final)
    state = abs(N - N_final) <= tol ? sol(N_final; continuity=:right) : sol(N)
    _ε_H(state[2], p.M_pl²) - one(T)
end

@inline function _inflation_end_residual_derivative(
    sol,
    p::_BackgroundParams{T},
    N::T,
) where {T<:AbstractFloat}
    N_final = T(sol.t[end])
    tol = _N_tolerance(T, T(sol.t[1]), N_final)
    state = abs(N - N_final) <= tol ? sol(N_final; continuity=:right) : sol(N)
    ϕ_N = state[2]
    ϕ_N * _ϕ_NN(p, state[1], ϕ_N) / p.M_pl²
end

function _safeguarded_newton_root(
    f,
    f′,
    N_lo::T,
    N_hi::T;
    N_tol::T,
    f_tol::T,
    max_iters::Int=128,
) where {T<:AbstractFloat}
    left = N_lo
    right = N_hi
    f_left = f(left)
    f_right = f(right)

    abs(f_left) <= f_tol && return left
    abs(f_right) <= f_tol && return right
    signbit(f_left) == signbit(f_right) && throw(ArgumentError(
        "Safeguarded Newton requires the root to be bracketed"))

    x = abs(f_right) < abs(f_left) ? right : left
    f_x = abs(f_right) < abs(f_left) ? f_right : f_left

    for _ in 1:max_iters
        abs(f_x) <= f_tol && return x
        abs(right - left) <= N_tol && return x

        f′_x = f′(x)
        use_newton = isfinite(f′_x) && f′_x != zero(T)
        N_candidate = x

        if use_newton
            N_candidate = x - f_x / f′_x
            use_newton = left < N_candidate < right
        end

        use_newton || (N_candidate = (left + right) / T(2))
        f_candidate = f(N_candidate)
        isfinite(f_candidate) || throw(ArgumentError(
            "Safeguarded Newton encountered a non-finite function value inside the bracket"))

        if abs(f_candidate) <= f_tol
            return N_candidate
        end

        if signbit(f_left) == signbit(f_candidate)
            left = N_candidate
            f_left = f_candidate
        else
            right = N_candidate
            f_right = f_candidate
        end

        if abs(f_right) < abs(f_left)
            x = right
            f_x = f_right
        else
            x = left
            f_x = f_left
        end
    end

    abs(f_right) < abs(f_left) ? right : left
end

function _inflation_end_bracket(
    sol,
    p::_BackgroundParams{T};
    f_tol::T,
) where {T<:AbstractFloat}
    N_nodes = T.(sol.t)
    length(N_nodes) >= 2 || throw(ArgumentError(
        "Need at least two ODE nodes to bracket the end of inflation"))

    N_hi = N_nodes[end]
    f_hi = _inflation_end_residual(sol, p, N_hi)
    abs(f_hi) <= f_tol && return (N_nodes[end - 1], N_hi)
    f_hi > zero(T) || throw(ArgumentError(
        "Continuous callback terminated without placing the final solution near ε_H = 1"))

    for i in (length(N_nodes) - 1):-1:1
        N_i = N_nodes[i]
        f_i = _inflation_end_residual(sol, p, N_i)
        abs(f_i) <= f_tol && return (N_i, N_i)
        f_i < zero(T) && return (N_i, N_hi)
    end

    throw(ArgumentError("Could not locate a pre-event ODE node with ε_H < 1"))
end

function _polished_inflation_end(
    sol,
    p::_BackgroundParams{T},
    config::BackgroundConfig{T},
) where {T<:AbstractFloat}
    f_tol = _root_residual_tolerance(T, config.rel_tol)
    N_lo, N_hi = _inflation_end_bracket(sol, p; f_tol=f_tol)
    N_lo == N_hi && return N_lo

    N_tol = _root_N_tolerance(T, N_lo, N_hi, config.rel_tol)
    f(N) = _inflation_end_residual(sol, p, N)
    f′(N) = _inflation_end_residual_derivative(sol, p, N)
    _safeguarded_newton_root(f, f′, N_lo, N_hi; N_tol=N_tol, f_tol=f_tol)
end

function _termination_symbol(sol, N_end_requested, p::_BackgroundParams, T::Type)
    sol.retcode === ReturnCode.Terminated && return :inflation_end
    sol.retcode === ReturnCode.Success && return :N_max_reached
    sol.retcode === ReturnCode.Unstable && return :domain_violation

    final_state = sol.u[end]
    _domain_violation(final_state, p) && return :domain_violation

    N_reached = T(sol.t[end])
    N_reached < N_end_requested ? :solver_failure : :N_max_reached
end
