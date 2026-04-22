@inline function _log_aH(bg::BackgroundSolution{T}, N::T) where {T<:AbstractFloat}
    state = bg.ode_solution(N)
    H²_value = _H²(bg.potential, state[1], state[2], bg.config.M_pl^2)
    N + log(H²_value) / T(2)
end

@inline function _log_crossing_residual(
    bg::BackgroundSolution{T},
    log_target::T,
    N::T,
) where {T<:AbstractFloat}
    log_target - _log_aH(bg, N)
end

@inline function _log_crossing_residual_derivative(
    bg::BackgroundSolution{T},
    N::T,
) where {T<:AbstractFloat}
    state = bg.ode_solution(N)
    _ε_H(state[2], bg.config.M_pl^2) - one(T)
end

function _solve_log_crossing(
    bg::BackgroundSolution{T},
    log_target::T;
    lower_error::AbstractString,
    upper_error::AbstractString,
) where {T<:AbstractFloat}
    N_lo, N_hi = N_span(bg)
    f(N) = _log_crossing_residual(bg, log_target, N)
    f_lo = f(N_lo)
    f_hi = f(N_hi)
    f_tol = _root_residual_tolerance(T, bg.config.rel_tol)

    f_lo < -f_tol && throw(ArgumentError(lower_error))
    f_hi > f_tol && throw(ArgumentError(upper_error))
    abs(f_lo) <= f_tol && return N_lo
    abs(f_hi) <= f_tol && return N_hi

    N_tol = _root_N_tolerance(T, N_lo, N_hi, bg.config.rel_tol)
    f′(N) = _log_crossing_residual_derivative(bg, N)
    _safeguarded_newton_root(f, f′, N_lo, N_hi; N_tol=N_tol, f_tol=f_tol)
end

function N_hc(bg::BackgroundSolution{T}, k::Real) where {T<:AbstractFloat}
    k_value = _checked_positive_finite(k, "k", T)
    log_target = log(k_value)
    N_lo, N_hi = N_span(bg)

    _solve_log_crossing(
        bg,
        log_target;
        lower_error=
            "k=$(k) is already superhorizon at N=$(N_lo); choose an earlier background start",
        upper_error=
            "k=$(k) does not cross the horizon before N=$(N_hi); extend the background interval",
    )
end

function N_at_k_over_aH(bg::BackgroundSolution{T}, k::Real, r::Real) where {T<:AbstractFloat}
    k_value = _checked_positive_finite(k, "k", T)
    r_value = _checked_positive_finite(r, "r", T)
    log_target = log(k_value) - log(r_value)
    N_lo, N_hi = N_span(bg)

    _solve_log_crossing(
        bg,
        log_target;
        lower_error=
            "k/(aH) is already below r=$(r) at N=$(N_lo); choose an earlier background start",
        upper_error=
            "k/(aH) does not reach r=$(r) before N=$(N_hi); extend the background interval",
    )
end
