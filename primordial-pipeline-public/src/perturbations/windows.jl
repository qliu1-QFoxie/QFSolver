@inline function _validate_mode_window_config(config::ModeConfig{T}) where {T<:AbstractFloat}
    isfinite(config.subhorizon_ratio) && config.subhorizon_ratio > one(T) || throw(ArgumentError(
        "subhorizon_ratio must be finite and strictly greater than 1"))
    isfinite(config.end_padding) && config.end_padding >= zero(T) || throw(ArgumentError(
        "end_padding must be finite and non-negative"))
    config.max_iters > 0 || throw(ArgumentError("max_iters must be positive"))
    isfinite(config.abs_tol) && config.abs_tol > zero(T) || throw(ArgumentError(
        "abs_tol must be finite and strictly positive"))
    isfinite(config.rel_tol) && config.rel_tol > zero(T) || throw(ArgumentError(
        "rel_tol must be finite and strictly positive"))
    isfinite(config.adiabatic_alpha1_tol) && config.adiabatic_alpha1_tol > zero(T) || throw(ArgumentError(
        "adiabatic_alpha1_tol must be finite and strictly positive"))
    isfinite(config.adiabatic_alpha2_tol) && config.adiabatic_alpha2_tol > zero(T) || throw(ArgumentError(
        "adiabatic_alpha2_tol must be finite and strictly positive"))
    isfinite(config.adiabatic_derivative_step) && config.adiabatic_derivative_step > zero(T) || throw(ArgumentError(
        "adiabatic_derivative_step must be finite and strictly positive"))
    isfinite(config.adiabatic_search_step) && config.adiabatic_search_step > zero(T) || throw(ArgumentError(
        "adiabatic_search_step must be finite and strictly positive"))
    nothing
end

@inline function _mode_stop(bg::BackgroundSolution{T}, end_padding::T) where {T<:AbstractFloat}
    N_lo, N_hi = N_span(bg)
    N_stop = N_hi - end_padding
    N_stop > N_lo || throw(ArgumentError(
        "end_padding=$(end_padding) leaves no room inside the solved background interval [$N_lo, $N_hi]"))
    N_stop
end

Base.@kwdef struct _ModeWindowStatus{T<:AbstractFloat}
    window::Union{Nothing,ModeWindow{T}}
    N_in_ratio::T
    N_in_used::T
    init_failure_reason::Union{Nothing,String}
end

@inline function _shared_adiabaticity_pass(
    bg::BackgroundSolution{T},
    k::T,
    N::T,
    config::ModeConfig{T},
) where {T<:AbstractFloat}
    scalar_metrics = scalar_adiabatic_metrics(
        bg,
        k,
        N;
        h=config.adiabatic_derivative_step,
        alpha1_tol=config.adiabatic_alpha1_tol,
        alpha2_tol=config.adiabatic_alpha2_tol,
    )
    tensor_metrics = tensor_adiabatic_metrics(
        bg,
        k,
        N;
        h=config.adiabatic_derivative_step,
        alpha1_tol=config.adiabatic_alpha1_tol,
        alpha2_tol=config.adiabatic_alpha2_tol,
    )
    scalar_metrics.passed && tensor_metrics.passed
end

@inline function _adiabatic_refine_tolerance(
    ::Type{T},
    N_lo::T,
    N_hi::T,
    config::ModeConfig{T},
) where {T<:AbstractFloat}
    max(
        sqrt(eps(T)) * max(abs(N_lo), abs(N_hi), one(T)),
        config.adiabatic_search_step / T(32),
    )
end

function _latest_adiabatic_start(
    bg::BackgroundSolution{T},
    k::T,
    N_ratio::T,
    config::ModeConfig{T},
) where {T<:AbstractFloat}
    N_lo, N_hi = N_span(bg)
    _shared_adiabaticity_pass(bg, k, N_ratio, config) && return (N_ratio, nothing)

    last_failing = N_ratio
    current = N_ratio
    first_passing = nothing

    while current > N_lo
        next_trial = max(N_lo, current - config.adiabatic_search_step)
        next_trial < current || break
        if _shared_adiabaticity_pass(bg, k, next_trial, config)
            first_passing = next_trial
            break
        end
        current = next_trial
        last_failing = current
    end

    if first_passing === nothing
        failure_reason = "no exact adiabatic plane-wave BD start exists in the solved background interval before N_in_ratio=$(N_ratio)"
        return (current, failure_reason)
    end

    if !config.adiabatic_refine
        return (first_passing, nothing)
    end

    left = first_passing
    right = last_failing
    tolerance = _adiabatic_refine_tolerance(T, N_lo, N_hi, config)
    while right - left > tolerance
        midpoint = (left + right) / T(2)
        if _shared_adiabaticity_pass(bg, k, midpoint, config)
            left = midpoint
        else
            right = midpoint
        end
    end

    (left, nothing)
end

"""
The ratio-based start `k/(aH)=subhorizon_ratio` provides the initial candidate.
When `init_requires_adiabatic=true`, the start point is moved earlier until the
exact canonical scalar and tensor adiabaticity tests both pass.
"""
function _mode_window_status(
    bg::BackgroundSolution{T},
    k::Real;
    config::ModeConfig{T}=ModeConfig{T}(),
) where {T<:AbstractFloat}
    _validate_mode_window_config(config)

    k_value = _checked_positive_finite(k, "k", T)
    N_ratio = N_at_k_over_aH(bg, k_value, config.subhorizon_ratio)
    N_cross = N_hc(bg, k_value)
    N_stop = _mode_stop(bg, config.end_padding)

    N_used, failure_reason = if config.init_requires_adiabatic
        _latest_adiabatic_start(bg, k_value, N_ratio, config)
    else
        (N_ratio, nothing)
    end

    if failure_reason !== nothing
        return _ModeWindowStatus{T}(
            window=nothing,
            N_in_ratio=N_ratio,
            N_in_used=N_used,
            init_failure_reason=failure_reason,
        )
    end

    N_used < N_cross || throw(ArgumentError(
        "subhorizon_ratio=$(config.subhorizon_ratio) must place mode initialization before horizon crossing"))
    N_stop > N_cross || throw(ArgumentError(
        "end_padding=$(config.end_padding) leaves no room after horizon crossing before the background endpoint"))

    _ModeWindowStatus{T}(
        window=ModeWindow{T}(N_used, N_cross, N_stop),
        N_in_ratio=N_ratio,
        N_in_used=N_used,
        init_failure_reason=nothing,
    )
end

function mode_window(
    bg::BackgroundSolution{T},
    k::Real;
    config::ModeConfig{T}=ModeConfig{T}(),
) where {T<:AbstractFloat}
    status = _mode_window_status(bg, k; config=config)
    status.window === nothing && throw(ArgumentError(status.init_failure_reason))
    status.window
end
