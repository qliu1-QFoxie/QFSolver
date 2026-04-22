@inline function _diagnostic_tolerance(::Type{T}, values...) where {T<:AbstractFloat}
    max(T(100) * eps(T) * max((abs(v) for v in values)..., one(T)), sqrt(eps(T)))
end

function _validate_endpoint_config(config::ExactEndpointConfig{T}) where {T<:AbstractFloat}
    _checked_nonnegative_finite(config.endpoint_offset, "endpoint_offset", T)
    _checked_positive_finite(config.requested_guard_width, "requested_guard_width", T)
    _validate_guard_samples(config.guard_samples)
    nothing
end

function build_exact_guard_window(
    window::ModeWindow{T},
    N_eval::Real,
    endpoint_config::ExactEndpointConfig{T},
) where {T<:AbstractFloat}
    _validate_endpoint_config(endpoint_config)
    N_eval_value = T(N_eval)
    requested_width = endpoint_config.requested_guard_width
    guard_start = max(N_eval_value - requested_width, window.N_hc)
    achieved_width = N_eval_value - guard_start
    achieved_width >= zero(T) || throw(ArgumentError(
        "achieved guard width must be non-negative"))

    tol = _diagnostic_tolerance(T, N_eval_value, guard_start, requested_width)
    diagnostics_complete = achieved_width > zero(T) && achieved_width >= requested_width - tol

    (
        guard_start_N=guard_start,
        guard_end_N=N_eval_value,
        requested_guard_width=requested_width,
        achieved_guard_width=achieved_width,
        diagnostics_complete=diagnostics_complete,
    )
end

@inline function frozen_digits_from_rel_drift(rel_drift::T) where {T<:AbstractFloat}
    isfinite(rel_drift) || return T(NaN)
    rel_drift < zero(T) && throw(ArgumentError("relative drift must be non-negative"))
    iszero(rel_drift) && return T(Inf)
    -log10(rel_drift)
end

function _diagnostics_nan_bundle(::Type{T}) where {T<:AbstractFloat}
    nan_value = T(NaN)
    (
        endpoint_log_slope=nan_value,
        endpoint_rel_drift=nan_value,
        max_guard_rel_drift=nan_value,
        max_guard_abs_log_slope=nan_value,
        frozen_digits=nan_value,
    )
end

function _scalar_observable_and_slope(
    solution::ScalarModeSolution{T},
    N::T,
) where {T<:AbstractFloat}
    (
        exact_scalar_endpoint_observable(solution, N),
        exact_scalar_endpoint_log_slope(solution, N),
    )
end

function _tensor_observable_and_slope(
    solution::TensorModeSolution{T},
    N::T,
) where {T<:AbstractFloat}
    (
        exact_tensor_endpoint_observable(solution, N),
        exact_tensor_endpoint_log_slope(solution, N),
    )
end

@inline function _exact_endpoint_log_slope(
    exact_result::ExactSolveResult{T},
    N::T,
) where {T<:AbstractFloat}
    if exact_result.kind == scalar
        return exact_scalar_endpoint_log_slope(exact_result.solution, N)
    elseif exact_result.kind == tensor
        return exact_tensor_endpoint_log_slope(exact_result.solution, N)
    end

    throw(ArgumentError("exact endpoint diagnostics requires scalar or tensor kind"))
end

@inline function _exact_endpoint_observable(
    exact_result::ExactSolveResult{T},
    N::T,
) where {T<:AbstractFloat}
    if exact_result.kind == scalar
        return exact_scalar_endpoint_observable(exact_result.solution, N)
    elseif exact_result.kind == tensor
        return exact_tensor_endpoint_observable(exact_result.solution, N)
    end

    throw(ArgumentError("exact endpoint diagnostics requires scalar or tensor kind"))
end

@inline function _exact_endpoint_observable_and_slope(
    exact_result::ExactSolveResult{T},
    N::T,
) where {T<:AbstractFloat}
    if exact_result.kind == scalar
        return _scalar_observable_and_slope(exact_result.solution, N)
    elseif exact_result.kind == tensor
        return _tensor_observable_and_slope(exact_result.solution, N)
    end

    throw(ArgumentError("exact endpoint diagnostics requires scalar or tensor kind"))
end

@inline _expected_diagnostics_failure(err) = err isa ArgumentError || err isa DomainError

function _incomplete_exact_endpoint_diagnostics(
    ::Type{T},
    guard;
    endpoint_log_slope::T,
    endpoint_config::ExactEndpointConfig{T},
    failure_reason::Union{Nothing,String},
) where {T<:AbstractFloat}
    undefined_metrics = _diagnostics_nan_bundle(T)
    ExactEndpointDiagnostics(
        requested_guard_width=guard.requested_guard_width,
        achieved_guard_width=guard.achieved_guard_width,
        guard_samples=endpoint_config.guard_samples,
        diagnostics_complete=false,
        diagnostics_failure_reason=failure_reason,
        guard_start_N=guard.guard_start_N,
        guard_end_N=guard.guard_end_N,
        endpoint_log_slope=endpoint_log_slope,
        endpoint_rel_drift=undefined_metrics.endpoint_rel_drift,
        max_guard_rel_drift=undefined_metrics.max_guard_rel_drift,
        max_guard_abs_log_slope=undefined_metrics.max_guard_abs_log_slope,
        frozen_digits=undefined_metrics.frozen_digits,
    )
end

function _evaluate_exact_endpoint_diagnostics_or_incomplete(
    exact_result::ExactSolveResult{T},
    endpoint_config::ExactEndpointConfig{T};
    endpoint_log_slope_evaluator=_exact_endpoint_log_slope,
    endpoint_observable_evaluator=_exact_endpoint_observable,
    sample_evaluator=_exact_endpoint_observable_and_slope,
) where {T<:AbstractFloat}
    _validate_endpoint_config(endpoint_config)

    window = exact_result.window
    N_eval = window.N_stop
    guard = build_exact_guard_window(window, N_eval, endpoint_config)

    endpoint_log_slope = try
        endpoint_log_slope_evaluator(exact_result, N_eval)
    catch err
        _expected_diagnostics_failure(err) || rethrow()
        return _incomplete_exact_endpoint_diagnostics(
            T,
            guard;
            endpoint_log_slope=_diagnostics_nan_bundle(T).endpoint_log_slope,
            endpoint_config=endpoint_config,
            failure_reason=sprint(showerror, err),
        )
    end

    if iszero(guard.achieved_guard_width)
        return _incomplete_exact_endpoint_diagnostics(
            T,
            guard;
            endpoint_log_slope=endpoint_log_slope,
            endpoint_config=endpoint_config,
            failure_reason=nothing,
        )
    end

    endpoint_value = try
        endpoint_observable_evaluator(exact_result, N_eval)
    catch err
        _expected_diagnostics_failure(err) || rethrow()
        return _incomplete_exact_endpoint_diagnostics(
            T,
            guard;
            endpoint_log_slope=endpoint_log_slope,
            endpoint_config=endpoint_config,
            failure_reason=sprint(showerror, err),
        )
    end

    sample_grid = collect(range(guard.guard_start_N, guard.guard_end_N; length=endpoint_config.guard_samples))
    max_rel_drift = zero(T)
    max_abs_log_slope = zero(T)
    first_rel_drift = zero(T)

    for (index, N_sample) in enumerate(sample_grid)
        observable_value, abs_log_slope = try
            sample_evaluator(exact_result, N_sample)
        catch err
            _expected_diagnostics_failure(err) || rethrow()
            return _incomplete_exact_endpoint_diagnostics(
                T,
                guard;
                endpoint_log_slope=endpoint_log_slope,
                endpoint_config=endpoint_config,
                failure_reason=sprint(showerror, err),
            )
        end

        rel_drift = abs(observable_value - endpoint_value) / endpoint_value
        if !(isfinite(rel_drift) && isfinite(abs_log_slope))
            return _incomplete_exact_endpoint_diagnostics(
                T,
                guard;
                endpoint_log_slope=endpoint_log_slope,
                endpoint_config=endpoint_config,
                failure_reason="sampled diagnostics must remain finite across the guard window",
            )
        end

        if index == 1
            first_rel_drift = rel_drift
        end

        max_rel_drift = max(max_rel_drift, rel_drift)
        max_abs_log_slope = max(max_abs_log_slope, abs_log_slope)
    end

    ExactEndpointDiagnostics(
        requested_guard_width=guard.requested_guard_width,
        achieved_guard_width=guard.achieved_guard_width,
        guard_samples=endpoint_config.guard_samples,
        diagnostics_complete=guard.diagnostics_complete,
        diagnostics_failure_reason=nothing,
        guard_start_N=guard.guard_start_N,
        guard_end_N=guard.guard_end_N,
        endpoint_log_slope=endpoint_log_slope,
        endpoint_rel_drift=first_rel_drift,
        max_guard_rel_drift=max_rel_drift,
        max_guard_abs_log_slope=max_abs_log_slope,
        frozen_digits=frozen_digits_from_rel_drift(max_rel_drift),
    )
end

function evaluate_exact_endpoint_diagnostics(
    exact_result::ExactSolveResult{T},
    endpoint_config::ExactEndpointConfig{T},
) where {T<:AbstractFloat}
    _evaluate_exact_endpoint_diagnostics_or_incomplete(exact_result, endpoint_config)
end
