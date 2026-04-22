@inline function _mode_window_status_or_failure(
    bg::BackgroundSolution{T},
    k::Real;
    config::ModeConfig{T},
) where {T<:AbstractFloat}
    try
        _mode_window_status(bg, k; config=config)
    catch err
        err isa ArgumentError || rethrow()
        N_lo, _ = N_span(bg)
        _ModeWindowStatus{T}(
            window=nothing,
            N_in_ratio=N_lo,
            N_in_used=N_lo,
            init_failure_reason=sprint(showerror, err),
        )
    end
end

@inline function _resolved_tau_profile(
    bg::BackgroundSolution,
    tau_profile,
)
    if tau_profile === nothing
        return build_τ(bg)
    end
    if hasproperty(tau_profile, :background) && getproperty(tau_profile, :background) !== bg
        throw(ArgumentError("tau_profile must belong to the supplied background"))
    end
    tau_profile
end

@inline function _resolved_exact_solve_config(
    mode_config::ModeConfig{T},
    endpoint_config::ExactEndpointConfig{T},
) where {T<:AbstractFloat}
    ModeConfig{T}(
        subhorizon_ratio=mode_config.subhorizon_ratio,
        end_padding=endpoint_config.endpoint_offset,
        abs_tol=mode_config.abs_tol,
        rel_tol=mode_config.rel_tol,
        max_iters=mode_config.max_iters,
        solver=mode_config.solver,
        dt_min=mode_config.dt_min,
        dt_max=mode_config.dt_max,
        init_requires_adiabatic=mode_config.init_requires_adiabatic,
        adiabatic_alpha1_tol=mode_config.adiabatic_alpha1_tol,
        adiabatic_alpha2_tol=mode_config.adiabatic_alpha2_tol,
        adiabatic_derivative_step=mode_config.adiabatic_derivative_step,
        adiabatic_search_step=mode_config.adiabatic_search_step,
        adiabatic_refine=mode_config.adiabatic_refine,
    )
end

@inline _failure_reason(err) = sprint(showerror, err)

@inline _expected_mode_solve_failure(err) = err isa ModeSolveFailure

@inline _expected_endpoint_observable_failure(err) = err isa ArgumentError || err isa DomainError

function _solve_exact_result_with_window(
    kind::PerturbationKind,
    bg::BackgroundSolution{T},
    k_value::T,
    window::ModeWindow{T},
    solve_config::ModeConfig{T};
    tau_profile=nothing,
) where {T<:AbstractFloat}
    solve_exact(
        kind,
        bg,
        k_value;
        precision=T,
        config=solve_config,
        window=window,
        tau_profile=tau_profile,
    )
end

function _exact_endpoint_observable(exact_result::ExactSolveResult)
    if exact_result.kind == scalar
        return :P_s, exact_scalar_endpoint_observable(
            exact_result.solution,
            exact_result.window.N_stop,
        )
    elseif exact_result.kind == tensor
        return :P_t, exact_tensor_endpoint_observable(
            exact_result.solution,
            exact_result.window.N_stop,
        )
    end

    error("exact endpoint observable requires scalar or tensor kind")
end

function _solve_exact_endpoint_mode_with_window(
    kind::PerturbationKind,
    bg::BackgroundSolution{T},
    k_value::T,
    window::ModeWindow{T},
    solve_config::ModeConfig{T},
    endpoint_config::ExactEndpointConfig{T};
    tau_profile=nothing,
) where {T<:AbstractFloat}
    exact_result = _solve_exact_result_with_window(
        kind,
        bg,
        k_value,
        window,
        solve_config;
        tau_profile=tau_profile,
    )
    observable_name, observable_value = _exact_endpoint_observable(exact_result)

    ExactEndpointModeResult(
        kind=kind,
        k=k_value,
        window=exact_result.window,
        N_eval=exact_result.window.N_stop,
        observable_name=observable_name,
        observable_value=observable_value,
        diagnostics=evaluate_exact_endpoint_diagnostics(exact_result, endpoint_config),
        initialization=exact_result.initialization,
    )
end

function _mode_failure_result(
    kind::PerturbationKind,
    k_value,
    initialization;
    failure_stage::Symbol,
    failure_reason,
)
    ExactEndpointModeFailure(
        kind=kind,
        k=k_value,
        initialization=initialization,
        failure_stage=failure_stage,
        failure_reason=String(failure_reason),
    )
end

function solve_exact_endpoint_mode(
    kind::PerturbationKind,
    bg::BackgroundSolution{T},
    k;
    mode_config::ModeConfig{T}=ModeConfig{T}(),
    endpoint_config::ExactEndpointConfig{T}=ExactEndpointConfig{T}(),
    tau_profile=nothing,
) where {T<:AbstractFloat}
    _validate_endpoint_config(endpoint_config)

    k_value = _checked_positive_finite(k, "k", T)
    solve_config = _resolved_exact_solve_config(mode_config, endpoint_config)
    window = mode_window(bg, k_value; config=solve_config)
    _solve_exact_endpoint_mode_with_window(
        kind,
        bg,
        k_value,
        window,
        solve_config,
        endpoint_config;
        tau_profile=tau_profile,
    )
end

function _mode_failure_result(
    kind::PerturbationKind,
    bg::BackgroundSolution{T},
    k_value::T,
    status::_ModeWindowStatus{T},
    config::ModeConfig{T},
    ;
    failure_stage::Symbol=:window,
    failure_reason=status.init_failure_reason,
) where {T<:AbstractFloat}
    _mode_failure_result(
        kind,
        k_value,
        _mode_initialization_diagnostics(
            kind,
            bg,
            k_value,
            status.N_in_used,
            config;
            failure_reason=failure_reason,
            N_in_ratio=status.N_in_ratio,
        ),
        ;
        failure_stage=failure_stage,
        failure_reason=something(failure_reason, "unknown mode failure"),
    )
end

function _solve_exact_endpoint_mode_or_failure(
    kind::PerturbationKind,
    bg::BackgroundSolution{T},
    k_value::T,
    status::_ModeWindowStatus{T},
    solve_config::ModeConfig{T},
    endpoint_config::ExactEndpointConfig{T};
    tau_profile=nothing,
) where {T<:AbstractFloat}
    window = status.window
    window === nothing && return _mode_failure_result(
        kind,
        bg,
        k_value,
        status,
        solve_config;
        failure_stage=:window,
    )

    exact_result = try
        _solve_exact_result_with_window(
            kind,
            bg,
            k_value,
            window,
            solve_config;
            tau_profile=tau_profile,
        )
    catch err
        _expected_mode_solve_failure(err) || rethrow()
        return _mode_failure_result(
            kind,
            bg,
            k_value,
            status,
            solve_config;
            failure_stage=:solve,
            failure_reason=_failure_reason(err),
        )
    end

    observable_name, observable_value = try
        _exact_endpoint_observable(exact_result)
    catch err
        _expected_endpoint_observable_failure(err) || rethrow()
        return _mode_failure_result(
            kind,
            k_value,
            exact_result.initialization;
            failure_stage=:endpoint_observable,
            failure_reason=_failure_reason(err),
        )
    end

    ExactEndpointModeResult(
        kind=kind,
        k=k_value,
        window=exact_result.window,
        N_eval=exact_result.window.N_stop,
        observable_name=observable_name,
        observable_value=observable_value,
        diagnostics=evaluate_exact_endpoint_diagnostics(exact_result, endpoint_config),
        initialization=exact_result.initialization,
    )
end

function solve_exact_endpoint_case(
    bg::BackgroundSolution{T},
    k_grid;
    mode_config::ModeConfig{T}=ModeConfig{T}(),
    endpoint_config::ExactEndpointConfig{T}=ExactEndpointConfig{T}(),
    tau_profile=nothing,
) where {T<:AbstractFloat}
    _validate_endpoint_config(endpoint_config)
    k_values = T[_checked_positive_finite(k, "k", T) for k in k_grid]
    _validate_requested_k_grid(k_values)

    tau_profile_value = _resolved_tau_profile(bg, tau_profile)
    solve_config = _resolved_exact_solve_config(mode_config, endpoint_config)

    ScalarModeResult = Union{ExactEndpointModeResult,ExactEndpointModeFailure}
    TensorModeResult = Union{ExactEndpointModeResult,ExactEndpointModeFailure}
    scalar_modes = Vector{ScalarModeResult}(undef, length(k_values))
    tensor_modes = Vector{TensorModeResult}(undef, length(k_values))

    for index in eachindex(k_values)
        k_value = k_values[index]
        status = _mode_window_status_or_failure(bg, k_value; config=solve_config)
        scalar_modes[index] = _solve_exact_endpoint_mode_or_failure(
            scalar,
            bg,
            k_value,
            status,
            solve_config,
            endpoint_config;
            tau_profile=tau_profile_value,
        )
        tensor_modes[index] = _solve_exact_endpoint_mode_or_failure(
            tensor,
            bg,
            k_value,
            status,
            solve_config,
            endpoint_config;
            tau_profile=tau_profile_value,
        )
    end

    _validate_case_alignment(k_values, scalar_modes, tensor_modes)
    ExactEndpointCaseResult(k_grid=k_values, scalar_modes=scalar_modes, tensor_modes=tensor_modes)
end
