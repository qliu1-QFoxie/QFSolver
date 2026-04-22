function _validate_physical_k_request(req::ExplicitPhysicalKIntervalRequest)
    _validate_requested_k_interval(req.k_min_phys_mpc, req.k_max_phys_mpc)
    _validate_dense_ppd(req.ppd)
    nothing
end

function _validate_physical_k_request(req::PivotCenteredPhysicalKRequest)
    _checked_positive_finite(
        req.half_width_log10,
        "half_width_log10",
        typeof(req.half_width_log10),
    )
    _validate_dense_ppd(req.ppd)
    nothing
end

@inline function _resolved_exact_request_endpoint_config(::Type{T}, endpoint_config) where {T<:AbstractFloat}
    _resolved_endpoint_config(T, endpoint_config)
end

@inline function _resolved_mode_config(::Type{T}, mode_config) where {T<:AbstractFloat}
    mode_config isa ModeConfig{T} || throw(ArgumentError(
        "mode_config must be ModeConfig{$T}"))
    mode_config
end

@inline function _resolved_endpoint_config(::Type{T}, endpoint_config) where {T<:AbstractFloat}
    if endpoint_config === nothing
        return ExactEndpointConfig{T}()
    end
    endpoint_config isa ExactEndpointConfig{T} || throw(ArgumentError(
        "endpoint_config must be nothing or ExactEndpointConfig{$T}"))
    endpoint_config
end

@inline function _resolved_workflow_k_calibration_config(::Type{T}, calibration_config) where {T<:AbstractFloat}
    config = _resolved_k_calibration_config(T, calibration_config)
    _validate_k_calibration_config(config)
    config.calibrate_k_axis || throw(ArgumentError(
        "physical-k pipeline requests require calibrate_k_axis=true"))
    config.export_physical_k || throw(ArgumentError(
        "physical-k pipeline requests require export_physical_k=true"))
    config
end

function _resolved_exact_request_mode_config(
    ::Type{T},
    mode_config,
    endpoint_config::ExactEndpointConfig{T},
) where {T<:AbstractFloat}
    if mode_config === nothing
        return ModeConfig{T}(
            subhorizon_ratio=ModeConfig{T}().subhorizon_ratio,
            end_padding=endpoint_config.endpoint_offset,
            abs_tol=ModeConfig{T}().abs_tol,
            rel_tol=ModeConfig{T}().rel_tol,
            max_iters=ModeConfig{T}().max_iters,
            solver=ModeConfig{T}().solver,
            dt_min=ModeConfig{T}().dt_min,
            dt_max=ModeConfig{T}().dt_max,
            init_requires_adiabatic=ModeConfig{T}().init_requires_adiabatic,
            adiabatic_alpha1_tol=ModeConfig{T}().adiabatic_alpha1_tol,
            adiabatic_alpha2_tol=ModeConfig{T}().adiabatic_alpha2_tol,
            adiabatic_derivative_step=ModeConfig{T}().adiabatic_derivative_step,
            adiabatic_search_step=ModeConfig{T}().adiabatic_search_step,
            adiabatic_refine=ModeConfig{T}().adiabatic_refine,
        )
    end

    mode_config isa ModeConfig{T} || throw(ArgumentError(
        "mode_config must be nothing or ModeConfig{$T}"))
    mode_config.end_padding == endpoint_config.endpoint_offset || throw(ArgumentError(
        "mode_config.end_padding must match endpoint_config.endpoint_offset in the exact primordial workflow"))
    mode_config
end

function _resolved_exact_request_k_request(
    ::Type{T},
    k_request,
    half_width_log10::Real,
    ppd::Real,
) where {T<:AbstractFloat}
    if k_request === nothing
        return PivotCenteredPhysicalKRequest(
            T(half_width_log10),
            Float64(ppd),
        )
    end

    k_request isa AbstractPhysicalKRequest || throw(ArgumentError(
        "k_request must be nothing or an AbstractPhysicalKRequest"))
    k_request
end

function _resolved_exact_request_calibration_config(
    ::Type{T},
    calibration_config,
    N_star::Real,
    k_pivot_phys_mpc::Real,
) where {T<:AbstractFloat}
    if calibration_config === nothing
        config = ExactKCalibrationConfig{T}(
            calibrate_k_axis=true,
            N_star=T(N_star),
            k_pivot_phys_mpc=T(k_pivot_phys_mpc),
            anchor_convention=:hubble_crossing,
            export_physical_k=true,
        )
        _validate_k_calibration_config(config)
        return config
    end

    _resolved_workflow_k_calibration_config(T, calibration_config)
end

function _resolve_requested_k_grids(
    calibration::ExactKCalibration{T},
    req::ExplicitPhysicalKIntervalRequest,
) where {T<:AbstractFloat}
    _validate_physical_k_request(req)
    k_min_phys = _checked_positive_finite(req.k_min_phys_mpc, "k_min_phys_mpc", T)
    k_max_phys = _checked_positive_finite(req.k_max_phys_mpc, "k_max_phys_mpc", T)
    physical_k_grid = build_exact_logk_grid(k_min_phys, k_max_phys; ppd=req.ppd)
    internal_k_grid = T[physical_to_internal_k(k_phys, calibration) for k_phys in physical_k_grid]
    _validate_requested_k_grid(internal_k_grid)
    return physical_k_grid, internal_k_grid
end

function _resolve_requested_k_grids(
    calibration::ExactKCalibration{T},
    req::PivotCenteredPhysicalKRequest,
) where {T<:AbstractFloat}
    _validate_physical_k_request(req)
    k_pivot_phys = calibration.k_pivot_phys_mpc
    half_width_log10 = _checked_positive_finite(req.half_width_log10, "half_width_log10", T)
    half_factor = exp(log(T(10)) * half_width_log10)
    k_min_phys = k_pivot_phys / half_factor
    k_max_phys = k_pivot_phys * half_factor
    physical_k_grid = build_exact_logk_grid(k_min_phys, k_max_phys; ppd=req.ppd)
    internal_k_grid = T[physical_to_internal_k(k_phys, calibration) for k_phys in physical_k_grid]
    _validate_requested_k_grid(internal_k_grid)
    return physical_k_grid, internal_k_grid
end
