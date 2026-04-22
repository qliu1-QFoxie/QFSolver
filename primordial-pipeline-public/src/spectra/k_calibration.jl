"""
    ExactKCalibrationConfig{T}

Configuration for externally anchoring the unified pipeline `k` axis to
physical units.

Internal `k` in this pipeline is defined relative to the internal scale-factor
normalization `a = exp(N)`. The inflationary solve fixes the primordial spectra
as functions of this internal `k`, but only up to one overall multiplicative
normalization of the horizontal axis.

`N_star` is therefore treated as an external anchoring input, not as a quantity
predicted by the inflation ODE alone. This implementation uses only the
standard Hubble-crossing pivot convention

`k_pivot_internal = a(N_end - N_star) * H(N_end - N_star)`.

In USR or other non-attractor regimes, `N_star` only chooses which internal
mode is identified with the physical pivot mode. It does **not** determine when
to read the amplitude. Late-time amplitudes remain those already produced by the
exact endpoint solver.
"""
Base.@kwdef struct ExactKCalibrationConfig{T<:AbstractFloat}
    calibrate_k_axis::Bool = false
    N_star::Union{Nothing,T} = nothing
    k_pivot_phys_mpc::T = T(0.05)
    anchor_convention::Symbol = :hubble_crossing
    export_physical_k::Bool = false
end

Base.@kwdef struct ExactKCalibration{T<:AbstractFloat}
    enabled::Bool
    N_star::Union{Nothing,T}
    N_end::T
    N_pivot::Union{Nothing,T}
    k_pivot_phys_mpc::Union{Nothing,T}
    k_pivot_internal::Union{Nothing,T}
    C_k::Union{Nothing,T}
    anchor_convention::Union{Nothing,Symbol}
end

Base.@kwdef struct ExactKAxisMetadata{T<:AbstractFloat}
    k_axis_type::String
    k_pivot_phys_mpc::Union{Nothing,T} = nothing
    N_star::Union{Nothing,T} = nothing
    N_pivot::Union{Nothing,T} = nothing
    k_pivot_internal::Union{Nothing,T} = nothing
    C_k::Union{Nothing,T} = nothing
    anchor_convention::Union{Nothing,String} = nothing
end

Base.@kwdef struct ExactEndpointCaseOutput{K,C,M}
    spectrum_case::C
    k_grid::Vector{K}
    metadata::M
end

@inline function _resolved_k_calibration_config(::Type{T}, calibration_config) where {T<:AbstractFloat}
    if calibration_config === nothing
        return ExactKCalibrationConfig{T}()
    end

    calibration_config isa ExactKCalibrationConfig || throw(ArgumentError(
        "calibration_config must be nothing or an ExactKCalibrationConfig"))

    ExactKCalibrationConfig{T}(
        calibrate_k_axis=calibration_config.calibrate_k_axis,
        N_star=calibration_config.N_star === nothing ? nothing : T(calibration_config.N_star),
        k_pivot_phys_mpc=T(calibration_config.k_pivot_phys_mpc),
        anchor_convention=calibration_config.anchor_convention,
        export_physical_k=calibration_config.export_physical_k,
    )
end

@inline function _k_calibration_domain_tolerance(::Type{T}, values...) where {T<:AbstractFloat}
    max(T(100) * eps(T) * max((abs(v) for v in values)..., one(T)), sqrt(eps(T)))
end

function _validate_k_calibration_config(config::ExactKCalibrationConfig{T}) where {T<:AbstractFloat}
    _checked_positive_finite(config.k_pivot_phys_mpc, "k_pivot_phys_mpc", T)

    if !config.calibrate_k_axis
        config.export_physical_k && throw(ArgumentError(
            "export_physical_k=true requires calibrate_k_axis=true"))
        return nothing
    end

    config.N_star === nothing && throw(ArgumentError(
        "calibrate_k_axis=true requires N_star to be provided"))
    _checked_positive_finite(config.N_star, "N_star", T)
    config.anchor_convention == :hubble_crossing || throw(ArgumentError(
        "unsupported anchor_convention=$(config.anchor_convention); only :hubble_crossing is supported"))
    config.export_physical_k || throw(ArgumentError(
        "physical-k pipeline requests require export_physical_k=true"))

    nothing
end

function compute_k_calibration(
    bg::BackgroundSolution{T},
    calibration_config=ExactKCalibrationConfig{T}(),
) where {T<:AbstractFloat}
    config = _resolved_k_calibration_config(T, calibration_config)
    _validate_k_calibration_config(config)

    if !config.calibrate_k_axis
        return ExactKCalibration{T}(
            enabled=false,
            N_star=nothing,
            N_end=bg.N_end,
            N_pivot=nothing,
            k_pivot_phys_mpc=nothing,
            k_pivot_internal=nothing,
            C_k=nothing,
            anchor_convention=nothing,
        )
    end

    N_lo, N_hi = N_span(bg)
    N_pivot = bg.N_end - config.N_star
    tol = _k_calibration_domain_tolerance(T, N_lo, N_hi, bg.N_end, N_pivot)

    if N_pivot < N_lo - tol || N_pivot > N_hi + tol
        throw(ArgumentError(
            "N_star=$(config.N_star) places N_pivot=$(N_pivot) outside the solved background interval [$N_lo, $N_hi]"))
    end

    N_pivot_value = clamp(N_pivot, N_lo, N_hi)
    k_pivot_internal = aH(bg, N_pivot_value)
    C_k = config.k_pivot_phys_mpc / k_pivot_internal
    isfinite(C_k) && C_k > zero(T) || throw(ArgumentError(
        "computed k-axis calibration constant must be finite and strictly positive"))

    ExactKCalibration{T}(
        enabled=true,
        N_star=config.N_star,
        N_end=bg.N_end,
        N_pivot=N_pivot_value,
        k_pivot_phys_mpc=config.k_pivot_phys_mpc,
        k_pivot_internal=k_pivot_internal,
        C_k=C_k,
        anchor_convention=config.anchor_convention,
    )
end

@inline function _require_enabled_k_calibration(calibration::ExactKCalibration)
    calibration.enabled || throw(ArgumentError(
        "this operation requires an enabled k calibration"))
    calibration
end

function internal_to_physical_k(
    k_internal::Real,
    calibration::ExactKCalibration{T},
) where {T<:AbstractFloat}
    calibration_value = _require_enabled_k_calibration(calibration)
    _checked_positive_finite(k_internal, "k_internal", T) * calibration_value.C_k
end

function physical_to_internal_k(
    k_phys::Real,
    calibration::ExactKCalibration{T},
) where {T<:AbstractFloat}
    calibration_value = _require_enabled_k_calibration(calibration)
    _checked_positive_finite(k_phys, "k_phys", T) / calibration_value.C_k
end

function _build_k_axis_metadata(calibration::ExactKCalibration{T}) where {T<:AbstractFloat}
    if !calibration.enabled
        return ExactKAxisMetadata{T}(k_axis_type="internal")
    end

    ExactKAxisMetadata{T}(
        k_axis_type="physical_mpc^-1",
        k_pivot_phys_mpc=calibration.k_pivot_phys_mpc,
        N_star=calibration.N_star,
        N_pivot=calibration.N_pivot,
        k_pivot_internal=calibration.k_pivot_internal,
        C_k=calibration.C_k,
        anchor_convention=calibration.anchor_convention === nothing ? nothing : String(calibration.anchor_convention),
    )
end

@inline function _k_axis_match_tolerance(::Type{T}, values...) where {T<:AbstractFloat}
    max(T(100) * eps(T) * max((abs(v) for v in values)..., one(T)), sqrt(eps(T)))
end

function _validate_calibration_matches_background(
    calibration::ExactKCalibration{T},
    bg::BackgroundSolution{T},
) where {T<:AbstractFloat}
    N_end_tol = _k_axis_match_tolerance(T, calibration.N_end, bg.N_end)
    abs(calibration.N_end - bg.N_end) <= N_end_tol || throw(ArgumentError(
        "resolved k calibration does not match the supplied background end time"))

    if !calibration.enabled
        return calibration
    end

    calibration.N_pivot === nothing && throw(ArgumentError(
        "enabled k calibration must store N_pivot"))
    calibration.k_pivot_internal === nothing && throw(ArgumentError(
        "enabled k calibration must store k_pivot_internal"))
    calibration.C_k === nothing && throw(ArgumentError(
        "enabled k calibration must store C_k"))
    calibration.k_pivot_phys_mpc === nothing && throw(ArgumentError(
        "enabled k calibration must store k_pivot_phys_mpc"))

    N_lo, N_hi = N_span(bg)
    pivot_tol = _k_axis_match_tolerance(T, N_lo, N_hi, calibration.N_pivot)
    if calibration.N_pivot < N_lo - pivot_tol || calibration.N_pivot > N_hi + pivot_tol
        throw(ArgumentError(
            "resolved k calibration pivot lies outside the supplied background interval"))
    end

    expected_k_pivot_internal = aH(bg, calibration.N_pivot)
    k_pivot_tol = _k_axis_match_tolerance(T, expected_k_pivot_internal, calibration.k_pivot_internal)
    abs(expected_k_pivot_internal - calibration.k_pivot_internal) <= k_pivot_tol || throw(ArgumentError(
        "resolved k calibration does not match the supplied background pivot mode"))

    expected_C_k = calibration.k_pivot_phys_mpc / expected_k_pivot_internal
    C_k_tol = _k_axis_match_tolerance(T, expected_C_k, calibration.C_k)
    abs(expected_C_k - calibration.C_k) <= C_k_tol || throw(ArgumentError(
        "resolved k calibration does not match the supplied background calibration constant"))

    calibration
end

function _expected_exported_physical_k_grid(
    case_result::ExactEndpointCaseResult,
    calibration::ExactKCalibration{T},
) where {T<:AbstractFloat}
    T[internal_to_physical_k(k, calibration) for k in case_result.k_grid]
end

function _validated_requested_physical_k_grid(
    case_result::ExactEndpointCaseResult,
    calibration::ExactKCalibration{T},
    requested_physical_k_grid,
) where {T<:AbstractFloat}
    calibration.enabled || throw(ArgumentError(
        "requested_physical_k_grid requires an enabled k calibration"))

    k_values = T[_checked_positive_finite(k, "requested_physical_k_grid", T) for k in requested_physical_k_grid]
    length(k_values) == length(case_result.k_grid) || throw(ArgumentError(
        "requested_physical_k_grid length must match solved exact spectrum length"))
    _validate_requested_k_grid(k_values)

    expected_k_values = _expected_exported_physical_k_grid(case_result, calibration)
    for index in eachindex(k_values)
        provided = k_values[index]
        expected = expected_k_values[index]
        tol = _k_axis_match_tolerance(T, provided, expected)
        abs(provided - expected) <= tol || throw(ArgumentError(
            "requested_physical_k_grid[$index]=$(provided) does not match the calibrated solved mode label $(expected)"))
    end

    k_values
end

function _export_exact_endpoint_case_with_calibration(
    case_result::ExactEndpointCaseResult,
    calibration::ExactKCalibration{T};
    requested_physical_k_grid=nothing,
) where {T<:AbstractFloat}
    exported_k_grid = if requested_physical_k_grid === nothing
        if calibration.enabled
            _expected_exported_physical_k_grid(case_result, calibration)
        else
            copy(case_result.k_grid)
        end
    else
        _validated_requested_physical_k_grid(
            case_result,
            calibration,
            requested_physical_k_grid,
        )
    end

    ExactEndpointCaseOutput(
        spectrum_case=case_result,
        k_grid=exported_k_grid,
        metadata=_build_k_axis_metadata(calibration),
    )
end

function export_exact_endpoint_case(
    case_result::ExactEndpointCaseResult,
    bg::BackgroundSolution{T};
    calibration_config=ExactKCalibrationConfig{T}(),
    requested_physical_k_grid=nothing,
) where {T<:AbstractFloat}
    calibration = compute_k_calibration(bg, calibration_config)
    _validate_calibration_matches_background(
        calibration,
        bg,
    )
    _export_exact_endpoint_case_with_calibration(
        case_result,
        calibration;
        requested_physical_k_grid=requested_physical_k_grid,
    )
end
