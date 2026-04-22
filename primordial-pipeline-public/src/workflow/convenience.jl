@inline function _resolved_request_potential(::Type{T}, potential::QuadraticPotential) where {T<:AbstractFloat}
    QuadraticPotential(T(potential.m))
end

@inline function _resolved_request_potential(::Type{T}, potential::InflatonPotential) where {T<:AbstractFloat}
    potential
end

function build_exact_physical_spectrum_request(
    potential::InflatonPotential,
    ϕ₀,
    ϕ_N₀;
    background_config::BackgroundConfig{T}=BackgroundConfig{Float64}(),
    k_request=nothing,
    N_star::Real=55.0,
    k_pivot_phys_mpc::Real=0.05,
    half_width_log10::Real=0.5,
    ppd::Real=64.0,
    mode_config=nothing,
    endpoint_config=nothing,
    calibration_config=nothing,
) where {T<:AbstractFloat}
    endpoint_config_value = _resolved_exact_request_endpoint_config(T, endpoint_config)
    mode_config_value = _resolved_exact_request_mode_config(T, mode_config, endpoint_config_value)
    k_request_value = _resolved_exact_request_k_request(T, k_request, half_width_log10, ppd)
    _validate_physical_k_request(k_request_value)
    calibration_config_value = _resolved_exact_request_calibration_config(
        T,
        calibration_config,
        N_star,
        k_pivot_phys_mpc,
    )
    potential_value = _resolved_request_potential(T, potential)
    ϕ₀_value = _checked_finite(ϕ₀, "ϕ₀", T)
    ϕ_N₀_value = _checked_finite(ϕ_N₀, "ϕ_N₀", T)

    PrimordialPipelineRequest(
        potential=potential_value,
        ϕ₀=ϕ₀_value,
        ϕ_N₀=ϕ_N₀_value,
        background_config=background_config,
        k_request=k_request_value,
        mode_config=mode_config_value,
        endpoint_config=endpoint_config_value,
        calibration_config=calibration_config_value,
    )
end

function solve_exact_physical_spectrum(
    potential::InflatonPotential,
    ϕ₀,
    ϕ_N₀;
    kwargs...,
)
    request = build_exact_physical_spectrum_request(
        potential,
        ϕ₀,
        ϕ_N₀;
        kwargs...,
    )
    solve_primordial_pipeline(request)
end

function build_class_primordial_export(
    result::PrimordialPipelineResult;
    kwargs...,
)
    build_class_primordial_export(result.spectrum_case; kwargs...)
end
