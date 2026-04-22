abstract type AbstractPhysicalKRequest end

Base.@kwdef struct ExplicitPhysicalKIntervalRequest{T<:AbstractFloat} <: AbstractPhysicalKRequest
    k_min_phys_mpc::T
    k_max_phys_mpc::T
    ppd::Float64
end

Base.@kwdef struct PivotCenteredPhysicalKRequest{T<:AbstractFloat} <: AbstractPhysicalKRequest
    half_width_log10::T
    ppd::Float64
end

Base.@kwdef struct PrimordialPipelineRequest{P,B,K,M,E,C}
    potential::P
    ϕ₀
    ϕ_N₀
    background_config::B
    k_request::K
    mode_config::M
    endpoint_config::E
    calibration_config::C
end

Base.@kwdef struct PrimordialPipelineResult{R,B,T,S}
    request::R
    background::B
    tau_profile::T
    spectrum_case::S
end
