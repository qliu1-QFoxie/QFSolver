function _resolve_tau_profile(
    bg::BackgroundSolution{T};
    tau_profile=nothing,
    tau_end::Union{Nothing,Real}=nothing,
    abs_tol::Union{Nothing,Real}=nothing,
    rel_tol::Union{Nothing,Real}=nothing,
) where {T<:AbstractFloat}
    if tau_profile === nothing
        return build_τ(
            bg;
            τ_end=tau_end,
            abs_tol=abs_tol,
            rel_tol=rel_tol,
        )
    end

    tau_profile.background === bg || throw(ArgumentError(
        "tau_profile must be built from the same background solution"))
    tau_profile
end

@inline function _bd_phase(k_value, tau_in, ::Type{T}) where {T<:AbstractFloat}
    exp(complex(zero(T), -(k_value * tau_in)))
end

function scalar_bd_initial_conditions(
    bg::BackgroundSolution{T},
    k::Real,
    N_in::Real;
    tau_profile=nothing,
    tau_end::Union{Nothing,Real}=nothing,
    abs_tol::Union{Nothing,Real}=nothing,
    rel_tol::Union{Nothing,Real}=nothing,
) where {T<:AbstractFloat}
    k_value = _checked_positive_finite(k, "k", T)
    N_value = _checked_finite(N_in, "N_in", T)
    tau_profile_value = _resolve_tau_profile(
        bg;
        tau_profile=tau_profile,
        tau_end=tau_end,
        abs_tol=abs_tol,
        rel_tol=rel_tol,
    )

    tau_in = τ(tau_profile_value, N_value)
    aH_in = aH(bg, N_value)
    phase = _bd_phase(k_value, tau_in, T)

    Q0 = exp(-N_value) * phase / sqrt(T(2) * k_value)
    dQ0 = (-one(T) - complex(zero(T), k_value / aH_in)) * Q0
    (Q0, dQ0)
end

function scalar_bd_initial_conditions(
    bg::BackgroundSolution{T},
    k::Real,
    window::ModeWindow{T};
    kwargs...,
) where {T<:AbstractFloat}
    scalar_bd_initial_conditions(bg, k, window.N_in; kwargs...)
end

function tensor_bd_initial_conditions(
    bg::BackgroundSolution{T},
    k::Real,
    N_in::Real;
    tau_profile=nothing,
    tau_end::Union{Nothing,Real}=nothing,
    abs_tol::Union{Nothing,Real}=nothing,
    rel_tol::Union{Nothing,Real}=nothing,
) where {T<:AbstractFloat}
    k_value = _checked_positive_finite(k, "k", T)
    N_value = _checked_finite(N_in, "N_in", T)
    tau_profile_value = _resolve_tau_profile(
        bg;
        tau_profile=tau_profile,
        tau_end=tau_end,
        abs_tol=abs_tol,
        rel_tol=rel_tol,
    )

    tau_in = τ(tau_profile_value, N_value)
    aH_in = aH(bg, N_value)
    phase = _bd_phase(k_value, tau_in, T)

    u0 = phase / sqrt(T(2) * k_value)
    du0 = -complex(zero(T), k_value / aH_in) * u0
    (u0, du0)
end

function tensor_bd_initial_conditions(
    bg::BackgroundSolution{T},
    k::Real,
    window::ModeWindow{T};
    kwargs...,
) where {T<:AbstractFloat}
    tensor_bd_initial_conditions(bg, k, window.N_in; kwargs...)
end
