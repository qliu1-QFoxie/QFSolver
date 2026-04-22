@inline ModeCoefficients(point::BackgroundPoint{T}) where {T<:AbstractFloat} =
    ModeCoefficients{T}(point.ε_H, point.aH, point.μ_Q)

function mode_coefficients(bg::BackgroundSolution{T}, N::Real) where {T<:AbstractFloat}
    ModeCoefficients(background_point(bg, N))
end

@inline scalar_mode_friction(coeffs::ModeCoefficients{T}) where {T<:AbstractFloat} =
    T(3) - coeffs.epsilon_H

@inline tensor_mode_friction(coeffs::ModeCoefficients{T}) where {T<:AbstractFloat} =
    one(T) - coeffs.epsilon_H

function scalar_mode_omega2(coeffs::ModeCoefficients{T}, k::Real) where {T<:AbstractFloat}
    k_value = _checked_positive_finite(k, "k", T)
    (k_value / coeffs.aH)^2 + coeffs.mu_Q
end

function tensor_mode_omega2(coeffs::ModeCoefficients{T}, k::Real) where {T<:AbstractFloat}
    k_value = _checked_positive_finite(k, "k", T)
    (k_value / coeffs.aH)^2 - (T(2) - coeffs.epsilon_H)
end

function scalar_canonical_omega2(coeffs::ModeCoefficients{T}, k::Real) where {T<:AbstractFloat}
    k_value = _checked_positive_finite(k, "k", T)
    k_value^2 + coeffs.aH^2 * (coeffs.mu_Q - (T(2) - coeffs.epsilon_H))
end

function tensor_canonical_omega2(coeffs::ModeCoefficients{T}, k::Real) where {T<:AbstractFloat}
    k_value = _checked_positive_finite(k, "k", T)
    k_value^2 - coeffs.aH^2 * (T(2) - coeffs.epsilon_H)
end

@inline function scalar_canonical_frequency_sq(
    bg::BackgroundSolution{T},
    k::Real,
    N::Real,
) where {T<:AbstractFloat}
    scalar_canonical_omega2(mode_coefficients(bg, N), k)
end

@inline function tensor_canonical_frequency_sq(
    bg::BackgroundSolution{T},
    k::Real,
    N::Real,
) where {T<:AbstractFloat}
    tensor_canonical_omega2(mode_coefficients(bg, N), k)
end

@inline function _adiabatic_infinity(::Type{T}) where {T<:AbstractFloat}
    T(Inf)
end

@inline function _adiabatic_failure_metrics(::Type{T}, omega_sq::T) where {T<:AbstractFloat}
    AdiabaticMetrics{T}(
        omega_sq=omega_sq,
        alpha1=_adiabatic_infinity(T),
        alpha2=_adiabatic_infinity(T),
        passed=false,
    )
end

@inline function _adiabatic_sample_omega(
    frequency_sq_fn,
    bg::BackgroundSolution{T},
    k_value::T,
    N::T,
) where {T<:AbstractFloat}
    omega_sq = frequency_sq_fn(bg, k_value, N)
    omega_sq > zero(T) || return nothing
    (omega_sq=omega_sq, omega=sqrt(omega_sq))
end

function _adiabatic_derivative_geometry(
    bg::BackgroundSolution{T},
    N::T,
    h::Real,
) where {T<:AbstractFloat}
    h_requested = _checked_positive_finite(h, "adiabatic_derivative_step", T)
    N_lo, N_hi = N_span(bg)
    backward_room = N - N_lo
    forward_room = N_hi - N
    central_room = min(backward_room, forward_room)

    if central_room > zero(T)
        return (:central, min(h_requested, central_room))
    end
    if forward_room > zero(T)
        return (:forward, min(h_requested, forward_room / T(2)))
    end
    if backward_room > zero(T)
        return (:backward, min(h_requested, backward_room / T(2)))
    end

    throw(ArgumentError(
        "adiabatic derivative evaluation requires a non-degenerate solved background interval"))
end

function _adiabatic_omega_derivatives(
    frequency_sq_fn,
    bg::BackgroundSolution{T},
    k_value::T,
    N_value::T;
    h::Real,
) where {T<:AbstractFloat}
    scheme, h_value = _adiabatic_derivative_geometry(bg, N_value, h)
    sample_0 = _adiabatic_sample_omega(frequency_sq_fn, bg, k_value, N_value)
    sample_0 === nothing && return nothing

    if scheme == :central
        sample_p = _adiabatic_sample_omega(frequency_sq_fn, bg, k_value, N_value + h_value)
        sample_m = _adiabatic_sample_omega(frequency_sq_fn, bg, k_value, N_value - h_value)
        (sample_p === nothing || sample_m === nothing) && return nothing
        dω_dN = (sample_p.omega - sample_m.omega) / (T(2) * h_value)
        d2ω_dN2 = (sample_p.omega - T(2) * sample_0.omega + sample_m.omega) / (h_value^2)
        return (omega_sq=sample_0.omega_sq, omega=sample_0.omega, dω_dN=dω_dN, d2ω_dN2=d2ω_dN2)
    end

    if scheme == :forward
        sample_p1 = _adiabatic_sample_omega(frequency_sq_fn, bg, k_value, N_value + h_value)
        sample_p2 = _adiabatic_sample_omega(frequency_sq_fn, bg, k_value, N_value + T(2) * h_value)
        (sample_p1 === nothing || sample_p2 === nothing) && return nothing
        dω_dN = (-T(3) * sample_0.omega + T(4) * sample_p1.omega - sample_p2.omega) / (T(2) * h_value)
        d2ω_dN2 = (sample_0.omega - T(2) * sample_p1.omega + sample_p2.omega) / (h_value^2)
        return (omega_sq=sample_0.omega_sq, omega=sample_0.omega, dω_dN=dω_dN, d2ω_dN2=d2ω_dN2)
    end

    sample_m1 = _adiabatic_sample_omega(frequency_sq_fn, bg, k_value, N_value - h_value)
    sample_m2 = _adiabatic_sample_omega(frequency_sq_fn, bg, k_value, N_value - T(2) * h_value)
    (sample_m1 === nothing || sample_m2 === nothing) && return nothing
    dω_dN = (T(3) * sample_0.omega - T(4) * sample_m1.omega + sample_m2.omega) / (T(2) * h_value)
    d2ω_dN2 = (sample_0.omega - T(2) * sample_m1.omega + sample_m2.omega) / (h_value^2)
    (omega_sq=sample_0.omega_sq, omega=sample_0.omega, dω_dN=dω_dN, d2ω_dN2=d2ω_dN2)
end

function _canonical_adiabatic_metrics(
    frequency_sq_fn,
    bg::BackgroundSolution{T},
    k::Real,
    N::Real;
    h::Real,
    alpha1_tol::Real,
    alpha2_tol::Real,
) where {T<:AbstractFloat}
    k_value = _checked_positive_finite(k, "k", T)
    N_value = _checked_finite(N, "N", T)
    alpha1_tol_value = _checked_positive_finite(alpha1_tol, "adiabatic_alpha1_tol", T)
    alpha2_tol_value = _checked_positive_finite(alpha2_tol, "adiabatic_alpha2_tol", T)

    derivative_data = _adiabatic_omega_derivatives(
        frequency_sq_fn,
        bg,
        k_value,
        N_value;
        h=h,
    )
    derivative_data === nothing && return _adiabatic_failure_metrics(
        T,
        frequency_sq_fn(bg, k_value, N_value),
    )

    coeffs = mode_coefficients(bg, N_value)
    omega_prime = coeffs.aH * derivative_data.dω_dN
    omega_double_prime = coeffs.aH^2 * (
        derivative_data.d2ω_dN2 + (one(T) - coeffs.epsilon_H) * derivative_data.dω_dN
    )

    alpha1 = abs(omega_prime / derivative_data.omega_sq)
    alpha2 = abs(omega_double_prime / (derivative_data.omega^3))
    passed = (
        derivative_data.omega_sq > zero(T) &&
        alpha1 <= alpha1_tol_value &&
        alpha2 <= alpha2_tol_value
    )

    AdiabaticMetrics{T}(
        omega_sq=derivative_data.omega_sq,
        alpha1=alpha1,
        alpha2=alpha2,
        passed=passed,
    )
end

function scalar_adiabatic_metrics(
    bg::BackgroundSolution{T},
    k::Real,
    N::Real;
    h::Real=ModeConfig{T}().adiabatic_derivative_step,
    alpha1_tol::Real=ModeConfig{T}().adiabatic_alpha1_tol,
    alpha2_tol::Real=ModeConfig{T}().adiabatic_alpha2_tol,
) where {T<:AbstractFloat}
    _canonical_adiabatic_metrics(
        scalar_canonical_frequency_sq,
        bg,
        k,
        N;
        h=h,
        alpha1_tol=alpha1_tol,
        alpha2_tol=alpha2_tol,
    )
end

function tensor_adiabatic_metrics(
    bg::BackgroundSolution{T},
    k::Real,
    N::Real;
    h::Real=ModeConfig{T}().adiabatic_derivative_step,
    alpha1_tol::Real=ModeConfig{T}().adiabatic_alpha1_tol,
    alpha2_tol::Real=ModeConfig{T}().adiabatic_alpha2_tol,
) where {T<:AbstractFloat}
    _canonical_adiabatic_metrics(
        tensor_canonical_frequency_sq,
        bg,
        k,
        N;
        h=h,
        alpha1_tol=alpha1_tol,
        alpha2_tol=alpha2_tol,
    )
end
