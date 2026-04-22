@inline function _scalar_mode_acceleration(
    Q,
    Q_N,
    coeffs::ModeCoefficients{T},
    k::Real,
) where {T<:AbstractFloat}
    -scalar_mode_friction(coeffs) * Q_N - scalar_mode_omega2(coeffs, k) * Q
end

@inline function _tensor_mode_acceleration(
    u,
    u_N,
    coeffs::ModeCoefficients{T},
    k::Real,
) where {T<:AbstractFloat}
    -tensor_mode_friction(coeffs) * u_N - tensor_mode_omega2(coeffs, k) * u
end

function scalar_mode_rhs(
    Q,
    Q_N,
    coeffs::ModeCoefficients{T},
    k::Real,
) where {T<:AbstractFloat}
    (Q_N, _scalar_mode_acceleration(Q, Q_N, coeffs, k))
end

function scalar_mode_rhs(
    bg::BackgroundSolution{T},
    N::Real,
    Q,
    Q_N,
    k::Real,
) where {T<:AbstractFloat}
    scalar_mode_rhs(Q, Q_N, mode_coefficients(bg, N), k)
end

function tensor_mode_rhs(
    u,
    u_N,
    coeffs::ModeCoefficients{T},
    k::Real,
) where {T<:AbstractFloat}
    (u_N, _tensor_mode_acceleration(u, u_N, coeffs, k))
end

function tensor_mode_rhs(
    bg::BackgroundSolution{T},
    N::Real,
    u,
    u_N,
    k::Real,
) where {T<:AbstractFloat}
    tensor_mode_rhs(u, u_N, mode_coefficients(bg, N), k)
end
