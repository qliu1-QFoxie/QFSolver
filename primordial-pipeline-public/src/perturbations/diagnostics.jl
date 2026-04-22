@inline _canonical_wronskian_target(::Type{T}) where {T<:AbstractFloat} = complex(zero(T), one(T))

@inline _raw_mode_wronskian(value, derivative) = value * conj(derivative) - conj(value) * derivative

@inline scalar_mode_wronskian(Q, Q_N) = _raw_mode_wronskian(Q, Q_N)

function scalar_mode_wronskian(solution::ScalarModeSolution, N::Real)
    state = scalar_mode_state(solution, N)
    scalar_mode_wronskian(state[1], state[2])
end

function scalar_mode_canonical_wronskian(
    coeffs::ModeCoefficients{T},
    N::Real,
    Q,
    Q_N,
) where {T<:AbstractFloat}
    N_value = _checked_finite(N, "N", T)
    exp(T(2) * N_value) * coeffs.aH * scalar_mode_wronskian(Q, Q_N)
end

function scalar_mode_canonical_wronskian(
    bg::BackgroundSolution{T},
    N::Real,
    Q,
    Q_N,
) where {T<:AbstractFloat}
    scalar_mode_canonical_wronskian(mode_coefficients(bg, N), N, Q, Q_N)
end

function scalar_mode_canonical_wronskian(solution::ScalarModeSolution{T}, N::Real) where {T<:AbstractFloat}
    coeffs = mode_coefficients(solution.background, N)
    state = scalar_mode_canonical_state(solution, N)
    coeffs.aH * scalar_mode_wronskian(state[1], state[2])
end

function scalar_mode_wronskian_residual(
    coeffs::ModeCoefficients{T},
    N::Real,
    Q,
    Q_N,
) where {T<:AbstractFloat}
    scalar_mode_canonical_wronskian(coeffs, N, Q, Q_N) - _canonical_wronskian_target(T)
end

function scalar_mode_wronskian_residual(
    bg::BackgroundSolution{T},
    N::Real,
    Q,
    Q_N,
) where {T<:AbstractFloat}
    scalar_mode_canonical_wronskian(bg, N, Q, Q_N) - _canonical_wronskian_target(T)
end

function scalar_mode_wronskian_residual(solution::ScalarModeSolution{T}, N::Real) where {T<:AbstractFloat}
    scalar_mode_canonical_wronskian(solution, N) - _canonical_wronskian_target(T)
end

@inline tensor_mode_wronskian(u, u_N) = _raw_mode_wronskian(u, u_N)

function tensor_mode_wronskian(solution::TensorModeSolution, N::Real)
    state = tensor_mode_state(solution, N)
    tensor_mode_wronskian(state[1], state[2])
end

function tensor_mode_canonical_wronskian(
    coeffs::ModeCoefficients{T},
    u,
    u_N,
) where {T<:AbstractFloat}
    coeffs.aH * tensor_mode_wronskian(u, u_N)
end

function tensor_mode_canonical_wronskian(
    bg::BackgroundSolution{T},
    N::Real,
    u,
    u_N,
) where {T<:AbstractFloat}
    tensor_mode_canonical_wronskian(mode_coefficients(bg, N), u, u_N)
end

function tensor_mode_canonical_wronskian(solution::TensorModeSolution{T}, N::Real) where {T<:AbstractFloat}
    state = tensor_mode_state(solution, N)
    tensor_mode_canonical_wronskian(solution.background, N, state[1], state[2])
end

function tensor_mode_wronskian_residual(
    coeffs::ModeCoefficients{T},
    u,
    u_N,
) where {T<:AbstractFloat}
    tensor_mode_canonical_wronskian(coeffs, u, u_N) - _canonical_wronskian_target(T)
end

function tensor_mode_wronskian_residual(
    bg::BackgroundSolution{T},
    N::Real,
    u,
    u_N,
) where {T<:AbstractFloat}
    tensor_mode_canonical_wronskian(bg, N, u, u_N) - _canonical_wronskian_target(T)
end

function tensor_mode_wronskian_residual(solution::TensorModeSolution{T}, N::Real) where {T<:AbstractFloat}
    tensor_mode_canonical_wronskian(solution, N) - _canonical_wronskian_target(T)
end
