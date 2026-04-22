@inline function _checked_positive_finite_observable(
    value::Real,
    name::AbstractString,
)
    isfinite(value) || throw(ArgumentError("$(name) must be finite"))
    value > zero(value) || throw(ArgumentError("$(name) must be strictly positive"))
    value
end

function exact_scalar_endpoint_observable(
    solution::ScalarModeSolution{T},
    N::Real,
) where {T<:AbstractFloat}
    N_value = T(N)
    k_value = solution.k
    bg = solution.background

    Q_k = scalar_mode_value(solution, N_value)
    point = background_point(bg, N_value)
    ϕ_N_value = point.ϕ_N
    ϕ_N_value != zero(T) || throw(ArgumentError(
        "scalar exact endpoint observable is undefined when ϕ_N == 0"))

    R_k = -Q_k / ϕ_N_value
    prefactor = k_value^3 / (T(2) * T(pi)^2)
    _checked_positive_finite_observable(prefactor * abs2(R_k), "scalar endpoint observable")
end

function exact_tensor_endpoint_observable(
    solution::TensorModeSolution{T},
    N::Real,
) where {T<:AbstractFloat}
    N_value = T(N)
    k_value = solution.k
    bg = solution.background

    u_k = tensor_mode_value(solution, N_value)
    a_value = a(bg, N_value)
    M_pl = bg.config.M_pl
    M_pl^2 > zero(T) || throw(ArgumentError(
        "tensor exact endpoint observable requires M_pl^2 > 0"))

    h_proxy = u_k / a_value
    prefactor = T(4) * k_value^3 / (T(pi)^2 * M_pl^2)
    _checked_positive_finite_observable(prefactor * abs2(h_proxy), "tensor endpoint observable")
end

function exact_scalar_endpoint_log_slope(
    solution::ScalarModeSolution{T},
    N::Real,
) where {T<:AbstractFloat}
    N_value = T(N)
    bg = solution.background

    Q_k = scalar_mode_value(solution, N_value)
    Q_N = scalar_mode_derivative(solution, N_value)
    point = background_point(bg, N_value)
    ϕ_N_value = point.ϕ_N
    ϕ_NN_value = point.ϕ_NN
    ϕ_N_value != zero(T) || throw(ArgumentError(
        "scalar endpoint log slope is undefined when ϕ_N == 0"))

    R_k = -Q_k / ϕ_N_value
    abs2_R = abs2(R_k)
    abs2_R > zero(T) || throw(ArgumentError(
        "scalar endpoint log slope is undefined when |R_k|^2 == 0"))

    R_N = -(Q_N * ϕ_N_value - Q_k * ϕ_NN_value) / (ϕ_N_value^2)
    slope = abs(T(2) * real(R_N * conj(R_k)) / abs2_R)
    isfinite(slope) || throw(ArgumentError("scalar endpoint log slope must be finite"))
    slope
end

function exact_tensor_endpoint_log_slope(
    solution::TensorModeSolution{T},
    N::Real,
) where {T<:AbstractFloat}
    N_value = T(N)
    bg = solution.background

    u_k = tensor_mode_value(solution, N_value)
    u_N = tensor_mode_derivative(solution, N_value)
    a_value = a(bg, N_value)
    h_proxy = u_k / a_value
    abs2_h = abs2(h_proxy)
    abs2_h > zero(T) || throw(ArgumentError(
        "tensor endpoint log slope is undefined when |u_k/a|^2 == 0"))

    h_proxy_N = (u_N - u_k) / a_value
    slope = abs(T(2) * real(h_proxy_N * conj(h_proxy)) / abs2_h)
    isfinite(slope) || throw(ArgumentError("tensor endpoint log slope must be finite"))
    slope
end
