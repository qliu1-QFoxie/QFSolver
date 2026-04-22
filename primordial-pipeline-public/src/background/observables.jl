@inline function _checked_N(bg::BackgroundSolution{T}, N::Real) where {T<:AbstractFloat}
    N_value = T(N)
    isfinite(N_value) || throw(ArgumentError("Background evaluation point N must be finite"))

    N_lo, N_hi = N_span(bg)
    tol = _N_tolerance(T, N_lo, N_hi)
    if N_value < N_lo - tol || N_value > N_hi + tol
        throw(DomainError(N,
            "Requested N=$(N) lies outside the solved background interval [$N_lo, $N_hi]"))
    end

    clamp(N_value, N_lo, N_hi)
end

@inline function _state(bg::BackgroundSolution, N)
    T = typeof(bg.config.M_pl)
    N_final = T(bg.ode_solution.t[end])
    tol = _N_tolerance(T, T(bg.config.N_start), N_final)
    abs(N - N_final) <= tol ? bg.ode_solution(N_final; continuity=:right) : bg.ode_solution(N)
end

@inline function _η_H(ϕ_N, ϕ_NN, M_pl)
    T = typeof(M_pl)
    ϕ_N_tol = sqrt(eps(T)) * M_pl
    abs(ϕ_N) > ϕ_N_tol || return T(NaN)
    T(2) * ϕ_NN / ϕ_N
end

@inline function _μ_Q(V_ϕϕ_value, H²_value, ε_H_value, ϕ_N_value, ϕ_NN_value, M_pl²)
    T = typeof(M_pl²)
    V_ϕϕ_value / H²_value -
        T(2) * ε_H_value * (T(3) - ε_H_value) -
        T(2) * ϕ_N_value * ϕ_NN_value / M_pl²
end

function background_point(bg::BackgroundSolution{T}, N::Real) where {T<:AbstractFloat}
    N_value = _checked_N(bg, N)
    state = _state(bg, N_value)
    M_pl = bg.config.M_pl
    M_pl² = M_pl^2

    ϕ_value = state[1]
    ϕ_N_value = state[2]

    V_value = V(bg.potential, ϕ_value)
    V_ϕ_value = V_ϕ(bg.potential, ϕ_value)
    V_ϕϕ_value = V_ϕϕ(bg.potential, ϕ_value)
    ε_H_value = _ε_H(ϕ_N_value, M_pl²)
    H²_value = _H²(bg.potential, ϕ_value, ϕ_N_value, M_pl²)
    H_value = sqrt(H²_value)
    ϕ_NN_value = -(T(3) - ε_H_value) * ϕ_N_value - V_ϕ_value / H²_value
    η_H_value = _η_H(ϕ_N_value, ϕ_NN_value, M_pl)
    a_value = exp(N_value)
    aH_value = a_value * H_value
    μ_Q_value = _μ_Q(V_ϕϕ_value, H²_value, ε_H_value, ϕ_N_value, ϕ_NN_value, M_pl²)

    BackgroundPoint(
        ϕ_value,
        ϕ_N_value,
        ϕ_NN_value,
        V_value,
        V_ϕ_value,
        V_ϕϕ_value,
        H²_value,
        H_value,
        ε_H_value,
        η_H_value,
        a_value,
        aH_value,
        μ_Q_value,
    )
end

ϕ(bg::BackgroundSolution, N::Real) = background_point(bg, N).ϕ
ϕ_N(bg::BackgroundSolution, N::Real) = background_point(bg, N).ϕ_N
ϕ_NN(bg::BackgroundSolution, N::Real) = background_point(bg, N).ϕ_NN

V(bg::BackgroundSolution, N::Real) = background_point(bg, N).V
V_ϕ(bg::BackgroundSolution, N::Real) = background_point(bg, N).V_ϕ
V_ϕϕ(bg::BackgroundSolution, N::Real) = background_point(bg, N).V_ϕϕ

H²(bg::BackgroundSolution, N::Real) = background_point(bg, N).H²
H(bg::BackgroundSolution, N::Real) = background_point(bg, N).H
ε_H(bg::BackgroundSolution, N::Real) = background_point(bg, N).ε_H
η_H(bg::BackgroundSolution, N::Real) = background_point(bg, N).η_H
a(bg::BackgroundSolution, N::Real) = background_point(bg, N).a
aH(bg::BackgroundSolution, N::Real) = background_point(bg, N).aH

μ_Q(bg::BackgroundSolution, N::Real) = background_point(bg, N).μ_Q
