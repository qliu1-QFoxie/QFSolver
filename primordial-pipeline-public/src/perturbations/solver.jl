struct _CheckedDenseSolution{T<:AbstractFloat,S}
    ode_solution::S
    N_lo::T
    N_hi::T
end

struct _ScalarCanonicalODESolution{T<:AbstractFloat,S,D}
    raw_solution::S
    t::Vector{T}
    u::Vector{Vector{Complex{T}}}
    destats::D
end

struct ScalarModeSolution{T<:AbstractFloat,B,S}
    background::B
    k::T
    window::ModeWindow{T}
    solution::S
    value_stop::Complex{T}
    derivative_stop::Complex{T}
end

struct TensorModeSolution{T<:AbstractFloat,B,S}
    background::B
    k::T
    window::ModeWindow{T}
    solution::S
    value_stop::Complex{T}
    derivative_stop::Complex{T}
end

@inline function _checked_mode_N(::Type{T}, N_lo, N_hi, N::Real) where {T<:AbstractFloat}
    N_value = T(N)
    isfinite(N_value) || throw(ArgumentError("mode evaluation point N must be finite"))

    tol = max(T(100) * eps(T) * max(abs(N_lo), abs(N_hi), one(T)), sqrt(eps(T)))
    if N_value < N_lo - tol || N_value > N_hi + tol
        throw(DomainError(N,
            "Requested N=$(N) lies outside the mode interval [$N_lo, $N_hi]"))
    end

    clamp(N_value, N_lo, N_hi)
end

@inline function (solution::_CheckedDenseSolution{T})(N::Real) where {T<:AbstractFloat}
    solution.ode_solution(_checked_mode_N(T, solution.N_lo, solution.N_hi, N))
end

@inline function (solution::_ScalarCanonicalODESolution{T})(N::Real) where {T<:AbstractFloat}
    _scalar_mode_complex_state(solution.raw_solution(T(N)))
end

@inline function _mode_solver_kwargs(dt_min::T, dt_max::T) where {T<:AbstractFloat}
    has_dt_min = dt_min > zero(T)
    has_dt_max = dt_max > zero(T)

    has_dt_min && has_dt_max && return (; dtmin=dt_min, dtmax=dt_max)
    has_dt_min && return (; dtmin=dt_min)
    has_dt_max && return (; dtmax=dt_max)
    (;)
end

@inline function _validate_mode_solver_config(config::ModeConfig{T}) where {T<:AbstractFloat}
    config.solver === :vern9 || throw(ArgumentError(
        "Only solver=:vern9 is supported for exact single-mode solves"))
    config.max_iters > 0 || throw(ArgumentError("max_iters must be positive"))
    nothing
end

@inline function _scalar_mode_state_scale_floor(
    ::Type{T},
    v0,
    v_N0,
) where {T<:AbstractFloat}
    sqrt(eps(T)) * max(abs(v0), abs(v_N0), one(T))
end

@inline function _scalar_mode_state_scales(
    v0,
    v_N0,
    config::ModeConfig{T},
) where {T<:AbstractFloat}
    floor = _scalar_mode_state_scale_floor(T, v0, v_N0)
    (
        v=max(abs(v0), floor),
        v_N=max(abs(v_N0), floor),
    )
end

@inline function _scalar_mode_component_abstol(
    scales,
    config::ModeConfig{T},
) where {T<:AbstractFloat}
    T[
        config.abs_tol * scales.v,
        config.abs_tol * scales.v,
        config.abs_tol * scales.v_N,
        config.abs_tol * scales.v_N,
    ]
end

@inline _scalar_mode_scale_max(scales) = max(scales.v, scales.v_N)

@inline function _scalar_mode_scaled_maxnorm(u::Number, scales)
    abs(u) / _scalar_mode_scale_max(scales)
end

@inline function _scalar_mode_scaled_maxnorm(u, scales)
    max(
        abs(u[1]) / scales.v,
        abs(u[2]) / scales.v,
        abs(u[3]) / scales.v_N,
        abs(u[4]) / scales.v_N,
    )
end

@inline _scalar_mode_internalnorm(scales) = (u, N) -> _scalar_mode_scaled_maxnorm(u, scales)

@inline function _scalar_mode_real_state(v::Complex{T}, v_N::Complex{T}) where {T<:AbstractFloat}
    T[real(v), imag(v), real(v_N), imag(v_N)]
end

@inline function _scalar_mode_complex_state(state::AbstractVector{T}) where {T<:AbstractFloat}
    Complex{T}[
        Complex{T}(state[1], state[2]),
        Complex{T}(state[3], state[4]),
    ]
end

function _wrap_scalar_mode_ode_solution(raw_solution, ::Type{T}) where {T<:AbstractFloat}
    _ScalarCanonicalODESolution{T,typeof(raw_solution),typeof(raw_solution.destats)}(
        raw_solution,
        T.(raw_solution.t),
        [_scalar_mode_complex_state(state) for state in raw_solution.u],
        raw_solution.destats,
    )
end

@inline function _validate_mode_window(window::ModeWindow{T}) where {T<:AbstractFloat}
    window.N_in < window.N_hc || throw(ArgumentError(
        "ModeWindow must satisfy N_in < N_hc"))
    window.N_hc < window.N_stop || throw(ArgumentError(
        "ModeWindow must satisfy N_hc < N_stop"))
    nothing
end

function _resolve_mode_window(
    bg::BackgroundSolution{T},
    k::Real;
    config::ModeConfig{T},
    window::Union{Nothing,ModeWindow{T}}=nothing,
) where {T<:AbstractFloat}
    if window === nothing
        return mode_window(bg, k; config=config)
    end

    _validate_mode_window(window)
    N_lo, N_hi = N_span(bg)
    window.N_in >= N_lo || throw(ArgumentError("ModeWindow.N_in lies before the solved background interval"))
    window.N_stop <= N_hi || throw(ArgumentError("ModeWindow.N_stop lies after the solved background interval"))
    if config.init_requires_adiabatic && !_shared_adiabaticity_pass(bg, T(k), window.N_in, config)
        throw(ArgumentError(
            "supplied ModeWindow.N_in does not satisfy the exact adiabaticity gate required for plane-wave BD initialization"))
    end
    window
end

@inline function _reference_ratio_start(
    bg::BackgroundSolution{T},
    k_value::T,
    config::ModeConfig{T},
    fallback_N::T,
) where {T<:AbstractFloat}
    try
        N_at_k_over_aH(bg, k_value, config.subhorizon_ratio)
    catch err
        err isa ArgumentError || rethrow()
        fallback_N
    end
end

@inline function _adiabatic_metrics_for_kind(
    kind::PerturbationKind,
    bg::BackgroundSolution{T},
    k_value::T,
    N_value::T,
    config::ModeConfig{T},
) where {T<:AbstractFloat}
    if kind == scalar
        return scalar_adiabatic_metrics(
            bg,
            k_value,
            N_value;
            h=config.adiabatic_derivative_step,
            alpha1_tol=config.adiabatic_alpha1_tol,
            alpha2_tol=config.adiabatic_alpha2_tol,
        )
    elseif kind == tensor
        return tensor_adiabatic_metrics(
            bg,
            k_value,
            N_value;
            h=config.adiabatic_derivative_step,
            alpha1_tol=config.adiabatic_alpha1_tol,
            alpha2_tol=config.adiabatic_alpha2_tol,
        )
    end

    throw(ArgumentError("unsupported perturbation kind for initialization diagnostics"))
end

function _mode_initialization_diagnostics(
    kind::PerturbationKind,
    bg::BackgroundSolution{T},
    k_value::T,
    N_in_used::T,
    config::ModeConfig{T};
    failure_reason::Union{Nothing,String}=nothing,
    N_in_ratio::Union{Nothing,T}=nothing,
) where {T<:AbstractFloat}
    N_ratio = N_in_ratio === nothing ? _reference_ratio_start(bg, k_value, config, N_in_used) : N_in_ratio
    metrics = _adiabatic_metrics_for_kind(kind, bg, k_value, N_in_used, config)

    ModeInitializationDiagnostics{T}(
        N_in_ratio=N_ratio,
        N_in_used=N_in_used,
        q_in=k_value / aH(bg, N_in_used),
        init_omega_sq=metrics.omega_sq,
        init_alpha1=metrics.alpha1,
        init_alpha2=metrics.alpha2,
        init_adiabatic_passed=metrics.passed,
        init_kind="plane_wave_bd",
        init_failure_reason=failure_reason,
    )
end

function _solve_mode_problem(
    kind::PerturbationKind,
    k::T,
    window::ModeWindow{T},
    rhs!,
    initial_state,
    tspan::Tuple{T,T},
    config::ModeConfig{T};
    abstol=config.abs_tol,
    internalnorm=nothing,
) where {T<:AbstractFloat}
    problem = ODEProblem(rhs!, initial_state, tspan)
    solution = if internalnorm === nothing
        solve(
            problem,
            Vern9();
            abstol=abstol,
            reltol=config.rel_tol,
            dense=true,
            save_everystep=true,
            maxiters=config.max_iters,
            _mode_solver_kwargs(config.dt_min, config.dt_max)...,
        )
    else
        solve(
            problem,
            Vern9();
            abstol=abstol,
            reltol=config.rel_tol,
            dense=true,
            save_everystep=true,
            maxiters=config.max_iters,
            internalnorm=internalnorm,
            _mode_solver_kwargs(config.dt_min, config.dt_max)...,
        )
    end

    solution.retcode == ReturnCode.Success || throw(ModeSolveFailure(
        kind,
        k,
        window,
        solution.retcode,
    ))
    solution
end

function _scalar_mode_canonical_initial_conditions(
    bg::BackgroundSolution{T},
    k::Real,
    window::ModeWindow{T};
    tau_profile=nothing,
    tau_end::Union{Nothing,Real}=nothing,
    abs_tol::Union{Nothing,Real}=nothing,
    rel_tol::Union{Nothing,Real}=nothing,
) where {T<:AbstractFloat}
    Q0, Q_N0 = scalar_bd_initial_conditions(
        bg,
        k,
        window;
        tau_profile=tau_profile,
        tau_end=tau_end,
        abs_tol=abs_tol,
        rel_tol=rel_tol,
    )
    a_in = exp(window.N_in)
    v0 = a_in * Q0
    v_N0 = v0 + a_in * Q_N0
    (v0, v_N0)
end

function _scalar_mode_canonical_rhs(
    bg::BackgroundSolution{T},
    N::Real,
    v,
    v_N,
    k::Real,
) where {T<:AbstractFloat}
    coeffs = mode_coefficients(bg, N)
    omega2 = scalar_mode_omega2(coeffs, k) - (T(2) - coeffs.epsilon_H)
    (v_N, -tensor_mode_friction(coeffs) * v_N - omega2 * v)
end

@inline function _checked_exact_precision_type(
    ::Type{T},
    precision::Type,
) where {T<:AbstractFloat}
    precision === T || throw(ArgumentError(
        "solve_exact requires precision=$(T) to match the solved background " *
        "and mode configuration precision; requested precision=$(precision)"))
    precision
end

function solve_scalar_mode(
    bg::BackgroundSolution{T},
    k::Real;
    config::ModeConfig{T}=ModeConfig{T}(),
    window::Union{Nothing,ModeWindow{T}}=nothing,
    tau_profile=nothing,
    tau_end::Union{Nothing,Real}=nothing,
    abs_tol::Union{Nothing,Real}=nothing,
    rel_tol::Union{Nothing,Real}=nothing,
) where {T<:AbstractFloat}
    _validate_mode_solver_config(config)
    k_value = _checked_positive_finite(k, "k", T)
    window_value = _resolve_mode_window(bg, k_value; config=config, window=window)
    v0, v_N0 = _scalar_mode_canonical_initial_conditions(
        bg,
        k_value,
        window_value;
        tau_profile=tau_profile,
        tau_end=tau_end,
        abs_tol=abs_tol,
        rel_tol=rel_tol,
    )
    initial_state = _scalar_mode_real_state(v0, v_N0)
    scales = _scalar_mode_state_scales(v0, v_N0, config)

    rhs! = function (du, u, p, N)
        state = _scalar_mode_complex_state(u)
        derivative = _scalar_mode_canonical_rhs(bg, N, state[1], state[2], k_value)
        du[1] = real(derivative[1])
        du[2] = imag(derivative[1])
        du[3] = real(derivative[2])
        du[4] = imag(derivative[2])
        nothing
    end

    raw_ode_solution = _solve_mode_problem(
        scalar,
        k_value,
        window_value,
        rhs!,
        initial_state,
        (window_value.N_in, window_value.N_stop),
        config;
        abstol=_scalar_mode_component_abstol(scales, config),
        internalnorm=_scalar_mode_internalnorm(scales),
    )
    ode_solution = _wrap_scalar_mode_ode_solution(raw_ode_solution, T)
    checked_solution = _CheckedDenseSolution{T,typeof(ode_solution)}(
        ode_solution,
        window_value.N_in,
        window_value.N_stop,
    )
    stop_state = checked_solution(window_value.N_stop)
    a_stop = exp(window_value.N_stop)
    Q_stop = stop_state[1] / a_stop
    Q_N_stop = (stop_state[2] - stop_state[1]) / a_stop

    ScalarModeSolution(
        bg,
        k_value,
        window_value,
        checked_solution,
        Complex{T}(Q_stop),
        Complex{T}(Q_N_stop),
    )
end

function solve_tensor_mode(
    bg::BackgroundSolution{T},
    k::Real;
    config::ModeConfig{T}=ModeConfig{T}(),
    window::Union{Nothing,ModeWindow{T}}=nothing,
    tau_profile=nothing,
    tau_end::Union{Nothing,Real}=nothing,
    abs_tol::Union{Nothing,Real}=nothing,
    rel_tol::Union{Nothing,Real}=nothing,
) where {T<:AbstractFloat}
    _validate_mode_solver_config(config)
    k_value = _checked_positive_finite(k, "k", T)
    window_value = _resolve_mode_window(bg, k_value; config=config, window=window)
    u0, u_N0 = tensor_bd_initial_conditions(
        bg,
        k_value,
        window_value;
        tau_profile=tau_profile,
        tau_end=tau_end,
        abs_tol=abs_tol,
        rel_tol=rel_tol,
    )

    rhs! = function (du, u, p, N)
        state = tensor_mode_rhs(bg, N, u[1], u[2], k_value)
        du[1] = state[1]
        du[2] = state[2]
        nothing
    end

    ode_solution = _solve_mode_problem(
        tensor,
        k_value,
        window_value,
        rhs!,
        Complex{T}[u0, u_N0],
        (window_value.N_in, window_value.N_stop),
        config,
    )
    checked_solution = _CheckedDenseSolution{T,typeof(ode_solution)}(
        ode_solution,
        window_value.N_in,
        window_value.N_stop,
    )
    stop_state = checked_solution(window_value.N_stop)

    TensorModeSolution(
        bg,
        k_value,
        window_value,
        checked_solution,
        Complex{T}(stop_state[1]),
        Complex{T}(stop_state[2]),
    )
end

function solve_exact(
    kind::PerturbationKind,
    bg::BackgroundSolution{T},
    k::Real;
    precision::Type{<:AbstractFloat}=T,
    config::ModeConfig{T}=ModeConfig{T}(),
    window::Union{Nothing,ModeWindow{T}}=nothing,
    tau_profile=nothing,
    tau_end::Union{Nothing,Real}=nothing,
    abs_tol::Union{Nothing,Real}=nothing,
    rel_tol::Union{Nothing,Real}=nothing,
) where {T<:AbstractFloat}
    _checked_exact_precision_type(T, precision)
    k_value = _checked_positive_finite(k, "k", T)

    solution = if kind == scalar
        solve_scalar_mode(
            bg,
            k_value;
            config=config,
            window=window,
            tau_profile=tau_profile,
            tau_end=tau_end,
            abs_tol=abs_tol,
            rel_tol=rel_tol,
        )
    elseif kind == tensor
        solve_tensor_mode(
            bg,
            k_value;
            config=config,
            window=window,
            tau_profile=tau_profile,
            tau_end=tau_end,
            abs_tol=abs_tol,
            rel_tol=rel_tol,
        )
    else
        throw(ArgumentError("solve_exact requires kind to be scalar or tensor"))
    end

    initialization = _mode_initialization_diagnostics(
        kind,
        bg,
        k_value,
        solution.window.N_in,
        config,
    )

    ExactSolveResult(config, solution, initialization)
end

@inline scalar_mode_canonical_state(solution::ScalarModeSolution, N::Real) = solution.solution(N)
@inline scalar_mode_canonical_value(solution::ScalarModeSolution, N::Real) =
    scalar_mode_canonical_state(solution, N)[1]
@inline scalar_mode_canonical_derivative(solution::ScalarModeSolution, N::Real) =
    scalar_mode_canonical_state(solution, N)[2]

function scalar_mode_state(solution::ScalarModeSolution{T}, N::Real) where {T<:AbstractFloat}
    N_value = _checked_mode_N(T, solution.window.N_in, solution.window.N_stop, N)
    state = scalar_mode_canonical_state(solution, N_value)
    a_value = exp(N_value)
    Complex{T}[state[1] / a_value, (state[2] - state[1]) / a_value]
end

@inline scalar_mode_value(solution::ScalarModeSolution, N::Real) = scalar_mode_state(solution, N)[1]
@inline scalar_mode_derivative(solution::ScalarModeSolution, N::Real) = scalar_mode_state(solution, N)[2]

@inline tensor_mode_state(solution::TensorModeSolution, N::Real) = solution.solution(N)
@inline tensor_mode_value(solution::TensorModeSolution, N::Real) = tensor_mode_state(solution, N)[1]
@inline tensor_mode_derivative(solution::TensorModeSolution, N::Real) = tensor_mode_state(solution, N)[2]
