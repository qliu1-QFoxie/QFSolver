function solve_background(
    potential::InflatonPotential,
    ϕ₀::Real,
    ϕ_N₀::Real;
    config::BackgroundConfig,
)
    _validate_config(config)

    T = typeof(config.M_pl)
    ϕ₀_T = T(ϕ₀)
    ϕ_N₀_T = T(ϕ_N₀)
    _validate_initial_state(potential, ϕ₀_T, ϕ_N₀_T, config)
    scales = _state_scales(ϕ₀_T, ϕ_N₀_T, config)
    abstol_vec = _component_abstol(scales, config)
    internalnorm = _internalnorm(scales)

    M_pl² = config.M_pl^2
    params = _BackgroundParams(potential, config.M_pl, M_pl²)
    N_end_requested = config.N_start + config.N_max

    u₀ = T[ϕ₀_T, ϕ_N₀_T]
    tspan = (config.N_start, N_end_requested)

    inflation_end_condition(u, N, integrator) =
        _ε_H(u[2], integrator.p.M_pl²) - one(T)
    inflation_end_callback = ContinuousCallback(
        inflation_end_condition,
        terminate!;
        affect_neg! = nothing,
        interp_points = 64,
        rootfind = RightRootFind,
    )

    problem = ODEProblem(_background_rhs!, u₀, tspan, params)
    solution = solve(
        problem,
        _solver_algorithm(config.solver);
        callback = inflation_end_callback,
        abstol = abstol_vec,
        reltol = config.rel_tol,
        internalnorm = internalnorm,
        dense = true,
        save_everystep = true,
        maxiters = config.max_iters,
        unstable_check = (dt, u, p, N) -> _domain_violation(u, p),
        _solver_kwargs(config.dt_min, config.dt_max)...,
    )

    termination = _termination_symbol(solution, N_end_requested, params, T)
    N_end = termination === :inflation_end ?
        _polished_inflation_end(solution, params, config) :
        T(solution.t[end])

    BackgroundSolution(solution, potential, config, N_end, termination)
end

@inline N_span(bg::BackgroundSolution) = (bg.config.N_start, bg.N_end)
@inline inflation_ended(bg::BackgroundSolution) = bg.termination === :inflation_end
@inline termination_reason(bg::BackgroundSolution) = bg.termination
