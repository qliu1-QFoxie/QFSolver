function solve_primordial_pipeline(req::PrimordialPipelineRequest)
    _validate_physical_k_request(req.k_request)
    bg = solve_background(
        req.potential,
        req.ϕ₀,
        req.ϕ_N₀;
        config=req.background_config,
    )
    tau_profile = build_τ(bg)

    mode_config = _resolved_mode_config(typeof(bg.config.M_pl), req.mode_config)
    endpoint_config = _resolved_endpoint_config(typeof(bg.config.M_pl), req.endpoint_config)
    calibration_config = _resolved_workflow_k_calibration_config(typeof(bg.config.M_pl), req.calibration_config)
    calibration = compute_k_calibration(bg, calibration_config)
    physical_k_grid, internal_k_grid = _resolve_requested_k_grids(calibration, req.k_request)
    exact_case = solve_exact_endpoint_case(
        bg,
        internal_k_grid;
        mode_config=mode_config,
        endpoint_config=endpoint_config,
        tau_profile=tau_profile,
    )
    exported_case = _export_exact_endpoint_case_with_calibration(
        exact_case,
        calibration;
        requested_physical_k_grid=physical_k_grid,
    )

    PrimordialPipelineResult(
        request=req,
        background=bg,
        tau_profile=tau_profile,
        spectrum_case=exported_case,
    )
end
