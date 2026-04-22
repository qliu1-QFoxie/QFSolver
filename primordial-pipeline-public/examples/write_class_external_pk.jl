using InflationPrimordialPipeline

const CLASS_EXTERNAL_PK_PP = InflationPrimordialPipeline

function _class_external_pk_usage()
    """
    Usage:
      julia --project=. examples/write_class_external_pk.jl \\
        --output primordial_external_pk.dat [--header] [--split-prefix PREFIX]

    Writes the canonical CLASS external_Pk primordial table for a representative
    calibrated quadratic-potential run.

    Options:
      --output PATH          Output path for the combined CLASS table.
      --header               Include a commented column header.
      --split-prefix PREFIX  Also write PREFIX_scalar.dat and PREFIX_tensor.dat.
      --k-min VALUE          Minimum physical k in Mpc^-1. Default: 0.03.
      --k-max VALUE          Maximum physical k in Mpc^-1. Default: 0.08.
      --ppd VALUE            Log-grid points per decade. Default: 16.
      -h, --help             Show this help text.
    """
end

function _require_class_external_pk_value(args, index, option)
    index < length(args) || throw(ArgumentError("$(option) requires a value"))
    args[index + 1]
end

function _parse_class_external_pk_args(args)
    output_path = nothing
    header = false
    split_prefix = nothing
    k_min = 0.03
    k_max = 0.08
    ppd = 16.0

    index = 1
    while index <= length(args)
        arg = args[index]

        if arg == "-h" || arg == "--help"
            println(_class_external_pk_usage())
            return nothing
        elseif arg == "--header"
            header = true
        elseif arg == "--output"
            output_path = _require_class_external_pk_value(args, index, arg)
            index += 1
        elseif startswith(arg, "--output=")
            output_path = arg[length("--output=") + 1:end]
        elseif arg == "--split-prefix"
            split_prefix = _require_class_external_pk_value(args, index, arg)
            index += 1
        elseif startswith(arg, "--split-prefix=")
            split_prefix = arg[length("--split-prefix=") + 1:end]
        elseif arg == "--k-min"
            k_min = parse(Float64, _require_class_external_pk_value(args, index, arg))
            index += 1
        elseif startswith(arg, "--k-min=")
            k_min = parse(Float64, arg[length("--k-min=") + 1:end])
        elseif arg == "--k-max"
            k_max = parse(Float64, _require_class_external_pk_value(args, index, arg))
            index += 1
        elseif startswith(arg, "--k-max=")
            k_max = parse(Float64, arg[length("--k-max=") + 1:end])
        elseif arg == "--ppd"
            ppd = parse(Float64, _require_class_external_pk_value(args, index, arg))
            index += 1
        elseif startswith(arg, "--ppd=")
            ppd = parse(Float64, arg[length("--ppd=") + 1:end])
        elseif startswith(arg, "-")
            throw(ArgumentError("unknown option: $(arg)"))
        elseif output_path === nothing
            output_path = arg
        else
            throw(ArgumentError("unexpected positional argument: $(arg)"))
        end

        index += 1
    end

    if output_path === nothing
        output_path = "primordial_external_pk.dat"
    end

    isfinite(k_min) && k_min > 0.0 || throw(ArgumentError("--k-min must be finite and positive"))
    isfinite(k_max) && k_max > k_min || throw(ArgumentError("--k-max must be finite and greater than --k-min"))
    isfinite(ppd) && ppd > 0.0 || throw(ArgumentError("--ppd must be finite and positive"))

    (
        output_path=String(output_path),
        header=header,
        split_prefix=split_prefix,
        k_min=k_min,
        k_max=k_max,
        ppd=ppd,
    )
end

function _representative_class_external_pk_result(; k_min::Float64, k_max::Float64, ppd::Float64)
    PP = CLASS_EXTERNAL_PK_PP
    potential = PP.QuadraticPotential(6.0e-6)
    phi_0 = 16.0
    phi_N_0 = -2.0 / phi_0

    PP.solve_exact_physical_spectrum(
        potential,
        phi_0,
        phi_N_0;
        background_config=PP.BackgroundConfig{Float64}(N_max=80.0),
        N_star=25.0,
        k_request=PP.ExplicitPhysicalKIntervalRequest(k_min, k_max, ppd),
        mode_config=PP.ModeConfig{Float64}(end_padding=0.0),
        endpoint_config=PP.ExactEndpointConfig{Float64}(),
    )
end

function run_class_external_pk_example(args=ARGS)
    options = _parse_class_external_pk_args(args)
    options === nothing && return nothing

    result = _representative_class_external_pk_result(
        k_min=options.k_min,
        k_max=options.k_max,
        ppd=options.ppd,
    )
    PP = CLASS_EXTERNAL_PK_PP
    class_export = PP.build_class_primordial_export(result)

    output_path = abspath(options.output_path)
    mkpath(dirname(output_path))
    table_path = PP.write_class_primordial_table(output_path, class_export; header=options.header)

    split_paths = nothing
    if options.split_prefix !== nothing
        split_prefix = abspath(options.split_prefix)
        mkpath(dirname(split_prefix))
        split_paths = PP.write_class_primordial_tables(split_prefix, class_export; header=options.header)
    end

    println("write_class_external_pk: ok")
    println("table = ", table_path)
    println("columns = ", class_export.P_t === nothing ? "k_mpc P_s" : "k_mpc P_s P_t")
    println("rows = ", length(class_export.k_mpc))
    println("k-range [Mpc^-1] = [", first(class_export.k_mpc), ", ", last(class_export.k_mpc), "]")

    if split_paths !== nothing
        println("scalar_table = ", split_paths.scalar_path)
        println("tensor_table = ", split_paths.tensor_path)
    end

    (table_path=table_path, split_paths=split_paths, class_export=class_export)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_class_external_pk_example(ARGS)
end
