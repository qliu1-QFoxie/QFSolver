Base.@kwdef struct ClassPrimordialExport{T<:AbstractFloat,M}
    k_mpc::Vector{T}
    P_s::Vector{T}
    P_t::Union{Nothing,Vector{T}} = nothing
    metadata::M
end

function _validated_class_export_case(case::ExactEndpointCaseOutput)
    case_result = case.spectrum_case
    case_result isa ExactEndpointCaseResult || throw(ArgumentError(
        "CLASS export requires an ExactEndpointCaseOutput carrying an ExactEndpointCaseResult"))
    _validate_case_alignment(
        case_result.k_grid,
        case_result.scalar_modes,
        case_result.tensor_modes,
    )

    length(case.k_grid) == length(case_result.k_grid) || throw(ArgumentError(
        "CLASS export requires the exported physical k grid to align with the solved exact spectrum length"))

    metadata = case.metadata
    hasproperty(metadata, :k_axis_type) || throw(ArgumentError(
        "CLASS export requires k-axis metadata with k_axis_type"))
    metadata.k_axis_type == "physical_mpc^-1" || throw(ArgumentError(
        "CLASS export requires a calibrated physical k axis in Mpc^-1; got k_axis_type=$(metadata.k_axis_type)"))

    case_result
end

function _class_export_sortperm(k_grid::AbstractVector{<:Real})
    length(k_grid) >= 2 || throw(ArgumentError(
        "CLASS export requires at least two physical k points"))

    T = promote_type((_float_type(typeof(k)) for k in k_grid)...)
    T <: AbstractFloat || throw(ArgumentError(
        "CLASS export requires a floating-point-compatible physical k grid"))

    k_values = T[_checked_positive_finite(k, "case.k_grid[$index]", T) for (index, k) in enumerate(k_grid)]
    permutation = sortperm(k_values)

    previous_k = nothing
    for rank in eachindex(permutation)
        sorted_index = permutation[rank]
        k_value = k_values[sorted_index]
        if previous_k !== nothing && !(previous_k < k_value)
            throw(ArgumentError(
                "CLASS export physical k grid must be strictly increasing after sorting; duplicate or invalid entries detected near k=$(k_value)"))
        end
        previous_k = k_value
    end

    return k_values, permutation
end

@inline function _class_export_observable_name(expected_kind::PerturbationKind)
    expected_kind == scalar && return :P_s
    expected_kind == tensor && return :P_t
    throw(ArgumentError("unsupported perturbation kind for CLASS export"))
end

@inline function _class_export_observable_label(expected_kind::PerturbationKind)
    expected_kind == scalar && return "scalar"
    expected_kind == tensor && return "tensor"
    throw(ArgumentError("unsupported perturbation kind for CLASS export"))
end

function _extract_class_export_observable(
    entry,
    expected_kind::PerturbationKind,
    index::Int,
    k_phys::Real,
    ::Type{T};
    strict::Bool,
) where {T<:AbstractFloat}
    label = _class_export_observable_label(expected_kind)
    expected_name = _class_export_observable_name(expected_kind)

    if entry isa ExactEndpointModeFailure
        throw(ArgumentError(
            "CLASS export requires a successful $(label) mode at row $(index) (k=$(k_phys) Mpc^-1); got failure_stage=$(entry.failure_stage), reason=$(entry.failure_reason)"))
    end

    entry isa ExactEndpointModeResult || throw(ArgumentError(
        "CLASS export requires ExactEndpointModeResult or ExactEndpointModeFailure entries; row $(index) $(label) entry has type $(typeof(entry))"))
    entry.kind == expected_kind || throw(ArgumentError(
        "CLASS export expected kind=$(expected_kind) at row $(index) (k=$(k_phys) Mpc^-1); got kind=$(entry.kind)"))

    if strict && entry.observable_name != expected_name
        throw(ArgumentError(
            "CLASS export expected observable_name=$(expected_name) at row $(index) (k=$(k_phys) Mpc^-1); got $(entry.observable_name)"))
    end

    _checked_positive_finite(entry.observable_value, "$(label) observable at row $(index)", T)
end

function build_class_primordial_export(
    case::ExactEndpointCaseOutput;
    include_tensor::Bool=true,
    strict::Bool=true,
)
    case_result = _validated_class_export_case(case)
    k_values, permutation = _class_export_sortperm(case.k_grid)
    T = eltype(k_values)

    k_sorted = Vector{T}(undef, length(permutation))
    P_s_sorted = Vector{T}(undef, length(permutation))
    P_t_sorted = include_tensor ? Vector{T}(undef, length(permutation)) : nothing

    for export_index in eachindex(permutation)
        source_index = permutation[export_index]
        k_phys = k_values[source_index]
        k_sorted[export_index] = k_phys
        P_s_sorted[export_index] = _extract_class_export_observable(
            case_result.scalar_modes[source_index],
            scalar,
            export_index,
            k_phys,
            T;
            strict=strict,
        )

        if include_tensor
            P_t_sorted[export_index] = _extract_class_export_observable(
                case_result.tensor_modes[source_index],
                tensor,
                export_index,
                k_phys,
                T;
                strict=strict,
            )
        end
    end

    ClassPrimordialExport(
        k_mpc=k_sorted,
        P_s=P_s_sorted,
        P_t=P_t_sorted,
        metadata=case.metadata,
    )
end

function _validate_class_export_lengths(class_export::ClassPrimordialExport)
    length(class_export.k_mpc) == length(class_export.P_s) || throw(ArgumentError(
        "CLASS export scalar payload length must match the physical k grid length"))

    if class_export.P_t !== nothing
        length(class_export.k_mpc) == length(class_export.P_t) || throw(ArgumentError(
            "CLASS export tensor payload length must match the physical k grid length"))
    end

    nothing
end

function _class_export_header(class_export::ClassPrimordialExport)
    class_export.P_t === nothing && return "# k_mpc P_s"
    "# k_mpc P_s P_t"
end

@inline function _class_scalar_header()
    "# k_mpc P_s"
end

@inline function _class_tensor_header()
    "# k_mpc P_t"
end

function _write_class_primordial_rows(io::IO, class_export::ClassPrimordialExport)
    _validate_class_export_lengths(class_export)

    if class_export.P_t === nothing
        for index in eachindex(class_export.k_mpc)
            println(io, string(class_export.k_mpc[index], " ", class_export.P_s[index]))
        end
        return nothing
    end

    for index in eachindex(class_export.k_mpc)
        println(io, string(class_export.k_mpc[index], " ", class_export.P_s[index], " ", class_export.P_t[index]))
    end

    nothing
end

function write_class_primordial_table(
    path::AbstractString,
    class_export::ClassPrimordialExport;
    header::Bool=false,
)
    output_path = abspath(path)
    open(output_path, "w") do io
        header && println(io, _class_export_header(class_export))
        _write_class_primordial_rows(io, class_export)
    end
    output_path
end

function write_class_combined_table(
    path::AbstractString,
    class_export::ClassPrimordialExport;
    header::Bool=false,
)
    write_class_primordial_table(path, class_export; header=header)
end

function write_class_scalar_table(
    path::AbstractString,
    class_export::ClassPrimordialExport;
    header::Bool=false,
)
    _validate_class_export_lengths(class_export)
    output_path = abspath(path)
    open(output_path, "w") do io
        header && println(io, _class_scalar_header())
        for index in eachindex(class_export.k_mpc)
            println(io, string(class_export.k_mpc[index], " ", class_export.P_s[index]))
        end
    end
    output_path
end

function write_class_tensor_table(
    path::AbstractString,
    class_export::ClassPrimordialExport;
    header::Bool=false,
)
    _validate_class_export_lengths(class_export)
    class_export.P_t === nothing && throw(ArgumentError(
        "tensor table export requires a ClassPrimordialExport with tensor data"))

    output_path = abspath(path)
    open(output_path, "w") do io
        header && println(io, _class_tensor_header())
        for index in eachindex(class_export.k_mpc)
            println(io, string(class_export.k_mpc[index], " ", class_export.P_t[index]))
        end
    end
    output_path
end

function write_class_primordial_tables(
    prefix::AbstractString,
    class_export::ClassPrimordialExport;
    header::Bool=false,
)
    scalar_path = write_class_scalar_table(string(prefix, "_scalar.dat"), class_export; header=header)
    tensor_path = class_export.P_t === nothing ? nothing :
        write_class_tensor_table(string(prefix, "_tensor.dat"), class_export; header=header)

    (scalar_path=scalar_path, tensor_path=tensor_path)
end
