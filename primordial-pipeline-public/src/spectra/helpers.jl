@inline _float_type(::Type{T}) where {T<:Real} = float(T)

@inline function _checked_nonnegative_finite(
    value::Real,
    name::AbstractString,
    ::Type{T},
) where {T<:AbstractFloat}
    converted = _checked_finite(value, name, T)
    converted >= zero(T) || throw(ArgumentError("$(name) must be non-negative"))
    converted
end

function _validate_requested_k_interval(k_min_req::Real, k_max_req::Real)
    T = promote_type(_float_type(typeof(k_min_req)), _float_type(typeof(k_max_req)))
    T <: AbstractFloat || throw(ArgumentError(
        "requested k interval must be floating-point-compatible"))
    k_min_value = _checked_positive_finite(k_min_req, "k_min_req", T)
    k_max_value = _checked_positive_finite(k_max_req, "k_max_req", T)
    k_min_value < k_max_value || throw(ArgumentError(
        "requested k interval must satisfy k_min_req < k_max_req"))
    (k_min_value, k_max_value)
end

function _validate_dense_ppd(ppd::Real)
    _checked_positive_finite(ppd, "ppd", Float64)
end

function build_exact_logk_grid(k_min_req::Real, k_max_req::Real; ppd::Real)
    k_min_value, k_max_value = _validate_requested_k_interval(k_min_req, k_max_req)
    ppd_value = _validate_dense_ppd(ppd)
    n_k = ceil(Int, ppd_value * log10(k_max_value / k_min_value)) + 1
    n_k >= 2 || throw(ArgumentError(
        "exact dense log-k grid builder requires at least two points"))
    collect(exp.(range(log(k_min_value), log(k_max_value); length=n_k)))
end

function _validate_requested_k_grid(k_grid::AbstractVector{<:Real})
    length(k_grid) >= 2 || throw(ArgumentError(
        "requested k grid must contain at least two points"))

    K = promote_type((_float_type(typeof(k)) for k in k_grid)...)
    K <: AbstractFloat || throw(ArgumentError(
        "requested k grid must be floating-point-compatible"))

    previous_k = nothing
    for (index, k) in enumerate(k_grid)
        k_value = _checked_positive_finite(k, "k_grid[$index]", K)
        if previous_k !== nothing && !(previous_k < k_value)
            throw(ArgumentError(
                "requested k grid must be strictly increasing; entries $(index - 1) and $(index) violate monotonicity"))
        end
        previous_k = k_value
    end

    nothing
end

function _validate_case_alignment(k_grid, scalar_modes, tensor_modes)
    length(scalar_modes) == length(k_grid) || throw(ArgumentError(
        "scalar mode array length must match requested k grid length"))
    length(tensor_modes) == length(k_grid) || throw(ArgumentError(
        "tensor mode array length must match requested k grid length"))

    for index in eachindex(k_grid)
        scalar_mode = scalar_modes[index]
        tensor_mode = tensor_modes[index]

        if scalar_mode !== nothing
            getproperty(scalar_mode, :kind) == scalar || throw(ArgumentError(
                "scalar mode array entry $(index) must have kind=scalar"))
            getproperty(scalar_mode, :k) == k_grid[index] || throw(ArgumentError(
                "scalar mode array entry $(index) is not aligned with the requested k grid"))
        end

        if tensor_mode !== nothing
            getproperty(tensor_mode, :kind) == tensor || throw(ArgumentError(
                "tensor mode array entry $(index) must have kind=tensor"))
            getproperty(tensor_mode, :k) == k_grid[index] || throw(ArgumentError(
                "tensor mode array entry $(index) is not aligned with the requested k grid"))
        end
    end

    nothing
end
