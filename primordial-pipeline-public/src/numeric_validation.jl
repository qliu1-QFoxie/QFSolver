@inline function _checked_finite(
    value::Real,
    name::AbstractString,
    ::Type{T},
) where {T<:AbstractFloat}
    converted = T(value)
    isfinite(converted) || throw(ArgumentError("$(name) must be finite"))
    converted
end

@inline function _checked_positive_finite(
    value::Real,
    name::AbstractString,
    ::Type{T},
) where {T<:AbstractFloat}
    converted = _checked_finite(value, name, T)
    converted > zero(T) || throw(ArgumentError("$(name) must be strictly positive"))
    converted
end
