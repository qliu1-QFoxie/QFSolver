@inline function _τ_tolerances(
    bg::BackgroundSolution{T};
    abs_tol::Union{Nothing,Real}=nothing,
    rel_tol::Union{Nothing,Real}=nothing,
) where {T<:AbstractFloat}
    abs_tol_value = abs_tol === nothing ? bg.config.abs_tol / T(10) : T(abs_tol)
    rel_tol_value = rel_tol === nothing ? bg.config.rel_tol / T(10) : T(rel_tol)

    isfinite(abs_tol_value) && abs_tol_value > zero(T) || throw(ArgumentError(
        "τ quadrature abs_tol must be finite and positive"))
    isfinite(rel_tol_value) && rel_tol_value > zero(T) || throw(ArgumentError(
        "τ quadrature rel_tol must be finite and positive"))

    (abs_tol=abs_tol_value, rel_tol=rel_tol_value)
end

@inline function _τ_integrand(bg::BackgroundSolution{T}, N) where {T<:AbstractFloat}
    one(T) / aH(bg, N)
end

function _integrate_inverse_aH(
    bg::BackgroundSolution{T},
    N_lo::T,
    N_hi::T;
    abs_tol::T,
    rel_tol::T,
) where {T<:AbstractFloat}
    N_lo == N_hi && return zero(T)
    integral, _ = quadgk(
        N -> _τ_integrand(bg, N),
        N_lo,
        N_hi;
        atol=abs_tol,
        rtol=rel_tol,
    )
    T(integral)
end

@inline function _τ_end_value(::Type{T}, τ_end::Union{Nothing,Real}) where {T<:AbstractFloat}
    τ_end === nothing ? zero(T) : T(τ_end)
end

function _compressed_N_nodes(N_nodes::Vector{T}) where {T<:AbstractFloat}
    length(N_nodes) <= 1 && return copy(N_nodes)

    tol = _N_tolerance(T, first(N_nodes), last(N_nodes))
    compressed = T[N_nodes[1]]

    for N in Iterators.drop(N_nodes, 1)
        N < compressed[end] - tol && throw(ArgumentError(
            "Background ODE nodes must be nondecreasing to build τ"))
        abs(N - compressed[end]) <= tol && continue
        push!(compressed, N)
    end

    compressed
end

function _τ_nodes(bg::BackgroundSolution{T}) where {T<:AbstractFloat}
    raw_nodes = _compressed_N_nodes(T.(bg.ode_solution.t))
    isempty(raw_nodes) && throw(ArgumentError("Cannot build τ profile from an empty background solution"))

    N_end = bg.N_end
    tol = _N_tolerance(T, first(raw_nodes), N_end)
    nodes = T[]

    for N in raw_nodes
        N <= N_end + tol || continue
        push!(nodes, min(N, N_end))
    end

    isempty(nodes) && throw(ArgumentError("Could not construct τ nodes up to the polished background end time"))

    if abs(nodes[end] - N_end) <= tol
        nodes[end] = N_end
    else
        push!(nodes, N_end)
    end

    nodes
end

function build_τ(
    bg::BackgroundSolution{T};
    τ_end::Union{Nothing,Real}=nothing,
    abs_tol::Union{Nothing,Real}=nothing,
    rel_tol::Union{Nothing,Real}=nothing,
) where {T<:AbstractFloat}
    tolerances = _τ_tolerances(bg; abs_tol=abs_tol, rel_tol=rel_tol)
    τ_end_value = _τ_end_value(T, τ_end)
    isfinite(τ_end_value) || throw(ArgumentError("τ_end must be finite"))

    N_nodes = _τ_nodes(bg)

    τ_nodes = Vector{T}(undef, length(N_nodes))
    τ_nodes[end] = τ_end_value

    for i in (length(N_nodes) - 1):-1:1
        Δτ = _integrate_inverse_aH(
            bg,
            N_nodes[i],
            N_nodes[i + 1];
            abs_tol=tolerances.abs_tol,
            rel_tol=tolerances.rel_tol,
        )
        τ_nodes[i] = τ_nodes[i + 1] - Δτ
    end

    τProfile(
        bg,
        N_nodes,
        τ_nodes,
        τ_end_value,
        tolerances.abs_tol,
        tolerances.rel_tol,
    )
end

function τ(τp::τProfile{T}, N::Real) where {T<:AbstractFloat}
    N_value = _checked_N(τp.background, N)
    idx = searchsortedlast(τp.N_nodes, N_value)
    idx = clamp(idx, firstindex(τp.N_nodes), lastindex(τp.N_nodes))

    τp.N_nodes[idx] == N_value && return τp.τ_nodes[idx]
    idx == lastindex(τp.N_nodes) && return τp.τ_nodes[end]

    τp.τ_nodes[idx] + _integrate_inverse_aH(
        τp.background,
        τp.N_nodes[idx],
        N_value;
        abs_tol=τp.abs_tol,
        rel_tol=τp.rel_tol,
    )
end

@inline function τ_N(τp::τProfile{T}, N::Real) where {T<:AbstractFloat}
    one(T) / aH(τp.background, N)
end
