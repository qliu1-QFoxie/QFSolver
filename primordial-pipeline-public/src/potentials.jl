abstract type InflatonPotential end

function V end
function V_ϕ end
function V_ϕϕ end

struct QuadraticPotential{T<:AbstractFloat} <: InflatonPotential
    m::T
end

QuadraticPotential(m::Real) = begin
    m_value = float(m)
    QuadraticPotential{typeof(m_value)}(m_value)
end

V(pot::QuadraticPotential, ϕ) = pot.m^2 * ϕ^2 / 2
V_ϕ(pot::QuadraticPotential, ϕ) = pot.m^2 * ϕ
V_ϕϕ(pot::QuadraticPotential, ϕ) = pot.m^2
