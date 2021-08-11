# ECOS orders differently than MOI the second and third dimension of the exponential cone

struct PermutedExponentialCone <: MOI.AbstractVectorSet
end
Base.copy(set::PermutedExponentialCone) = set

const MOIB = MOI.Bridges

struct PermutedExponentialBridge{T,F} <:
       MOIB.Constraint.SetMapBridge{T,PermutedExponentialCone,MOI.ExponentialCone,F,F}
    constraint::MOI.ConstraintIndex{F,PermutedExponentialCone}
end

function MOIB.Constraint.concrete_bridge_type(
    ::Type{<:PermutedExponentialBridge{T}},
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{MOI.ExponentialCone},
) where {T}
    return PermutedExponentialBridge{T,F}
end

function MOIB.map_set(::Type{<:PermutedExponentialBridge}, ::MOI.ExponentialCone)
    return PermutedExponentialCone()
end

function MOIB.inverse_map_set(
    ::Type{<:PermutedExponentialBridge},
    ::PermutedExponentialCone,
)
    return MOI.ExponentialCone()
end

function MOIB.map_function(::Type{<:PermutedExponentialBridge{T}}, func) where {T}
    scalars = MOIU.eachscalar(func)
    return scalars[[1, 3, 2]]
end

# Swapping indices 2 <-> 3 is an involution (it is its own inverse)

function MOIB.inverse_map_function(BT::Type{<:PermutedExponentialBridge}, func)
    return MOIB.map_function(BT, func)
end

# It is also symmetric (it is its own transpose)

function MOIB.adjoint_map_function(BT::Type{<:PermutedExponentialBridge}, func)
    return MOIB.map_function(BT, func)
end

function MOIB.inverse_adjoint_map_function(
    BT::Type{<:PermutedExponentialBridge},
    func,
)
    return MOIB.map_function(BT, func)
end
