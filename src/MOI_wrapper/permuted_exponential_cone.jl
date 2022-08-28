# Copyright (c) 2014: Joao Felipe Santos, Iain Dunning, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# ECOS orders differently than MOI the second and third dimension of the
# exponential cone

struct PermutedExponentialCone <: MOI.AbstractVectorSet end

Base.copy(set::PermutedExponentialCone) = set
function MOI.Utilities.set_with_dimension(::Type{PermutedExponentialCone}, dim)
    if dim != 3
        error("Cannot create a `PermutedExponentialCone` with dimension $dim.")
    end
    return PermutedExponentialCone()
end

struct PermutedExponentialBridge{T,F} <: MOI.Bridges.Constraint.SetMapBridge{
    T,
    PermutedExponentialCone,
    MOI.ExponentialCone,
    F,
    F,
}
    constraint::MOI.ConstraintIndex{F,PermutedExponentialCone}
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:PermutedExponentialBridge{T}},
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{MOI.ExponentialCone},
) where {T}
    return PermutedExponentialBridge{T,F}
end

function MOI.Bridges.map_set(
    ::Type{<:PermutedExponentialBridge},
    ::MOI.ExponentialCone,
)
    return PermutedExponentialCone()
end

function MOI.Bridges.inverse_map_set(
    ::Type{<:PermutedExponentialBridge},
    ::PermutedExponentialCone,
)
    return MOI.ExponentialCone()
end

function MOI.Bridges.map_function(
    ::Type{<:PermutedExponentialBridge{T}},
    func,
) where {T}
    scalars = MOI.Utilities.eachscalar(func)
    return scalars[[1, 3, 2]]
end

# Swapping indices 2 <-> 3 is an involution (it is its own inverse)

function MOI.Bridges.inverse_map_function(
    BT::Type{<:PermutedExponentialBridge},
    func,
)
    return MOI.Bridges.map_function(BT, func)
end

# It is also symmetric (it is its own transpose)

function MOI.Bridges.adjoint_map_function(
    BT::Type{<:PermutedExponentialBridge},
    func,
)
    return MOI.Bridges.map_function(BT, func)
end

function MOI.Bridges.inverse_adjoint_map_function(
    BT::Type{<:PermutedExponentialBridge},
    func,
)
    return MOI.Bridges.map_function(BT, func)
end
