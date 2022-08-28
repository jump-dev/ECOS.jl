# Copyright (c) 2014: Joao Felipe Santos, Iain Dunning, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

using Test

@testset "Test the direct interface" begin
    include("c_wrapper.jl")
end

@testset "MathOptInterface" begin
    include("MOI_wrapper.jl")
end
