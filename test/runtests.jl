#############################################################################
# ECOS.jl
# Wrapper around the ECOS solver https://github.com/ifa-ethz/ecos
# See http://github.com/JuliaOpt/ECOS.jl
#############################################################################
# test/runtests.jl
# Test the ECOS.jl solver wrapper
#############################################################################

using ECOS

using Compat
using Compat.Test

@testset "Test the direct interface" begin
    include("direct.jl")
end

@testset "Test passing options" begin
    include("options.jl")
end

@testset "MathProgBase" begin
    include("MPB_wrapper.jl")
end

@testset "MathOptInterface" begin
    include("MOI_wrapper.jl")
end
