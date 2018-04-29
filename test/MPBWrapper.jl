#############################################################################
# ECOS.jl
# Wrapper around the ECOS solver https://github.com/ifa-ethz/ecos
# See http://github.com/JuliaOpt/ECOS.jl
#############################################################################
# test/MPBWrapper.jl
# Test the MathProgBase.jl interface for the ECOS.jl solver wrapper
#############################################################################

@testset "Test the MPB wrapper with linprog" begin
    include("mpb_linear.jl")
end

import ECOS
@testset "Run the conic interface test from MathProgBase.jl" begin
    include(joinpath(Pkg.dir("MathProgBase"),"test","conicinterface.jl"))
    coniclineartest(ECOS.ECOSSolver(), duals=true)
    conicSOCtest(ECOS.ECOSSolver(), duals=true)
    conicEXPtest(ECOS.ECOSSolver(), duals=true)

    include(joinpath(Pkg.dir("MathProgBase"),"test","quadprog.jl"))
    socptest(ECOS.ECOSSolver())
end
