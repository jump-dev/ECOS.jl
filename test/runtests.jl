#############################################################################
# ECOS.jl
# Wrapper around the ECOS solver https://github.com/ifa-ethz/ecos
# See http://github.com/JuliaOpt/ECOS.jl
#############################################################################
# test/runtests.jl
# Test the ECOS.jl solver wrapper
#############################################################################

tests = ["direct.jl",
         "options.jl",
         "mpb_linear.jl"]

println("Running tests:")

for curtest in tests
    println(" Test: $(curtest)")
    include(curtest)
end

# Run the conic interface test from MathProgBase.jl
import ECOS
include(joinpath(Pkg.dir("MathProgBase"),"test","conicinterface.jl"))
coniclineartest(ECOS.ECOSSolver(), duals=true)
conicSOCtest(ECOS.ECOSSolver(), duals=true)
conicEXPtest(ECOS.ECOSSolver(), duals=true)

include(joinpath(Pkg.dir("MathProgBase"),"test","quadprog.jl"))
socptest(ECOS.ECOSSolver())
