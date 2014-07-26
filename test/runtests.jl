#############################################################################
# ECOS.jl
# Wrapper around the ECOS solver https://github.com/ifa-ethz/ecos
# See http://github.com/JuliaOpt/ECOS.jl
#############################################################################
# test/runtests.jl
# Test the ECOS.jl solver wrapper
#############################################################################

tests = ["direct.jl",
         "mpb.jl",
         "mpb_conic.jl"]

println("Running tests:")

for curtest in tests
    println(" Test: $(curtest)")
    include(curtest)
end