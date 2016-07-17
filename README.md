# ECOS.jl

[![Build Status](https://travis-ci.org/JuliaOpt/ECOS.jl.svg?branch=master)](https://travis-ci.org/JuliaOpt/ECOS.jl)
[![Coverage Status](https://img.shields.io/coveralls/JuliaOpt/ECOS.jl.svg)](https://coveralls.io/r/JuliaOpt/ECOS.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/bnvddmeevtrmjyc2/branch/master)](https://ci.appveyor.com/project/mlubin/ecos-jl/branch/master)

[![ECOS](http://pkg.julialang.org/badges/ECOS_0.4.svg)](http://pkg.julialang.org/?pkg=ECOS&ver=0.4)
[![ECOS](http://pkg.julialang.org/badges/ECOS_0.5.svg)](http://pkg.julialang.org/?pkg=ECOS&ver=0.5)

Julia wrapper for the [ECOS](https://github.com/embotech/ecos) embeddable conic optimization interior point solver.

## Installation

You can install ECOS.jl through the Julia package manager:
```julia
julia> Pkg.add("ECOS")
```

ECOS.jl will automatically setup the ECOS solver itself:
 - On Linux it will build from source
 - On OS X it will download a binary via [Homebrew.jl].
 - On Windows it will download a binary.

## Usage

The ECOS interface is completely wrapped. ECOS functions corresponding to the C API are available as `ECOS.setup`, `ECOS.solve`, `ECOS.cleanup`, and `ECOS.ver` (these are not exported from the module). Function arguments are extensively documented in the source, and an example of usage can be found in `test/direct.jl`.

ECOS.jl also supports the [JuliaOpt] **[MathProgBase]** standard solver interface.
Thanks to this support ECOS can be used as a solver with both the **[JuMP]** and **[Convex.jl]** modeling languages.

All ECOS solver options can be set through the direct interface and through MathProgBase.
The list of options is defined the [`ecos.h` header](https://github.com/embotech/ecos/blob/master/include/ecos.h), which we reproduce here:
```julia
gamma          # scaling the final step length
delta          # regularization parameter
eps            # regularization threshold
feastol        # primal/dual infeasibility tolerance
abstol         # absolute tolerance on duality gap
reltol         # relative tolerance on duality gap
feastol_inacc  # primal/dual infeasibility relaxed tolerance
abstol_inacc   # absolute relaxed tolerance on duality gap
reltol_inacc   # relative relaxed tolerance on duality gap
nitref         # number of iterative refinement steps
maxit          # maximum number of iterations
verbose        # verbosity bool for PRINTLEVEL < 3
```
To use these settings you can either pass them as keyword arguments to `setup` (direct interface) or as arguments to the `ECOSSolver` constructor (MathProgBase interface), e.g.
```julia
# Direct
my_prob = ECOS.setup(n, m, ..., c, h, b; maxit=10, feastol=1e-5)
# MathProgBase (with JuMP)
m = Model(solver=ECOS.ECOSSolver(maxit=10, feastol=1e-5))
```

### JuMP example

This example shows how we can model a simple knapsack problem with JuMP and use ECOS to solve it.

```julia
using JuMP
using ECOS

items  = [:Gold, :Silver, :Bronze]
values = Dict(:Gold => 5.0,  :Silver => 3.0,  :Bronze => 1.0)
weight = Dict(:Gold => 2.0,  :Silver => 1.5,  :Bronze => 0.3)

m = Model(solver=ECOSSolver())
@variable(m, 0 <= take[items] <= 1)  # Define a variable for each item
@objective(m, Max, sum{ values[item] * take[item], item in items})
@constraint(m, sum{ weight[item] * take[item], item in items} <= 3)
solve(m)

println(getvalue(take))
# take
# [  Gold] = 0.9999999680446406
# [Silver] = 0.46666670881026834
# [Bronze] = 0.9999999633898735
```

---

`ECOS.jl` is licensed under the MIT License (see LICENSE.md), but note that ECOS itself is GPL v3.

[MathProgBase]: https://github.com/JuliaOpt/MathProgBase.jl
[JuMP]: https://github.com/JuliaOpt/JuMP.jl
[Convex.jl]: https://github.com/JuliaOpt/Convex.jl
[Homebrew.jl]: https://github.com/JuliaLang/Homebrew.jl
[JuliaOpt]: http://juliaopt.org
