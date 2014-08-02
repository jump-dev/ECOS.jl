# ECOS.jl

[![Build Status](https://travis-ci.org/JuliaOpt/ECOS.jl.svg?branch=master)](https://travis-ci.org/JuliaOpt/ECOS.jl)
[![Coverage Status](https://img.shields.io/coveralls/JuliaOpt/ECOS.jl.svg)](https://coveralls.io/r/JuliaOpt/ECOS.jl)

Julia wrapper for the [ECOS](https://github.com/ifa-ethz/ecos) embeddable second-order cone problem (SOCP) interior point solver.

## Installation

You can install ECOS.jl through the Julia package manager:
```julia
julia> Pkg.add("ECOS")
```

ECOS.jl will automatically setup the ECOS solver itself:
 - On Linux it will build from source
 - On OS X it will download a binary via [Homebrew.jl].
 - On Windows it will download a binary. [There is currently an issue with the 64-bit version of ECOS on Windows.](https://github.com/JuliaOpt/ECOS.jl/pull/8).

## Usage

The ECOS interface is completely wrapped: the package exports the functions `setup`, `solve`, and `cleanup`; it provides but does not export `ECOS_ver`. Function arguments are extensively documented in the source, and an example of usage can be found in `test/direct.jl`.

ECOS.jl also supports the [JuliaOpt] **[MathProgBase]** standard solver interface. This interface can be used to solve linear programs using `loadproblem!` (see `test/mpb_lin.jl`) and SOCPs through `loadconicproblem!` (see `test/mpb_conic.jl`). Thanks to this support ECOS can be used as a solver with both the **[JuMP]** and (in the very near future) **[CVX.jl]** modeling languages.

### JuMP example

This example shows how we can model a simple knapsack problem with JuMP and use ECOS to solve it.

```julia
using JuMP
import ECOS
m = Model(solver=ECOS.ECOSSolver())

items  = [:Gold, :Silver, :Bronze]
values = [:Gold => 5.0,  :Silver => 3.0,  :Bronze => 1.0]
weight = [:Gold => 2.0,  :Silver => 1.5,  :Bronze => 0.3]

@defVar(m, 0 <= take[items] <= 1)  # Define a variable for each item
@setObjective(m, Max, sum{ values[item] * take[item], item in items})
@addConstraint(m, sum{ weight[item] * take[item], item in items} <= 3)
solve(m)
println(getValue(take))
# take
# [  Gold] = 0.9999999680446406
# [Silver] = 0.46666670881026834
# [Bronze] = 0.9999999633898735
```

---

`ECOS.jl` is licensed under the MIT License (see LICENSE.md), but note that ECOS itself is GPL v3.

[MathProgBase]: https://github.com/JuliaOpt/MathProgBase.jl
[JuMP]: https://github.com/JuliaOpt/JuMP.jl
[CVX.jl]: https://github.com/cvxgrp/CVX.jl
[Homebrew.jl]: https://github.com/JuliaLang/Homebrew.jl
[JuliaOpt]: http://juliaopt.org
