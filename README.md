# ECOS.jl

[![Build Status](https://travis-ci.org/JuliaOpt/ECOS.jl.svg?branch=master)](https://travis-ci.org/JuliaOpt/ECOS.jl)

Julia wrapper for the [ECOS](https://github.com/ifa-ethz/ecos) second-order cone problem (SOCP) interior-point solver.

ECOS.jl will automatically setup the ECOS solver itself:
 - On Linux it will build from source
 - OS X and Windows it will download a binary.

The package exports the functions `setup`, `solve`, and `cleanup`. It provides but does not export `ECOS_ver`. Function arguments are extensively documented in the source, and an example usage is in `test/direct.jl`.

`ECOS.jl` also supports the **[MathProgBase]** standard interface. In particular it can be used to solve linear programs using `loadproblem!` (see `test/mpb_lin.jl`) and SOCPs through `loadconicproblem!` (see `test/mpb_conic.jl`). ECOS can be used as a solver with **[JuMP]** and in the near future **[CVX.jl]**.

`ECOS.jl` is licensed under the MIT License (see LICENSE.md) but ECOS itself is GPL v3.

[MathProgBase]: https://github.com/JuliaOpt/MathProgBase.jl
[JuMP]: https://github.com/JuliaOpt/JuMP.jl
[CVX.jl]: https://github.com/cvxgrp/CVX.jl
