# ECOS.jl

| **Build Status** |
|:----------------:|
| [![Build Status][build-img]][build-url] [![Build Status][winbuild-img]][winbuild-url] |
| [![Coveralls branch][coveralls-img]][coveralls-url] [![Codecov branch][codecov-img]][codecov-url] |

Julia wrapper for the [ECOS](https://github.com/embotech/ecos) embeddable conic optimization interior point solver.

## Installation

You can install ECOS.jl through the Julia package manager:
```julia
julia> Pkg.add("ECOS")
```

ECOS.jl will automatically install and setup the ECOS solver itself using [BinaryProvider.jl](https://github.com/JuliaPackaging/BinaryProvider.jl).

## Custom Installation

After ECOS.jl is installed and built, you can replace the installed binary dependencies with custom builds by overwritting the binaries and libraries in ECOS.jl's `deps/usr` folder (e.g. in Julia v0.6 `$HOME/.julia/v0.6/ECOS/deps/usr`).

Note that the custom binaries will not be overwritten by subsequent builds of the currently installed version of ECOS.jl. However, if ECOS.jl is updated and the update includes new BinaryProvider versions of the ECOS binaries, then the custom binaries will be overwritten by the new BinaryProvider versions.

## Usage

The ECOS interface is completely wrapped. ECOS functions corresponding to the C API are available as `ECOS.setup`, `ECOS.solve`, `ECOS.cleanup`, and `ECOS.ver` (these are not exported from the module). Function arguments are extensively documented in the source, and an example of usage can be found in `test/direct.jl`.

ECOS.jl also supports the **[MathOptInterface](https://github.com/JuliaOpt/MathOptInterface.jl)** standard solver interface.
Thanks to this support ECOS can be used as a solver with both the **[JuMP]** and **[Convex.jl]** modeling languages.

All ECOS solver options can be set through the direct interface and through MathOptInterface.
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
To use these settings you can either pass them as keyword arguments to `setup`
(direct interface) or as arguments to the `ECOS.Optimizer` constructor
(MathOptInterface interface), e.g.
```julia
# Direct
my_prob = ECOS.setup(n, m, ..., c, h, b; maxit=10, feastol=1e-5)
# MathOptInterface (with JuMP)
model = Model(with_optimizer(ECOS.Optimizer, maxit=10, feastol=1e-5))
```

### JuMP example

This example shows how we can model a simple knapsack problem with JuMP and use ECOS to solve it.

```julia
using JuMP
using ECOS

items  = [:Gold, :Silver, :Bronze]
values = Dict(:Gold => 5.0,  :Silver => 3.0,  :Bronze => 1.0)
weight = Dict(:Gold => 2.0,  :Silver => 1.5,  :Bronze => 0.3)

model = Model(with_optimizer(ECOS.Optimizer))
@variable(model, 0 <= take[items] <= 1)  # Define a variable for each item
@objective(model, Max, sum(values[item] * take[item] for item in items))
@constraint(model, sum(weight[item] * take[item] for item in items) <= 3)
optimize!(model)

println(value(take))
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

[build-img]: https://travis-ci.org/JuliaOpt/ECOS.jl.svg?branch=master
[build-url]: https://travis-ci.org/JuliaOpt/ECOS.jl
[winbuild-img]: https://ci.appveyor.com/api/projects/status/n0c8b6t1w39jho6d/branch/master?svg=true
[winbuild-url]: https://ci.appveyor.com/project/JuliaOpt/ecos-jl/branch/master
[coveralls-img]: https://coveralls.io/repos/github/JuliaOpt/ECOS.jl/badge.svg?branch=master
[coveralls-url]: https://coveralls.io/github/JuliaOpt/ECOS.jl?branch=master
[codecov-img]: http://codecov.io/github/JuliaOpt/ECOS.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/JuliaOpt/ECOS.jl?branch=master
