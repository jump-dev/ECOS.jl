# ECOS.jl

[![Build Status](https://github.com/jump-dev/ECOS.jl/workflows/CI/badge.svg?branch=master)](https://github.com/jump-dev/ECOS.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/jump-dev/ECOS.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jump-dev/ECOS.jl)

Julia wrapper for the [ECOS](https://github.com/embotech/ecos) embeddable conic
optimization interior point solver.

The wrapper has two components:
 * a thin wrapper around the complete C API
 * an iterface to [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl)

## Installation

Install ECOS.jl using `Pkg.add`:
```julia
import Pkg; Pkg.add("ECOS")
```

In addition to installing the ECOS.jl package, this will also download and
install the ECOS binaries. (You do not need to install ECOS separately.)

To use a custom binary, read the [Custom solver binaries](https://jump.dev/JuMP.jl/stable/developers/custom_solver_binaries/)
section of the JuMP documentation.

### License

`ECOS.jl` is licensed under the MIT License (see LICENSE.md), but note that ECOS
itself is GPL v3.

## Use with JuMP

TO use ECOS with [JuMP](https://github.com/jump-dev/JuMP.jl), use
`ECOS.Optimizer`:
```julia
using JuMP, ECOS
model = Model(ECOS.Optimizer)
set_optimizer_attribute(model, "maxit", 100)
```

## Options

The list of options is defined the [`ecos.h` header](https://github.com/embotech/ecos/blob/master/include/ecos.h),
which we reproduce here:
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
