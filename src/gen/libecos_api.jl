# Julia wrapper for header: cone.h
# Automatically generated using Clang.jl


function bring2cone(C, r, s)
    ccall((:bring2cone, ecos), Cvoid, (Ptr{cone}, Ptr{pfloat}, Ptr{pfloat}), C, r, s)
end

function unitInitialization(C, s, z, scaling)
    ccall((:unitInitialization, ecos), Cvoid, (Ptr{cone}, Ptr{pfloat}, Ptr{pfloat}, pfloat), C, s, z, scaling)
end

function updateScalings(C, s, z, lambda, mu)
    ccall((:updateScalings, ecos), idxint, (Ptr{cone}, Ptr{pfloat}, Ptr{pfloat}, Ptr{pfloat}, pfloat), C, s, z, lambda, mu)
end

function evalSymmetricBarrierValue(siter, ziter, tauIter, kapIter, C, D)
    ccall((:evalSymmetricBarrierValue, ecos), pfloat, (Ptr{pfloat}, Ptr{pfloat}, pfloat, pfloat, Ptr{cone}, pfloat), siter, ziter, tauIter, kapIter, C, D)
end

function scale(z, C, lambda)
    ccall((:scale, ecos), Cvoid, (Ptr{pfloat}, Ptr{cone}, Ptr{pfloat}), z, C, lambda)
end

function scale2add(x, y, C)
    ccall((:scale2add, ecos), Cvoid, (Ptr{pfloat}, Ptr{pfloat}, Ptr{cone}), x, y, C)
end

function unscale(lambda, C, z)
    ccall((:unscale, ecos), Cvoid, (Ptr{pfloat}, Ptr{cone}, Ptr{pfloat}), lambda, C, z)
end

function conicProduct(u, v, C, w)
    ccall((:conicProduct, ecos), pfloat, (Ptr{pfloat}, Ptr{pfloat}, Ptr{cone}, Ptr{pfloat}), u, v, C, w)
end

function conicDivision(u, v, C, w)
    ccall((:conicDivision, ecos), Cvoid, (Ptr{pfloat}, Ptr{pfloat}, Ptr{cone}, Ptr{pfloat}), u, v, C, w)
end

function getSOCDetails(soc, conesize, eta_square, d1, u0, u1, v1, q)
    ccall((:getSOCDetails, ecos), Cvoid, (Ptr{socone}, Ptr{idxint}, Ptr{pfloat}, Ptr{pfloat}, Ptr{pfloat}, Ptr{pfloat}, Ptr{pfloat}, Ptr{Ptr{pfloat}}), soc, conesize, eta_square, d1, u0, u1, v1, q)
end

function unstretch(n, p, C, Pinv, Px, dx, dy, dz)
    ccall((:unstretch, ecos), Cvoid, (idxint, idxint, Ptr{cone}, Ptr{idxint}, Ptr{pfloat}, Ptr{pfloat}, Ptr{pfloat}, Ptr{pfloat}), n, p, C, Pinv, Px, dx, dy, dz)
end
# Julia wrapper for header: ctrlc.h
# Automatically generated using Clang.jl

# Julia wrapper for header: data.h
# Automatically generated using Clang.jl

# Julia wrapper for header: ecos.h
# Automatically generated using Clang.jl


function ECOS_setup(n, m, p, l, ncones, q, nex, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h, b)
    ccall((:ECOS_setup, ecos), Ptr{pwork}, (idxint, idxint, idxint, idxint, idxint, Ptr{idxint}, idxint, Ptr{pfloat}, Ptr{idxint}, Ptr{idxint}, Ptr{pfloat}, Ptr{idxint}, Ptr{idxint}, Ptr{pfloat}, Ptr{pfloat}, Ptr{pfloat}), n, m, p, l, ncones, q, nex, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h, b)
end

function expConeLineSearch(w, dtau, dkappa, affine)
    ccall((:expConeLineSearch, ecos), pfloat, (Ptr{pwork}, pfloat, pfloat, idxint), w, dtau, dkappa, affine)
end

function ECOS_solve(w)
    ccall((:ECOS_solve, ecos), idxint, (Ptr{pwork},), w)
end

function ECOS_cleanup(w, keepvars)
    ccall((:ECOS_cleanup, ecos), Cvoid, (Ptr{pwork}, idxint), w, keepvars)
end

function ECOS_ver()
    ccall((:ECOS_ver, ecos), Cstring, ())
end

function ecos_updateDataEntry_h(w, idx, value)
    ccall((:ecos_updateDataEntry_h, ecos), Cvoid, (Ptr{pwork}, idxint, pfloat), w, idx, value)
end

function ecos_updateDataEntry_c(w, idx, value)
    ccall((:ecos_updateDataEntry_c, ecos), Cvoid, (Ptr{pwork}, idxint, pfloat), w, idx, value)
end

function ECOS_updateData(w, Gpr, Apr, c, h, b)
    ccall((:ECOS_updateData, ecos), Cvoid, (Ptr{pwork}, Ptr{pfloat}, Ptr{pfloat}, Ptr{pfloat}, Ptr{pfloat}, Ptr{pfloat}), w, Gpr, Apr, c, h, b)
end
# Julia wrapper for header: ecos_bb.h
# Automatically generated using Clang.jl


function ECOS_BB_setup(n, m, p, l, ncones, q, nex, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h, b, num_bool_vars, bool_vars_idx, num_int_vars, int_vars_idx, stgs)
    ccall((:ECOS_BB_setup, ecos), Ptr{ecos_bb_pwork}, (idxint, idxint, idxint, idxint, idxint, Ptr{idxint}, idxint, Ptr{pfloat}, Ptr{idxint}, Ptr{idxint}, Ptr{pfloat}, Ptr{idxint}, Ptr{idxint}, Ptr{pfloat}, Ptr{pfloat}, Ptr{pfloat}, idxint, Ptr{idxint}, idxint, Ptr{idxint}, Ptr{settings_bb}), n, m, p, l, ncones, q, nex, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h, b, num_bool_vars, bool_vars_idx, num_int_vars, int_vars_idx, stgs)
end

function ECOS_BB_solve(prob)
    ccall((:ECOS_BB_solve, ecos), idxint, (Ptr{ecos_bb_pwork},), prob)
end

function ECOS_BB_cleanup(prob, num_vars_keep)
    ccall((:ECOS_BB_cleanup, ecos), Cvoid, (Ptr{ecos_bb_pwork}, idxint), prob, num_vars_keep)
end

function updateDataEntry_h(w, idx, value)
    ccall((:updateDataEntry_h, ecos), Cvoid, (Ptr{ecos_bb_pwork}, idxint, pfloat), w, idx, value)
end

function updateDataEntry_c(w, idx, value)
    ccall((:updateDataEntry_c, ecos), Cvoid, (Ptr{ecos_bb_pwork}, idxint, pfloat), w, idx, value)
end

function get_default_ECOS_BB_settings()
    ccall((:get_default_ECOS_BB_settings, ecos), Ptr{settings_bb}, ())
end

function get_bool_node_id(idx, prob)
    ccall((:get_bool_node_id, ecos), Cstring, (idxint, Ptr{ecos_bb_pwork}), idx, prob)
end

function get_int_node_id(idx, prob)
    ccall((:get_int_node_id, ecos), Ptr{pfloat}, (idxint, Ptr{ecos_bb_pwork}), idx, prob)
end

function abs_2(number)
    ccall((:abs_2, ecos), pfloat, (pfloat,), number)
end

function pfloat_round(number)
    ccall((:pfloat_round, ecos), pfloat, (pfloat,), number)
end

function pfloat_ceil(number, integer_tol)
    ccall((:pfloat_ceil, ecos), pfloat, (pfloat, pfloat), number, integer_tol)
end

function pfloat_floor(number, integer_tol)
    ccall((:pfloat_floor, ecos), pfloat, (pfloat, pfloat), number, integer_tol)
end

function float_eqls(a, b, integer_tol)
    ccall((:float_eqls, ecos), idxint, (pfloat, pfloat, pfloat), a, b, integer_tol)
end
# Julia wrapper for header: equil.h
# Automatically generated using Clang.jl


function set_equilibration(w)
    ccall((:set_equilibration, ecos), Cvoid, (Ptr{pwork},), w)
end

function unset_equilibration(w)
    ccall((:unset_equilibration, ecos), Cvoid, (Ptr{pwork},), w)
end
# Julia wrapper for header: expcone.h
# Automatically generated using Clang.jl


function evalExpHessian(w, v, mu)
    ccall((:evalExpHessian, ecos), Cvoid, (Ptr{pfloat}, Ptr{pfloat}, pfloat), w, v, mu)
end

function evalExpGradient(w, g)
    ccall((:evalExpGradient, ecos), Cvoid, (Ptr{pfloat}, Ptr{pfloat}), w, g)
end

function evalBarrierValue(siter, ziter, fc, nexc)
    ccall((:evalBarrierValue, ecos), pfloat, (Ptr{pfloat}, Ptr{pfloat}, idxint, idxint), siter, ziter, fc, nexc)
end

function scaleToAddExpcone(y, x, expcones, nexc, fc)
    ccall((:scaleToAddExpcone, ecos), Cvoid, (Ptr{pfloat}, Ptr{pfloat}, Ptr{expcone}, idxint, idxint), y, x, expcones, nexc, fc)
end

function evalExpPrimalFeas(s, nexc)
    ccall((:evalExpPrimalFeas, ecos), idxint, (Ptr{pfloat}, idxint), s, nexc)
end

function evalExpDualFeas(s, nexc)
    ccall((:evalExpDualFeas, ecos), idxint, (Ptr{pfloat}, idxint), s, nexc)
end
# Julia wrapper for header: glblopts.h
# Automatically generated using Clang.jl

# Julia wrapper for header: kkt.h
# Automatically generated using Clang.jl


function kkt_factor(KKT, eps, delta)
    ccall((:kkt_factor, ecos), idxint, (Ptr{kkt}, pfloat, pfloat), KKT, eps, delta)
end

function kkt_solve(KKT, A, G, Pb, dx, dy, dz, n, p, m, C, isinit, nitref)
    ccall((:kkt_solve, ecos), idxint, (Ptr{kkt}, Ptr{spmat}, Ptr{spmat}, Ptr{pfloat}, Ptr{pfloat}, Ptr{pfloat}, Ptr{pfloat}, idxint, idxint, idxint, Ptr{cone}, idxint, idxint), KKT, A, G, Pb, dx, dy, dz, n, p, m, C, isinit, nitref)
end

function kkt_update(PKP, P, C)
    ccall((:kkt_update, ecos), Cvoid, (Ptr{spmat}, Ptr{idxint}, Ptr{cone}), PKP, P, C)
end

function kkt_init(PKP, P, C)
    ccall((:kkt_init, ecos), Cvoid, (Ptr{spmat}, Ptr{idxint}, Ptr{cone}), PKP, P, C)
end
# Julia wrapper for header: spla.h
# Automatically generated using Clang.jl


function sparseMV(A, x, y, a, newVector)
    ccall((:sparseMV, ecos), Cvoid, (Ptr{spmat}, Ptr{pfloat}, Ptr{pfloat}, idxint, idxint), A, x, y, a, newVector)
end

function sparseMtVm(A, x, y, newVector, skipDiagonal)
    ccall((:sparseMtVm, ecos), Cvoid, (Ptr{spmat}, Ptr{pfloat}, Ptr{pfloat}, idxint, idxint), A, x, y, newVector, skipDiagonal)
end

function vadd(n, x, y)
    ccall((:vadd, ecos), Cvoid, (idxint, Ptr{pfloat}, Ptr{pfloat}), n, x, y)
end

function vsubscale(n, a, x, y)
    ccall((:vsubscale, ecos), Cvoid, (idxint, pfloat, Ptr{pfloat}, Ptr{pfloat}), n, a, x, y)
end

function norm2(v, n)
    ccall((:norm2, ecos), pfloat, (Ptr{pfloat}, idxint), v, n)
end

function norminf(v, n)
    ccall((:norminf, ecos), pfloat, (Ptr{pfloat}, idxint), v, n)
end

function eddot(n, x, y)
    ccall((:eddot, ecos), pfloat, (idxint, Ptr{pfloat}, Ptr{pfloat}), n, x, y)
end
# Julia wrapper for header: splamm.h
# Automatically generated using Clang.jl


function ecoscreateSparseMatrix(m, n, nnz, jc, ir, pr)
    ccall((:ecoscreateSparseMatrix, ecos), Ptr{spmat}, (idxint, idxint, idxint, Ptr{idxint}, Ptr{idxint}, Ptr{pfloat}), m, n, nnz, jc, ir, pr)
end

function newSparseMatrix(m, n, nnz)
    ccall((:newSparseMatrix, ecos), Ptr{spmat}, (idxint, idxint, idxint), m, n, nnz)
end

function freeSparseMatrix(M)
    ccall((:freeSparseMatrix, ecos), Cvoid, (Ptr{spmat},), M)
end

function transposeSparseMatrix(M, MtoMt)
    ccall((:transposeSparseMatrix, ecos), Ptr{spmat}, (Ptr{spmat}, Ptr{idxint}), M, MtoMt)
end

function permuteSparseSymmetricMatrix(A, pinv, C, PK)
    ccall((:permuteSparseSymmetricMatrix, ecos), Cvoid, (Ptr{spmat}, Ptr{idxint}, Ptr{spmat}, Ptr{idxint}), A, pinv, C, PK)
end

function pinv(n, p, pinv)
    ccall((:pinv, ecos), Cvoid, (idxint, Ptr{idxint}, Ptr{idxint}), n, p, pinv)
end

function copySparseMatrix(A)
    ccall((:copySparseMatrix, ecos), Ptr{spmat}, (Ptr{spmat},), A)
end

function printDenseMatrix(M, dim1, dim2, name)
    ccall((:printDenseMatrix, ecos), Cvoid, (Ptr{pfloat}, idxint, idxint, Cstring), M, dim1, dim2, name)
end

function printDenseMatrix_i(M, dim1, dim2, name)
    ccall((:printDenseMatrix_i, ecos), Cvoid, (Ptr{idxint}, idxint, idxint, Cstring), M, dim1, dim2, name)
end

function printSparseMatrix(M)
    ccall((:printSparseMatrix, ecos), Cvoid, (Ptr{spmat},), M)
end

function dumpSparseMatrix(M, fn)
    ccall((:dumpSparseMatrix, ecos), Cvoid, (Ptr{spmat}, Cstring), M, fn)
end

function dumpDenseMatrix(M, dim1, dim2, fn)
    ccall((:dumpDenseMatrix, ecos), Cvoid, (Ptr{pfloat}, Cint, Cint, Cstring), M, dim1, dim2, fn)
end

function dumpDenseMatrix_i(M, dim1, dim2, fn)
    ccall((:dumpDenseMatrix_i, ecos), Cvoid, (Ptr{idxint}, Cint, Cint, Cstring), M, dim1, dim2, fn)
end
# Julia wrapper for header: wright_omega.h
# Automatically generated using Clang.jl


function wrightOmega(z)
    ccall((:wrightOmega, ecos), pfloat, (pfloat,), z)
end
