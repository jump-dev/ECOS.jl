# Copyright (c) 2014: ECOS.jl contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# Julia wrapper for header: cone.h
# Automatically generated using Clang.jl

function bring2cone(C, r, s)
    return ccall(
        (:bring2cone, ecos),
        Cvoid,
        (Ptr{cone}, Ptr{pfloat}, Ptr{pfloat}),
        C,
        r,
        s,
    )
end

function unitInitialization(C, s, z, scaling)
    return ccall(
        (:unitInitialization, ecos),
        Cvoid,
        (Ptr{cone}, Ptr{pfloat}, Ptr{pfloat}, pfloat),
        C,
        s,
        z,
        scaling,
    )
end

function updateScalings(C, s, z, lambda, mu)
    return ccall(
        (:updateScalings, ecos),
        idxint,
        (Ptr{cone}, Ptr{pfloat}, Ptr{pfloat}, Ptr{pfloat}, pfloat),
        C,
        s,
        z,
        lambda,
        mu,
    )
end

function evalSymmetricBarrierValue(siter, ziter, tauIter, kapIter, C, D)
    return ccall(
        (:evalSymmetricBarrierValue, ecos),
        pfloat,
        (Ptr{pfloat}, Ptr{pfloat}, pfloat, pfloat, Ptr{cone}, pfloat),
        siter,
        ziter,
        tauIter,
        kapIter,
        C,
        D,
    )
end

function scale(z, C, lambda)
    return ccall(
        (:scale, ecos),
        Cvoid,
        (Ptr{pfloat}, Ptr{cone}, Ptr{pfloat}),
        z,
        C,
        lambda,
    )
end

function scale2add(x, y, C)
    return ccall(
        (:scale2add, ecos),
        Cvoid,
        (Ptr{pfloat}, Ptr{pfloat}, Ptr{cone}),
        x,
        y,
        C,
    )
end

function unscale(lambda, C, z)
    return ccall(
        (:unscale, ecos),
        Cvoid,
        (Ptr{pfloat}, Ptr{cone}, Ptr{pfloat}),
        lambda,
        C,
        z,
    )
end

function conicProduct(u, v, C, w)
    return ccall(
        (:conicProduct, ecos),
        pfloat,
        (Ptr{pfloat}, Ptr{pfloat}, Ptr{cone}, Ptr{pfloat}),
        u,
        v,
        C,
        w,
    )
end

function conicDivision(u, v, C, w)
    return ccall(
        (:conicDivision, ecos),
        Cvoid,
        (Ptr{pfloat}, Ptr{pfloat}, Ptr{cone}, Ptr{pfloat}),
        u,
        v,
        C,
        w,
    )
end

function getSOCDetails(soc, conesize, eta_square, d1, u0, u1, v1, q)
    return ccall(
        (:getSOCDetails, ecos),
        Cvoid,
        (
            Ptr{socone},
            Ptr{idxint},
            Ptr{pfloat},
            Ptr{pfloat},
            Ptr{pfloat},
            Ptr{pfloat},
            Ptr{pfloat},
            Ptr{Ptr{pfloat}},
        ),
        soc,
        conesize,
        eta_square,
        d1,
        u0,
        u1,
        v1,
        q,
    )
end

function unstretch(n, p, C, Pinv, Px, dx, dy, dz)
    return ccall(
        (:unstretch, ecos),
        Cvoid,
        (
            idxint,
            idxint,
            Ptr{cone},
            Ptr{idxint},
            Ptr{pfloat},
            Ptr{pfloat},
            Ptr{pfloat},
            Ptr{pfloat},
        ),
        n,
        p,
        C,
        Pinv,
        Px,
        dx,
        dy,
        dz,
    )
end
# Julia wrapper for header: ctrlc.h
# Automatically generated using Clang.jl

# Julia wrapper for header: data.h
# Automatically generated using Clang.jl

# Julia wrapper for header: ecos.h
# Automatically generated using Clang.jl

function ECOS_setup(
    n,
    m,
    p,
    l,
    ncones,
    q,
    nex,
    Gpr,
    Gjc,
    Gir,
    Apr,
    Ajc,
    Air,
    c,
    h,
    b,
)
    return ccall(
        (:ECOS_setup, ecos),
        Ptr{pwork},
        (
            idxint,
            idxint,
            idxint,
            idxint,
            idxint,
            Ptr{idxint},
            idxint,
            Ptr{pfloat},
            Ptr{idxint},
            Ptr{idxint},
            Ptr{pfloat},
            Ptr{idxint},
            Ptr{idxint},
            Ptr{pfloat},
            Ptr{pfloat},
            Ptr{pfloat},
        ),
        n,
        m,
        p,
        l,
        ncones,
        q,
        nex,
        Gpr,
        Gjc,
        Gir,
        Apr,
        Ajc,
        Air,
        c,
        h,
        b,
    )
end

function expConeLineSearch(w, dtau, dkappa, affine)
    return ccall(
        (:expConeLineSearch, ecos),
        pfloat,
        (Ptr{pwork}, pfloat, pfloat, idxint),
        w,
        dtau,
        dkappa,
        affine,
    )
end

function ECOS_solve(w)
    return ccall((:ECOS_solve, ecos), idxint, (Ptr{pwork},), w)
end

function ECOS_cleanup(w, keepvars)
    return ccall(
        (:ECOS_cleanup, ecos),
        Cvoid,
        (Ptr{pwork}, idxint),
        w,
        keepvars,
    )
end

function ECOS_ver()
    return ccall((:ECOS_ver, ecos), Cstring, ())
end

function ecos_updateDataEntry_h(w, idx, value)
    return ccall(
        (:ecos_updateDataEntry_h, ecos),
        Cvoid,
        (Ptr{pwork}, idxint, pfloat),
        w,
        idx,
        value,
    )
end

function ecos_updateDataEntry_c(w, idx, value)
    return ccall(
        (:ecos_updateDataEntry_c, ecos),
        Cvoid,
        (Ptr{pwork}, idxint, pfloat),
        w,
        idx,
        value,
    )
end

function ECOS_updateData(w, Gpr, Apr, c, h, b)
    return ccall(
        (:ECOS_updateData, ecos),
        Cvoid,
        (
            Ptr{pwork},
            Ptr{pfloat},
            Ptr{pfloat},
            Ptr{pfloat},
            Ptr{pfloat},
            Ptr{pfloat},
        ),
        w,
        Gpr,
        Apr,
        c,
        h,
        b,
    )
end
# Julia wrapper for header: ecos_bb.h
# Automatically generated using Clang.jl

function ECOS_BB_setup(
    n,
    m,
    p,
    l,
    ncones,
    q,
    nex,
    Gpr,
    Gjc,
    Gir,
    Apr,
    Ajc,
    Air,
    c,
    h,
    b,
    num_bool_vars,
    bool_vars_idx,
    num_int_vars,
    int_vars_idx,
    stgs,
)
    return ccall(
        (:ECOS_BB_setup, ecos),
        Ptr{ecos_bb_pwork},
        (
            idxint,
            idxint,
            idxint,
            idxint,
            idxint,
            Ptr{idxint},
            idxint,
            Ptr{pfloat},
            Ptr{idxint},
            Ptr{idxint},
            Ptr{pfloat},
            Ptr{idxint},
            Ptr{idxint},
            Ptr{pfloat},
            Ptr{pfloat},
            Ptr{pfloat},
            idxint,
            Ptr{idxint},
            idxint,
            Ptr{idxint},
            Ptr{settings_bb},
        ),
        n,
        m,
        p,
        l,
        ncones,
        q,
        nex,
        Gpr,
        Gjc,
        Gir,
        Apr,
        Ajc,
        Air,
        c,
        h,
        b,
        num_bool_vars,
        bool_vars_idx,
        num_int_vars,
        int_vars_idx,
        stgs,
    )
end

function ECOS_BB_solve(prob)
    return ccall((:ECOS_BB_solve, ecos), idxint, (Ptr{ecos_bb_pwork},), prob)
end

function ECOS_BB_cleanup(prob, num_vars_keep)
    return ccall(
        (:ECOS_BB_cleanup, ecos),
        Cvoid,
        (Ptr{ecos_bb_pwork}, idxint),
        prob,
        num_vars_keep,
    )
end

function updateDataEntry_h(w, idx, value)
    return ccall(
        (:updateDataEntry_h, ecos),
        Cvoid,
        (Ptr{ecos_bb_pwork}, idxint, pfloat),
        w,
        idx,
        value,
    )
end

function updateDataEntry_c(w, idx, value)
    return ccall(
        (:updateDataEntry_c, ecos),
        Cvoid,
        (Ptr{ecos_bb_pwork}, idxint, pfloat),
        w,
        idx,
        value,
    )
end

function get_default_ECOS_BB_settings()
    return ccall((:get_default_ECOS_BB_settings, ecos), Ptr{settings_bb}, ())
end

function get_bool_node_id(idx, prob)
    return ccall(
        (:get_bool_node_id, ecos),
        Cstring,
        (idxint, Ptr{ecos_bb_pwork}),
        idx,
        prob,
    )
end

function get_int_node_id(idx, prob)
    return ccall(
        (:get_int_node_id, ecos),
        Ptr{pfloat},
        (idxint, Ptr{ecos_bb_pwork}),
        idx,
        prob,
    )
end

function abs_2(number)
    return ccall((:abs_2, ecos), pfloat, (pfloat,), number)
end

function pfloat_round(number)
    return ccall((:pfloat_round, ecos), pfloat, (pfloat,), number)
end

function pfloat_ceil(number, integer_tol)
    return ccall(
        (:pfloat_ceil, ecos),
        pfloat,
        (pfloat, pfloat),
        number,
        integer_tol,
    )
end

function pfloat_floor(number, integer_tol)
    return ccall(
        (:pfloat_floor, ecos),
        pfloat,
        (pfloat, pfloat),
        number,
        integer_tol,
    )
end

function float_eqls(a, b, integer_tol)
    return ccall(
        (:float_eqls, ecos),
        idxint,
        (pfloat, pfloat, pfloat),
        a,
        b,
        integer_tol,
    )
end
# Julia wrapper for header: equil.h
# Automatically generated using Clang.jl

function set_equilibration(w)
    return ccall((:set_equilibration, ecos), Cvoid, (Ptr{pwork},), w)
end

function unset_equilibration(w)
    return ccall((:unset_equilibration, ecos), Cvoid, (Ptr{pwork},), w)
end
# Julia wrapper for header: expcone.h
# Automatically generated using Clang.jl

function evalExpHessian(w, v, mu)
    return ccall(
        (:evalExpHessian, ecos),
        Cvoid,
        (Ptr{pfloat}, Ptr{pfloat}, pfloat),
        w,
        v,
        mu,
    )
end

function evalExpGradient(w, g)
    return ccall(
        (:evalExpGradient, ecos),
        Cvoid,
        (Ptr{pfloat}, Ptr{pfloat}),
        w,
        g,
    )
end

function evalBarrierValue(siter, ziter, fc, nexc)
    return ccall(
        (:evalBarrierValue, ecos),
        pfloat,
        (Ptr{pfloat}, Ptr{pfloat}, idxint, idxint),
        siter,
        ziter,
        fc,
        nexc,
    )
end

function scaleToAddExpcone(y, x, expcones, nexc, fc)
    return ccall(
        (:scaleToAddExpcone, ecos),
        Cvoid,
        (Ptr{pfloat}, Ptr{pfloat}, Ptr{expcone}, idxint, idxint),
        y,
        x,
        expcones,
        nexc,
        fc,
    )
end

function evalExpPrimalFeas(s, nexc)
    return ccall(
        (:evalExpPrimalFeas, ecos),
        idxint,
        (Ptr{pfloat}, idxint),
        s,
        nexc,
    )
end

function evalExpDualFeas(s, nexc)
    return ccall(
        (:evalExpDualFeas, ecos),
        idxint,
        (Ptr{pfloat}, idxint),
        s,
        nexc,
    )
end
# Julia wrapper for header: glblopts.h
# Automatically generated using Clang.jl

# Julia wrapper for header: kkt.h
# Automatically generated using Clang.jl

function kkt_factor(KKT, eps, delta)
    return ccall(
        (:kkt_factor, ecos),
        idxint,
        (Ptr{kkt}, pfloat, pfloat),
        KKT,
        eps,
        delta,
    )
end

function kkt_solve(KKT, A, G, Pb, dx, dy, dz, n, p, m, C, isinit, nitref)
    return ccall(
        (:kkt_solve, ecos),
        idxint,
        (
            Ptr{kkt},
            Ptr{spmat},
            Ptr{spmat},
            Ptr{pfloat},
            Ptr{pfloat},
            Ptr{pfloat},
            Ptr{pfloat},
            idxint,
            idxint,
            idxint,
            Ptr{cone},
            idxint,
            idxint,
        ),
        KKT,
        A,
        G,
        Pb,
        dx,
        dy,
        dz,
        n,
        p,
        m,
        C,
        isinit,
        nitref,
    )
end

function kkt_update(PKP, P, C)
    return ccall(
        (:kkt_update, ecos),
        Cvoid,
        (Ptr{spmat}, Ptr{idxint}, Ptr{cone}),
        PKP,
        P,
        C,
    )
end

function kkt_init(PKP, P, C)
    return ccall(
        (:kkt_init, ecos),
        Cvoid,
        (Ptr{spmat}, Ptr{idxint}, Ptr{cone}),
        PKP,
        P,
        C,
    )
end
# Julia wrapper for header: spla.h
# Automatically generated using Clang.jl

function sparseMV(A, x, y, a, newVector)
    return ccall(
        (:sparseMV, ecos),
        Cvoid,
        (Ptr{spmat}, Ptr{pfloat}, Ptr{pfloat}, idxint, idxint),
        A,
        x,
        y,
        a,
        newVector,
    )
end

function sparseMtVm(A, x, y, newVector, skipDiagonal)
    return ccall(
        (:sparseMtVm, ecos),
        Cvoid,
        (Ptr{spmat}, Ptr{pfloat}, Ptr{pfloat}, idxint, idxint),
        A,
        x,
        y,
        newVector,
        skipDiagonal,
    )
end

function vadd(n, x, y)
    return ccall(
        (:vadd, ecos),
        Cvoid,
        (idxint, Ptr{pfloat}, Ptr{pfloat}),
        n,
        x,
        y,
    )
end

function vsubscale(n, a, x, y)
    return ccall(
        (:vsubscale, ecos),
        Cvoid,
        (idxint, pfloat, Ptr{pfloat}, Ptr{pfloat}),
        n,
        a,
        x,
        y,
    )
end

function norm2(v, n)
    return ccall((:norm2, ecos), pfloat, (Ptr{pfloat}, idxint), v, n)
end

function norminf(v, n)
    return ccall((:norminf, ecos), pfloat, (Ptr{pfloat}, idxint), v, n)
end

function eddot(n, x, y)
    return ccall(
        (:eddot, ecos),
        pfloat,
        (idxint, Ptr{pfloat}, Ptr{pfloat}),
        n,
        x,
        y,
    )
end
# Julia wrapper for header: splamm.h
# Automatically generated using Clang.jl

function ecoscreateSparseMatrix(m, n, nnz, jc, ir, pr)
    return ccall(
        (:ecoscreateSparseMatrix, ecos),
        Ptr{spmat},
        (idxint, idxint, idxint, Ptr{idxint}, Ptr{idxint}, Ptr{pfloat}),
        m,
        n,
        nnz,
        jc,
        ir,
        pr,
    )
end

function newSparseMatrix(m, n, nnz)
    return ccall(
        (:newSparseMatrix, ecos),
        Ptr{spmat},
        (idxint, idxint, idxint),
        m,
        n,
        nnz,
    )
end

function freeSparseMatrix(M)
    return ccall((:freeSparseMatrix, ecos), Cvoid, (Ptr{spmat},), M)
end

function transposeSparseMatrix(M, MtoMt)
    return ccall(
        (:transposeSparseMatrix, ecos),
        Ptr{spmat},
        (Ptr{spmat}, Ptr{idxint}),
        M,
        MtoMt,
    )
end

function permuteSparseSymmetricMatrix(A, pinv, C, PK)
    return ccall(
        (:permuteSparseSymmetricMatrix, ecos),
        Cvoid,
        (Ptr{spmat}, Ptr{idxint}, Ptr{spmat}, Ptr{idxint}),
        A,
        pinv,
        C,
        PK,
    )
end

function pinv(n, p, pinv)
    return ccall(
        (:pinv, ecos),
        Cvoid,
        (idxint, Ptr{idxint}, Ptr{idxint}),
        n,
        p,
        pinv,
    )
end

function copySparseMatrix(A)
    return ccall((:copySparseMatrix, ecos), Ptr{spmat}, (Ptr{spmat},), A)
end

function printDenseMatrix(M, dim1, dim2, name)
    return ccall(
        (:printDenseMatrix, ecos),
        Cvoid,
        (Ptr{pfloat}, idxint, idxint, Cstring),
        M,
        dim1,
        dim2,
        name,
    )
end

function printDenseMatrix_i(M, dim1, dim2, name)
    return ccall(
        (:printDenseMatrix_i, ecos),
        Cvoid,
        (Ptr{idxint}, idxint, idxint, Cstring),
        M,
        dim1,
        dim2,
        name,
    )
end

function printSparseMatrix(M)
    return ccall((:printSparseMatrix, ecos), Cvoid, (Ptr{spmat},), M)
end

function dumpSparseMatrix(M, fn)
    return ccall((:dumpSparseMatrix, ecos), Cvoid, (Ptr{spmat}, Cstring), M, fn)
end

function dumpDenseMatrix(M, dim1, dim2, fn)
    return ccall(
        (:dumpDenseMatrix, ecos),
        Cvoid,
        (Ptr{pfloat}, Cint, Cint, Cstring),
        M,
        dim1,
        dim2,
        fn,
    )
end

function dumpDenseMatrix_i(M, dim1, dim2, fn)
    return ccall(
        (:dumpDenseMatrix_i, ecos),
        Cvoid,
        (Ptr{idxint}, Cint, Cint, Cstring),
        M,
        dim1,
        dim2,
        fn,
    )
end
# Julia wrapper for header: wright_omega.h
# Automatically generated using Clang.jl

function wrightOmega(z)
    return ccall((:wrightOmega, ecos), pfloat, (pfloat,), z)
end
