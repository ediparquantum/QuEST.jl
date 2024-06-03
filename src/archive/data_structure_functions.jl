# QuEST_jl/src/base/data_structure_functions.jl
#

import .._QUEST_LIB

function createQuESTEnv() :: QuEST_Types.QuESTEnv

    if _QUEST_LIB[]==Ptr{Cvoid}(0)
        QuEST_init()
    else
        @assert 4*ccall(:getQuEST_PREC,Cint,(),) == sizeof(Qreal)
    end

    return ccall(:createQuESTEnv, QuEST_Types.QuESTEnv, () )

end

####################################################################################################

function createCloneQureg(qureg ::QuEST_Types.Qureg, env ::QuEST_Types.QuESTEnv) ::QuEST_Types.Qureg
    return ccall(:createCloneQureg, QuEST_Types.Qureg, (QuEST_Types.Qureg,QuEST_Types.QuESTEnv), qureg, env)
end

function createComplexMatrixN(numQubits ::Integer) ::QuEST_Types.ComplexMatrixN
    @assert 1 ≤ numQubits ≤ 50

    return ccall(:createComplexMatrixN, QuEST_Types.ComplexMatrixN,
                 (Cint,),
                 numQubits)
end

function createDensityQureg(numQubits ::Integer,
                            env       ::QuEST_Types.QuESTEnv) ::QuEST_Types.Qureg
    @assert 1 ≤ numQubits ≤ 50
    return ccall(:createDensityQureg, QuEST_Types.Qureg,
                 (Cint,     QuEST_Types.QuESTEnv),
                 numQubits, env)
end

function createDiagonalOp(numQubits     :: Integer,
                          env           :: QuEST_Types.QuESTEnv)    :: QuEST_Types.DiagonalOp
    return ccall(:createDiagonalOp, QuEST_Types.DiagonalOp,
                 (Cint,     QuEST_Types.QuESTEnv),
                 numQubits, env)
end

function createPauliHamil(numQubits     :: Integer,
                          numSumTerms   :: Integer)     :: QuEST_Types.PauliHamil
    return ccall(:createPauliHamil, QuEST_Types.PauliHamil,
                 (Cint, Cint),
                 numQubits, numSumTerms)
end


function createPauliHamilFromFile(filename ::String)   :: QuEST_Types.PauliHamil
    return ccall(:createPauliHamilFromFile, QuEST_Types.PauliHamil,
                 (Cstring,),
                 filename)
end

function createQureg(numQubits ::Integer,
                     env       ::QuEST_Types.QuESTEnv) ::QuEST_Types.Qureg
    @assert 1 ≤ numQubits ≤ 50

    return ccall(:createQureg, QuEST_Types.Qureg,
                 (Cint,     QuEST_Types.QuESTEnv),
                 numQubits, env)
end

function destroyComplexMatrixN(M ::QuEST_Types.ComplexMatrixN) ::Nothing
    ccall(:destroyComplexMatrixN, Cvoid,
          (QuEST_Types.ComplexMatrixN,),
          M)
    nothing
end

function destroyDiagonalOp(op       :: QuEST_Types.DiagonalOp,
                           env      :: QuEST_Types.QuESTEnv)    :: Nothing
    ccall(:destroyDiagonalOp, Cvoid,
          (QuEST_Types.DiagonalOp, QuEST_Types.QuESTEnv),
          op,                      env)
    return nothing
end

function destroyPauliHamil(hamil        ::QuEST_Types.PauliHamil)   :: Nothing
    ccall(:destroyPauliHamil, Cvoid, (QuEST_Types.PauliHamil,), hamil)
    return nothing
end

function destroyQuESTEnv(env ::QuEST_Types.QuESTEnv) ::Nothing
    ccall(:destroyQuESTEnv, Cvoid, (QuEST_Types.QuESTEnv,), env)
    return nothing
end

function destroyQureg(qureg ::QuEST_Types.Qureg, env ::QuEST_Types.QuESTEnv) ::Nothing
    ccall(:destroyQureg, Cvoid, (QuEST_Types.Qureg,QuEST_Types.QuESTEnv), qureg, env)
    return nothing
end

@doc raw"
Function `make_QuEST_matrix(M ::Matrix{Qreal}) ::ComplexMatrixN`
Convenience function for creating (`createComplexMatrixN()`) and filling a QuEST matrix data structure (`struct ComplexMatrixN`).
### Input
* `M` must be a 2ⁿ×2ⁿ matrix for 1 ≤ 𝑛 ≤ 50 (haha) [`@assert`]
### Return value
* An 𝑛 qubit matrix.
"
function make_QuEST_matrix(M ::Matrix{Base.Complex{Qreal}}) ::QuEST_Types.ComplexMatrixN
    (R,C) = size(M)
    @assert R==C
    @assert R ≥ 2
    numQubits = Int(round( log2(R) ))
    @assert 2^numQubits == R
    @assert 1 ≤ numQubits ≤ 50

    MQ = createComplexMatrixN(numQubits)
    for c = 1:R                       # Julia is column major
        for r = 1:R
            re_MQ_r = unsafe_load(MQ.real,r)
            im_MQ_r = unsafe_load(MQ.imag,r)
            unsafe_store!(re_MQ_r, real(M[r,c]), c)
            unsafe_store!(im_MQ_r, imag(M[r,c]), c)
        end
    end
    return MQ
end

@doc raw"
Function
     fill_ComplexMatrix!(M ::ComplexMatrixN, F ::Function) ::Nothing
Fills the entries of `M`: ∀(k,ℓ): 𝑀[k,ℓ] = 𝐹(k,ℓ).
(Indices are in 1, ..., 2ⁿ.)
### Input
* `M` must have been ''created'' (`createComplexMatrixN()`)
### Output
* Entries of `M` are overwritten.
"
function fill_ComplexMatrix!(M ::QuEST_Types.ComplexMatrixN, M_ ::Function) ::Nothing
    N = 2^M.numQubits
    for k = 1:N
        re_M_k = unsafe_load(M.real,k)
        im_M_k = unsafe_load(M.imag,k)
        for ℓ = 1:N
            unsafe_store!(re_M_k, real(M_(k,ℓ)), ℓ)
            unsafe_store!(im_M_k, imag(M_(k,ℓ)), ℓ)
        end
    end
    nothing
end

function initComplexMatrixN(m       ::QuEST_Types.ComplexMatrixN,
                            real_   ::Vector{Qreal},
                            imag_   ::Vector{Qreal})        :: Nothing

    @assert 1<<(2*m.numQubits) == length(real_) == length(imag_)
    ccall(:initComplexMatrixN, Cvoid, (QuEST_Types.ComplexMatrixN, Ptr{Qreal}, Ptr{Qreal}), m, real_, imag_)
    return nothing
end

function initDiagonalOp(op      :: QuEST_Types.DiagonalOp,
                        real_   :: Vector{Qreal},
                        imag_   :: Vector{Qreal})  ::Nothing

    @assert 1<<op.numQubits == length(real_) == length(imag_)
    ccall(:initDiagonalOp, Cvoid, (QuEST_Types.DiagonalOp, Ptr{Qreal}, Ptr{Qreal}), op, real_, imag_)
    return nothing
end

function initPauliHamil(hamil       :: QuEST_Types.PauliHamil,
                        coeffs      :: Vector{Qreal},
                        codes       :: Vector{QuEST_Types.pauliOpType})    ::Nothing

    @assert length(codes) == hamil.numSumTerms*hamil.numQubits
    ccall(:initPauliHamil, Cvoid, (QuEST_Types.PauliHamil, Ptr{Qreal}, Ptr{QuEST_Types.pauliOpType}), hamil, coeffs, codes)
    return nothing
end

function setDiagonalOpElems(op          ::QuEST_Types.DiagonalOp,
    startInd    ::Integer,
    real_       ::Vector{Qreal},
    imag_       ::Vector{Qreal},
    numElems    ::Integer)  ::Nothing
ccall(:setDiagonalOpElems,
Cvoid,
(QuEST_Types.DiagonalOp, Clonglong, Ptr{Qreal}, Ptr{Qreal}, Clonglong),
op,
Clonglong(startInd),
real_,
imag_,
Clonglong(numElems))
return nothing
end

function syncDiagonalOp(op        ::QuEST_Types.DiagonalOp)     ::Nothing
    ccall(:syncDiagonalOp, Cvoid, (QuEST_Types.DiagonalOp,), op)
    return nothing
end
