##################################################################
# Filename  : functions_with_issues.jl
# Author    : Jonathan Miller
# Date      : 2024-01-11
# Aim       : List tests that are not working and functions from tartU 
#             that were used and give credit
#           : check_ComplexMatrixN
#           : test_DiagonalOp
#           : test_controlledMultiQubitUnitary
#           :
##################################################################



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




function check_ComplexMatrixN()
    num_qubits = 3
    M = createComplexMatrixN(num_qubits)
    @test M.numQubits == num_qubits
    M_j_real = rand(Float64, 2^num_qubits * 2^num_qubits)
    M_j_imag = rand(Float64, 2^num_qubits * 2^num_qubits)
    initComplexMatrixN(M, M_j_real, M_j_imag)
        real_ptrs = unsafe_wrap(Vector{Ptr{qreal}}, M.real, 2^i)
        imag_ptrs = unsafe_wrap(Vector{Ptr{qreal}}, M.imag, 2^i)
        for r=1:2^i
            real_row = unsafe_wrap(Vector{qreal}, real_ptrs[r], 2^i)
            imag_row = unsafe_wrap(Vector{qreal}, imag_ptrs[r], 2^i)
            for c=1:2^i
                @test M_j_real[(r-1)*2^i + c] ≈ real_row[c] atol = tolerance
                @test M_j_imag[(r-1)*2^i + c] ≈ imag_row[c] atol = tolerance
            end
        end
        QuEST.destroyComplexMatrixN(M)
   

    for i=1:10
        M_j = rand(Base.Complex{qreal}, 2^i, 2^i)
        M = QuEST.make_QuEST_matrix(M_j)
        real_ptrs = unsafe_wrap(Vector{Ptr{qreal}}, M.real, 2^i)
        imag_ptrs = unsafe_wrap(Vector{Ptr{qreal}}, M.imag, 2^i)
        for r=1:2^i
            real_row = unsafe_wrap(Vector{qreal}, real_ptrs[r], 2^i)
            imag_row = unsafe_wrap(Vector{qreal}, imag_ptrs[r], 2^i)
            for c=1:2^i
                @test real(M_j[r, c]) ≈ real_row[c] atol = tolerance
                @test imag(M_j[r, c]) ≈ imag_row[c] atol = tolerance
            end
        end
        QuEST.destroyComplexMatrixN(M)
    end

    for i=1:10
        m_func(r, c) = 0.12121r +0.343434c*im
        M = QuEST.createComplexMatrixN(i)
        QuEST.fill_ComplexMatrix!(M, m_func)
        real_ptrs = unsafe_wrap(Vector{Ptr{qreal}}, M.real, 2^i)
        imag_ptrs = unsafe_wrap(Vector{Ptr{qreal}}, M.imag, 2^i)
        for r=1:2^i
            real_row = unsafe_wrap(Vector{qreal}, real_ptrs[r], 2^i)
            imag_row = unsafe_wrap(Vector{qreal}, imag_ptrs[r], 2^i)
            for c=1:2^i
                @test real(m_func(r, c)) ≈ real_row[c] atol = tolerance
                @test imag(m_func(r, c)) ≈ imag_row[c] atol = tolerance
            end
        end
        QuEST.destroyComplexMatrixN(M)
    end
end



function test_DiagonalOp()
    env= QuEST.createQuESTEnv()
    for i=1:10
        num_qubits = rand(1:12)
        op = QuEST.createDiagonalOp(num_qubits, env)
        @test op.numQubits == num_qubits
        op_j = rand(Complex{qreal}, 2^num_qubits)
        QuEST.initDiagonalOp(op, real(op_j), imag(op_j))
        reals = unsafe_wrap(Vector{qreal}, op.real, 2^num_qubits)
        imags = unsafe_wrap(Vector{qreal}, op.imag, 2^num_qubits)
        for e = 1:2^num_qubits
            @test reals[e] ≈ real(op_j[e]) atol = tolerance
            @test imags[e] ≈ imag(op_j[e]) atol = tolerance
        end

        
        ind1 = rand(1:2^num_qubits)
        ind2 = rand(1:2^num_qubits)
        start_ind = min(ind1, ind2)
        len_arr = abs(ind1-ind2)+1
        op_j = rand(Complex{qreal}, len_arr)
        test_re = real(op_j)
        test_im = imag(op_j)
        QuEST.setDiagonalOpElems(op, start_ind-1, test_re, test_im, len_arr)
        for ind = start_ind:start_ind+len_arr-1
            @test reals[ind] ≈ test_re[ind-start_ind+1] atol = tolerance
            @test imags[ind] ≈ test_im[ind-start_ind+1] atol = tolerance
        end
        QuEST.syncDiagonalOp(op)
        QuEST.destroyDiagonalOp(op, env)
    end
    QuEST.destroyQuESTEnv(env)
end



function test_controlledMultiQubitUnitary()
    env= QuEST.createQuESTEnv()
    for t=1:10
        numQubits = rand(2:12)
        num_targs = rand(1:numQubits-1)

        targs = Cint.([x for x in 1:num_targs])
        targs_rev = Cint.([x for x in num_targs:-1:1])
        control = 0

        
        M_j = rand(Haar(2), 2^num_targs)
        if qreal == Float32
            M_j = Matrix{Complex{qreal}}(M_j)
        end
        M = QuEST.make_QuEST_matrix(M_j)
    
        qureg = QuEST.createQureg(numQubits, env)
        
        QuEST.pauliX(qureg, 0)
        
        QuEST.controlledMultiQubitUnitary(qureg, control, targs, M)
        
        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)

        X=Matrix([0 1.0;1 0])
        res = numQubits-num_targs-1
        Idm = 1.0*Matrix(I,2^res, 2^res)
        U=kron(Idm, M_j, X)
        v = zeros(2^numQubits)
        v[1]=1.0
        v=U*v
        #print(reals)
        #print(imags)
        for ind = 1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end
        QuEST.destroyComplexMatrixN(M)
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
    
end    

function compactUnitary(qureg       ::QuEST_Types.Qureg,
    targetQubit ::Integer,
    α           ::Complex{Qreal},
    β           ::Complex{Qreal})   ::Nothing

alpha = QuEST_Types.Complex(real(α),imag(α))
beta  = QuEST_Types.Complex(real(β),imag(β))

ccall(:compactUnitary, Cvoid,
(QuEST_Types.Qureg, Cint,        QuEST_Types.Complex, QuEST_Types.Complex),
qureg,              targetQubit, alpha,               beta)
end

function controlledCompactUnitary(qureg          ::QuEST_Types.Qureg,
              controlQubit   ::Integer,
              targetQubit    ::Integer,
              α              ::Complex{Qreal},
              β              ::Complex{Qreal} ) ::Nothing

alpha = QuEST_Types.Complex(real(α),imag(α))
beta  = QuEST_Types.Complex(real(β),imag(β))

ccall(:controlledCompactUnitary,
Cvoid,
(QuEST_Types.Qureg, Cint, Cint, QuEST_Types.Complex, QuEST_Types.Complex),
qureg,
Cint(controlQubit),
Cint(targetQubit),
alpha,
beta)
return nothing
end

function controlledMultiQubitUnitary(qureg   ::QuEST_Types.Qureg,
                 ctrl    ::Integer,
                 targs   ::Vector{QubitIdx},
                 u       ::QuEST_Types.ComplexMatrixN) ::Nothing
@assert 0 ≤ ctrl < qureg.numQubitsRepresented
@assert all( tidx -> 0 ≤ targs[tidx] < qureg.numQubitsRepresented,
1:length(targs) )
@assert ctrl ∉ targs

@assert u.numQubits < qureg.numQubitsRepresented

ccall(:controlledMultiQubitUnitary, Cvoid,
(QuEST_Types.Qureg, Cint, Ptr{Cint}, Cint,          QuEST_Types.ComplexMatrixN),
qureg,              ctrl, targs,     length(targs), u)
end

function _quest_mtx_2(U ::Matrix{Complex{Qreal}}) ::QuEST_Types.ComplexMatrix2
    @assert size(U) == (2,2)
    u = QuEST_Types.ComplexMatrix2(
        ( (real(U[1,1]), real(U[1,2])), (real(U[2,1]), real(U[2,2])) ),
        ( (imag(U[1,1]), imag(U[1,2])), (imag(U[2,1]), imag(U[2,2])) )  )
    return u
end

function _quest_mtx_4(U ::Matrix{Complex{Qreal}}) ::QuEST_Types.ComplexMatrix4
    @assert size(U) == (4,4)
    u = QuEST_Types.ComplexMatrix4( ( ( real(U[1,1]), real(U[1,2]), real(U[1,3]), real(U[1,4]) ),
                          ( real(U[2,1]), real(U[2,2]), real(U[2,3]), real(U[2,4]) ),
                          ( real(U[3,1]), real(U[3,2]), real(U[3,3]), real(U[3,4]) ),
                          ( real(U[4,1]), real(U[4,2]), real(U[4,3]), real(U[4,4]) ) ),
                        ( ( imag(U[1,1]), imag(U[1,2]), imag(U[1,3]), imag(U[1,4]) ),
                          ( imag(U[2,1]), imag(U[2,2]), imag(U[2,3]), imag(U[2,4]) ),
                          ( imag(U[3,1]), imag(U[3,2]), imag(U[3,3]), imag(U[3,4]) ),
                          ( imag(U[4,1]), imag(U[4,2]), imag(U[4,3]), imag(U[4,4]) ) )  )
end

function controlledTwoQubitUnitary(qureg         ::QuEST_Types.Qureg,
                                   controlQubit  ::Integer,
                                   targetQubit1  ::Integer,
                                   targetQubit2  ::Integer,
                                   U             ::Matrix{Complex{Qreal}}) ::Nothing

    @assert 0 ≤ controlQubit < qureg.numQubitsRepresented
    @assert 0 ≤ targetQubit1 < qureg.numQubitsRepresented
    @assert 0 ≤ targetQubit2 < qureg.numQubitsRepresented
    @assert targetQubit1 != targetQubit2
    @assert controlQubit ∉ [targetQubit1,targetQubit2]

    @assert size(U) == (4, 4)
    u = _quest_mtx_4(U)
    ccall(:controlledTwoQubitUnitary,
          Cvoid,
          (QuEST_Types.Qureg, Cint, Cint, Cint, QuEST_Types.ComplexMatrix4),
          qureg,
          Cint(controlQubit),
          Cint(targetQubit1),
          Cint(targetQubit2),
          u)
return nothing
end


function controlledUnitary(qureg        ::QuEST_Types.Qureg,
                           controlQubit ::Integer,
                           targetQubit  ::Integer,
                           U            ::Matrix{Complex{Qreal}}) ::Nothing

    @assert 0 ≤ controlQubit < qureg.numQubitsRepresented
    @assert 0 ≤ targetQubit  < qureg.numQubitsRepresented
    @assert controlQubit != targetQubit

    u = _quest_mtx_2(U)

    ccall(:controlledUnitary,  Cvoid,
          (QuEST_Types.Qureg, Cint,          Cint,          QuEST_Types.ComplexMatrix2),
          qureg,              controlQubit,  targetQubit,   u)
end

# From TartuQC 
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




function getEnvironmentString(env ::QuEST_Types.QuESTEnv, qureg ::QuEST_Types.Qureg) ::String
    cstr = Vector{Cchar}(undef,232)
    ccall(:getEnvironmentString, Cvoid, (QuEST_Types.QuESTEnv, QuEST_Types.Qureg, Ptr{Cchar}), env, qureg, cstr)

    return unsafe_string(pointer(cstr))
end


function seedQuEST(seedarray ::Base.Vector{Culong}) ::Nothing

    @assert  ! isempty(seedarray)

    ccall(:seedQuEST,
          Cvoid, (Ptr{Culong}, Cint),
          seedarray,
          Cint(length(seedarray)))
    return nothing
end         




function mixMultiQubitKrausMap(qureg        :: QuEST_Types.Qureg,
    targetQubits :: Vector{QubitIdx},
    ops          :: Vector{QuEST_Types.ComplexMatrixN}) ::Nothing

local numQubits ::Cint = length(targetQubits)

@assert all( opidx -> ops[opidx].numQubits==numQubits ,
1:length(ops) )

ccall(:mixMultiQubitKrausMap, Cvoid,
(QuEST_Types.Qureg, Ptr{Cint},    Cint,      Ptr{QuEST_Types.ComplexMatrixN}, Cint),
qureg,              targetQubits, numQubits, ops,                             length(ops))

return nothing
end

function measureWithStats(qureg             ::QuEST_Types.Qureg,
    measureQubit      ::Integer)            :: Tuple{Int,Qreal}

outcomeProb = Ref{Qreal}(-1)

outcome = ccall(:measureWithStats, Cint,
(QuEST_Types.Qureg, Cint,            Ref{Qreal}),
qureg,              measureQubit,    outcomeProb)

return outcome, outcomeProb[]

end

function applyMatrixN(qureg             :: QuEST_Types.Qureg,
    targetQubits      :: Vector{QubitIdx},
    u                 :: QuEST_Types.ComplexMatrixN)   ::Nothing

local numTargs  ::Cint = length(targetQubits)
@assert numTargs == u.numQubits

ccall(:applyMatrixN, Cvoid,
(QuEST_Types.Qureg, Ptr{QubitIdx}, Cint,      QuEST_Types.ComplexMatrixN),
qureg,              targetQubits,  numTargs,  u)
end

function applyMultiControlledMatrixN(qureg          :: QuEST_Types.Qureg,
                   ctrlQubits     :: Vector{QubitIdx},
                   targQubits     :: Vector{QubitIdx},
                   u              :: QuEST_Types.ComplexMatrixN)  ::Nothing

local numCtrls ::Cint = length(ctrlQubits)
local numTargs ::Cint = length(targQubits)

@assert numTargs == u.numQubits
@assert isdisjoint(ctrlQubits,targQubits)

ccall(:applyMultiControlledMatrixN, Cvoid,
(QuEST_Types.Qureg, Ptr{Cint},  Cint,      Ptr{Cint},  Cint,     QuEST_Types.ComplexMatrixN),
qureg,              ctrlQubits, numCtrls,  targQubits, numTargs, u)
end

function writeRecordedQASMToFile(qureg    ::QuEST_Types.Qureg,
    filename ::String) ::Nothing
@assert begin
ios = open(filename, "w")
close(ios) === nothing
end

ccall(:writeRecordedQASMToFile, Cvoid, (QuEST_Types.Qureg, Cstring), qureg, filename)
return nothing
end

function initStateFromAmps(qureg     ::QuEST_Types.Qureg,
    amps_real ::Vector{Qreal},
    amps_imag ::Vector{Qreal}      ) ::Nothing

@assert length(amps_real) == qureg.numAmpsTotal
@assert length(amps_imag) == qureg.numAmpsTotal

ccall(:initStateFromAmps,
Cvoid,
(QuEST_Types.Qureg, Ptr{Qreal}, Ptr{Qreal}),
qureg,              amps_real,  amps_imag)
end

function setAmps(qureg        ::QuEST_Types.Qureg,
    startIdx     ::Integer,
    amps_real    ::Vector{Qreal},
    amps_imag    ::Vector{Qreal},
    numAmps      ::Clonglong)   :: Nothing

@assert numAmps             == length(amps_real)
@assert numAmps             == length(amps_imag)
@assert startIdx + numAmps  <= qureg.numAmpsTotal

ccall(:setAmps, Cvoid,
(QuEST_Types.Qureg, Clonglong, Ptr{Qreal}, Ptr{Qreal}, Clonglong),
qureg,              startIdx,  amps_real,  amps_imag,  numAmps)
end


const __ℂ     = QuEST_Types.Complex
const __Qureg = QuEST_Types.Qureg

function setWeightedQureg(factor1    ::Complex{Qreal},
                          qureg1     ::QuEST_Types.Qureg,
                          factor2    ::Complex{Qreal},
                          qureg2     ::QuEST_Types.Qureg,
                          factorOut  ::Complex{Qreal},
                          quregOut   ::QuEST_Types.Qureg)  ::Nothing

    @assert Bool(qureg1.isDensityMatrix) == Bool(qureg2.isDensityMatrix) == Bool(quregOut.isDensityMatrix)
    @assert qureg1.numQubitsRepresented  == qureg2.numQubitsRepresented  == quregOut.numQubitsRepresented

    local factor1_quest   = __ℂ(real(factor1  ),imag(factor1  ))
    local factor2_quest   = __ℂ(real(factor2  ),imag(factor2  ))
    local factorOut_quest = __ℂ(real(factorOut),imag(factorOut))

    ccall(:setWeightedQureg, Cvoid,
          (__ℂ,          __Qureg,  __ℂ,           __Qureg, __ℂ,             __Qureg),
          factor1_quest, qureg1,   factor2_quest, qureg2,  factorOut_quest, quregOut)
end

function controlledRotateAroundAxis(qureg        ::QuEST_Types.Qureg,
    controlQubit ::Integer,
    targetQubit  ::Integer,
    angle        ::Qreal,
    axis         ::Union{NTuple{3,Qreal},Vector{Qreal}})  ::Nothing

@assert 0 ≤ controlQubit < qureg.numQubitsRepresented
@assert 0 ≤ targetQubit  < qureg.numQubitsRepresented
@assert controlQubit != targetQubit

@assert length(axis) == 3

q_axis = QuEST_Types.Vector(axis[1], axis[2], axis[3])

ccall(:controlledRotateAroundAxis,  Cvoid,
(QuEST_Types.Qureg, Cint,           Cint,         Qreal, QuEST_Types.Vector),
qureg,              controlQubit,   targetQubit,  angle, q_axis)
end

