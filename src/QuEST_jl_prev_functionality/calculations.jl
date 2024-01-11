# QuEST_jl/src/base/calculations.jl
#

function calcDensityInnerProduct(ϱ1 ::QuEST_Types.Qureg, ϱ2 ::QuEST_Types.Qureg) ::Qreal

    w = ccall(:calcDensityInnerProduct,
              Qreal,
              (QuEST_Types.Qureg, QuEST_Types.Qureg),
              ϱ1,
              ϱ2)

    return Qreal(w)
end

function calcExpecDiagonalOp(qureg      ::QuEST_Types.Qureg, op         ::QuEST_Types.DiagonalOp)  ::QuEST_Types.Complex

    ret= ccall(:calcExpecDiagonalOp,
               QuEST_Types.Complex,
               (QuEST_Types.Qureg, QuEST_Types.DiagonalOp),
               qureg,
               op)

    return ret
end

function calcExpecPauliHamil(qureg           ::QuEST_Types.Qureg,
                             hamil           ::QuEST_Types.PauliHamil,
                             workspace       ::QuEST_Types.Qureg)       Qreal

    ret = ccall(:calcExpecPauliHamil, Qreal, (QuEST_Types.Qureg, QuEST_Types.PauliHamil, QuEST_Types.Qureg), qureg, hamil, workspace)

    return ret

end

function calcExpecPauliProd(qureg        ::QuEST_Types.Qureg,
                            targetQubits ::Vector{QubitIdx},
                            pauliCodes   ::Vector{QuEST_Types.pauliOpType},
                            workspace    ::QuEST_Types.Qureg)          ::Qreal

    @assert length(targetQubits) == length(pauliCodes)
    #@assert all( σ -> 0 ≤ σ ≤ 3,   pauliCodes )

    expval = ccall(:calcExpecPauliProd,
                   Qreal,
                   (QuEST_Types.Qureg, Ptr{QubitIdx}, Ptr{QuEST_Types.pauliOpType}, Cint,                 QuEST_Types.Qureg),
                   qureg,              targetQubits,  pauliCodes,                   length(targetQubits), workspace)

    return expval
end

function calcExpecPauliSum(qureg         ::QuEST_Types.Qureg,
                           allPauliCodes ::Vector{QuEST_Types.pauliOpType},
                           termCoeffs    ::Vector{Qreal},
                           workspace     ::QuEST_Types.Qureg)          ::Float64

    @assert length(allPauliCodes) ==  length(termCoeffs) * getNumQubits(qureg)
    #@assert all( σ -> 0 ≤ σ ≤ 3,  allPauliCodes )

    ex = ccall(:calcExpecPauliSum,
               Qreal,
               (QuEST_Types.Qureg, Ptr{QuEST_Types.pauliOpType}, Ptr{Qreal}, Cint, QuEST_Types.Qureg),
               qureg,
               allPauliCodes,
               termCoeffs,
               Cint(length(termCoeffs)),
               workspace)

    return ex

end

function calcFidelity(qureg ::QuEST_Types.Qureg, pureState ::QuEST_Types.Qureg) ::Qreal

    fi = ccall(:calcFidelity, Qreal, (QuEST_Types.Qureg, QuEST_Types.Qureg), qureg, pureState)
    return fi

end

function calcHilbertSchmidtDistance(a ::QuEST_Types.Qureg, b ::QuEST_Types.Qureg) ::Qreal

    hsd = ccall(:calcHilbertSchmidtDistance, Qreal, (QuEST_Types.Qureg, QuEST_Types.Qureg), a, b)
    return hsd

end

function calcInnerProduct(bra ::QuEST_Types.Qureg, ket ::QuEST_Types.Qureg) ::Complex{Qreal}

    w = ccall(:calcInnerProduct, QuEST_Types.Complex, (QuEST_Types.Qureg, QuEST_Types.Qureg), bra, ket)
    return Complex{Qreal}(w.real,w.imag)

end

function calcProbOfOutcome(qureg        ::QuEST_Types.Qureg,
                           measureQubit ::Integer,
                           outcome      ::Integer)   ::Qreal

    p = ccall(:calcProbOfOutcome,
              Qreal,
              (QuEST_Types.Qureg, QubitIdx,      Cint),
              qureg,              measureQubit,  outcome)
    return p
end

function calcPurity(qureg ::QuEST_Types.Qureg) ::Qreal

    pu = ccall(:calcPurity, Qreal, (QuEST_Types.Qureg,), qureg)
    return pu

end

function calcTotalProb(qureg ::QuEST_Types.Qureg) ::Qreal

    one = ccall(:calcTotalProb, Qreal, (QuEST_Types.Qureg,), qureg)
    return one

end

function getAmp(qureg ::QuEST_Types.Qureg,  idx ::Integer) ::Complex{Qreal}

    α = ccall(:getAmp, QuEST_Types.Complex,
              (QuEST_Types.Qureg, Clonglong),
              qureg,              idx)
    return Complex{Qreal}(α.real,α.imag)

end

function getDensityAmp(qureg ::QuEST_Types.Qureg, row ::Integer, col ::Integer) ::Complex{Qreal}

    α = ccall(:getDensityAmp,
              QuEST_Types.Complex,
              (QuEST_Types.Qureg, Clonglong, Clonglong),
              qureg,
              Clonglong(row),
              Clonglong(col))

    return Complex{Qreal}(α.real, α.imag)

end

function getImagAmp(qureg      ::QuEST_Types.Qureg,
                    index      ::Integer)      ::Qreal

    ret = ccall(:getImagAmp, Qreal, (QuEST_Types.Qureg, Clonglong), qureg, Clonglong(index))
    return ret

end

function getNumAmps(qureg ::QuEST_Types.Qureg) ::Clonglong

    return ccall(:getNumAmps, Clonglong, (QuEST_Types.Qureg,), qureg)

end

function getNumQubits(qureg ::QuEST_Types.Qureg) ::QubitIdx

    return ccall(:getNumQubits, Cint, (QuEST_Types.Qureg,), qureg)

end

function getProbAmp(qureg ::QuEST_Types.Qureg,
                    idx   ::Integer) :: Qreal

    p = ccall(:getProbAmp, Qreal, (QuEST_Types.Qureg, Clonglong), qureg, Clonglong(idx))
    return p

end

function getRealAmp(qureg ::QuEST_Types.Qureg,
                    idx   ::Integer) :: Qreal

    p = ccall(:getRealAmp, Qreal, (QuEST_Types.Qureg, Clonglong), qureg, Clonglong(idx))
    return p

end

#EOF
