# QuEST_jl/src/base/operators.jl
#

function applyDiagonalOp(qureg          :: QuEST_Types.Qureg,
                         op             :: QuEST_Types.DiagonalOp)  ::Nothing

    ccall(:applyDiagonalOp, Cvoid, (QuEST_Types.Qureg, QuEST_Types.DiagonalOp), qureg, op)
    nothing;
end

function applyMatrix2(qureg             :: QuEST_Types.Qureg,
                      targetQubit       :: Integer,
                      u                 :: QuEST_Types.ComplexMatrix2)  ::Nothing

    ccall(:applyMatrix2, Cvoid,
          (QuEST_Types.Qureg, Cint,        QuEST_Types.ComplexMatrix2),
          qureg,              targetQubit, u)
end

function applyMatrix4(qureg             :: QuEST_Types.Qureg,
                      targetQubit1      :: Integer,
                      targetQubit2      :: Integer,
                      u                 :: QuEST_Types.ComplexMatrix4)   ::Nothing

    ccall(:applyMatrix4,  Cvoid,
          (QuEST_Types.Qureg, Cint,          Cint,         QuEST_Types.ComplexMatrix4),
          qureg,              targetQubit1,  targetQubit2, u)
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

function applyPauliHamil(inQureg        :: QuEST_Types.Qureg,
                         hamil          :: QuEST_Types.PauliHamil,
                         outQureg       :: QuEST_Types.Qureg)        ::Nothing

    ccall(:applyPauliHamil, Cvoid,
          (QuEST_Types.Qureg, QuEST_Types.PauliHamil, QuEST_Types.Qureg),
          inQureg,            hamil,                  outQureg)
end

function applyPauliSum(inQureg        ::QuEST_Types.Qureg,
                       allPauliCodes  ::Vector{QuEST_Types.pauliOpType},
                       termCoeffs     ::Vector{Qreal},
                       numSumTerms    ::Integer,
                       outQureg       ::QuEST_Types.Qureg)          ::Nothing

    @assert length(allPauliCodes) == numSumTerms * inQureg.numQubitsRepresented

    ccall(:applyPauliSum, Cvoid,
          (QuEST_Types.Qureg, Ptr{QuEST_Types.pauliOpType}, Ptr{Qreal}, Cint,        QuEST_Types.Qureg),
          inQureg,            allPauliCodes,                termCoeffs, numSumTerms, outQureg)
end

function applyTrotterCircuit(qureg          :: QuEST_Types.Qureg,
                             hamil          :: QuEST_Types.PauliHamil,
                             time           :: Qreal,
                             order          :: Integer,
                             reps           :: Integer)          ::Nothing

    ccall(:applyTrotterCircuit, Cvoid,
          (QuEST_Types.Qureg, QuEST_Types.PauliHamil, Qreal, Cint,  Cint),
          qureg,              hamil,                  time,  order, reps)
end

#EOF
