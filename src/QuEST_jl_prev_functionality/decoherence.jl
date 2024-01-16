 

function mixKrausMap(qureg          :: QuEST_Types.Qureg,
                     targetQubit    :: Integer,
                     ops            :: Vector{QuEST_Types.ComplexMatrix2}) ::Nothing

    ccall(:mixKrausMap, Cvoid,
          (QuEST_Types.Qureg, Cint,         Ptr{QuEST_Types.ComplexMatrix2}, Cint),
          qureg,              targetQubit,  ops,                             length(ops))

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

function mixPauli(qureg         :: QuEST_Types.Qureg,
                  targetQubit   :: Integer,
                  probX         :: Qreal,
                  probY         :: Qreal,
                  probZ         :: Qreal)   ::Nothing

    ccall(:mixPauli, Cvoid,
          (QuEST_Types.Qureg, Cint,          Qreal, Qreal, Qreal),
          qureg,              targetQubit,   probX, probY, probZ)

    return nothing
end

function mixTwoQubitDephasing(qureg         :: QuEST_Types.Qureg,
                              qubit1        :: Integer,
                              qubit2        :: Integer,
                              prob          :: Qreal)     ::Nothing

    ccall(:mixTwoQubitDephasing, Cvoid,
          (QuEST_Types.Qureg, Cint,   Cint,   Qreal),
          qureg,              qubit1, qubit2, prob)

    return nothing
end

function mixTwoQubitDepolarising(qureg         :: QuEST_Types.Qureg,
                                 qubit1        :: Integer,
                                 qubit2        :: Integer,
                                 prob          :: Qreal)     ::Nothing

    ccall(:mixTwoQubitDepolarising, Cvoid,
    (QuEST_Types.Qureg, Cint,   Cint,    Qreal),
    qureg,              qubit1, qubit2,  prob)

    return nothing
end

function mixTwoQubitKrausMap(qureg         :: QuEST_Types.Qureg,
                             qubit1        :: Integer,
                             qubit2        :: Integer,
                             ops           :: Vector{QuEST_Types.ComplexMatrix4}) ::Nothing

    ccall(:mixTwoQubitKrausMap, Cvoid,
          (QuEST_Types.Qureg, Cint,   Cint,   Ptr{QuEST_Types.ComplexMatrix4}, Cint),
          qureg,              qubit1, qubit2, ops,                             length(ops))

    return nothing
end

#EOF
