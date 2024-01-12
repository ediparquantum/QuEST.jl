

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









function multiControlledPhaseShift(qureg         ::QuEST_Types.Qureg,
                                   controlQubits ::Vector{QubitIdx},
                                   angle         ::Qreal)         ::Nothing

    ccall(:multiControlledPhaseShift, Cvoid,
          (QuEST_Types.Qureg, Ptr{Cint},     Cint,                  Qreal),
          qureg,              controlQubits, length(controlQubits), angle)
end

function multiControlledTwoQubitUnitary(qureg           ::QuEST_Types.Qureg,
                                        controlQubits   ::Vector{QubitIdx},
                                        targetQubit1    ::Integer,
                                        targetQubit2    ::Integer,
                                        U               ::Matrix{Complex{Qreal}}) ::Nothing

    @assert 0 ≤ length(controlQubits) ≤ qureg.numQubitsRepresented -2
    @assert 0 ≤ targetQubit1 < qureg.numQubitsRepresented
    @assert 0 ≤ targetQubit2 < qureg.numQubitsRepresented
    @assert targetQubit1 != targetQubit2
    @assert targetQubit1 ∉ controlQubits
    @assert targetQubit2 ∉ controlQubits

    u = _quest_mtx_4(U)

    ccall(:multiControlledTwoQubitUnitary, Cvoid,
          (QuEST_Types.Qureg, Ptr{Cint},     Cint,                  Cint,          Cint,           QuEST_Types.ComplexMatrix4),
          qureg,              controlQubits, length(controlQubits), targetQubit1,  targetQubit2,   u)
end

function multiControlledUnitary(qureg         ::QuEST_Types.Qureg,
                                controlQubits ::Vector{QubitIdx},
                                targetQubit   ::Integer,
                                U             ::Matrix{Complex{Qreal}}) ::Nothing

    u = _quest_mtx_2(U)

    ccall(:multiControlledUnitary,  Cvoid,
          (QuEST_Types.Qureg, Ptr{Cint},     Cint,                  Cint,          QuEST_Types.ComplexMatrix2),
          qureg,              controlQubits, length(controlQubits), targetQubit,   u)
end

function multiQubitUnitary(qureg ::QuEST_Types.Qureg,
                           targs ::Vector{QubitIdx},
                           u     ::QuEST_Types.ComplexMatrixN)  ::Nothing

    @assert u.numQubits ≤ qureg.numQubitsRepresented

    ccall(:multiQubitUnitary, Cvoid,
          (QuEST_Types.Qureg, Ptr{Cint}, Cint,           QuEST_Types.ComplexMatrixN),
          qureg,              targs,     length(targs),  u)
end

function multiRotatePauli(qureg         ::QuEST_Types.Qureg,
                          targetQubits  ::Vector{QubitIdx},
                          targetPaulis  ::Vector{QuEST_Types.pauliOpType},
                          angle         ::Qreal)          ::Nothing

    @assert length(targetQubits) == length(targetPaulis)
    #@assert all( σ -> 0 ≤ σ ≤ 3,   targetPaulis )

    ccall(:multiRotatePauli, Cvoid,
          (QuEST_Types.Qureg, Ptr{Cint},     Ptr{QuEST_Types.pauliOpType}, Cint,                 Qreal),
          qureg,              targetQubits,  targetPaulis,                 length(targetPaulis), angle)
end

function multiRotateZ(qureg      ::QuEST_Types.Qureg,
                      qubits     ::Vector{QubitIdx},
                      angle      ::Qreal)         ::Nothing
    ccall(:multiRotateZ, Cvoid,
          (QuEST_Types.Qureg, Ptr{Cint}, Cint,           Qreal),
          qureg,              qubits,    length(qubits), angle)
end

function multiStateControlledUnitary(qureg           ::QuEST_Types.Qureg,
                                     controlQubits   ::Vector{QubitIdx},
                                     controlState    ::Vector{Cint},
                                     targetQubit     ::Integer,
                                     U               ::Matrix{Complex{Qreal}}) ::Nothing

    @assert length(controlQubits) == length(controlState)

    u = _quest_mtx_2(U)

    ccall(:multiStateControlledUnitary, Cvoid,
          (QuEST_Types.Qureg, Ptr{Cint},      Ptr{Cint},    Cint,                  Cint,        QuEST_Types.ComplexMatrix2),
          qureg,              controlQubits,  controlState, length(controlQubits), targetQubit, u)
end

function pauliX(qureg       ::QuEST_Types.Qureg,
                targetQubit ::Integer)         :: Nothing

    ccall(:pauliX, Cvoid,
          (QuEST_Types.Qureg, Cint),
          qureg,              targetQubit)
end

function pauliY(qureg       ::QuEST_Types.Qureg,
                targetQubit ::Integer)         :: Nothing

    ccall(:pauliY, Cvoid,
          (QuEST_Types.Qureg, Cint),
          qureg,              targetQubit)
end
function pauliZ(qureg ::QuEST_Types.Qureg, targetQubit ::Integer):: Nothing

    ccall(:pauliZ,  Cvoid,
          (QuEST_Types.Qureg, Cint),
          qureg,              targetQubit)
end

function phaseShift(qureg       ::QuEST_Types.Qureg,
                    targetQubit ::Integer,
                    angle       ::Qreal) ::Nothing

    ccall(:phaseShift, Cvoid,
          (QuEST_Types.Qureg, Cint,        Qreal),
          qureg,              targetQubit, angle)
end

function rotateAroundAxis(qureg         ::QuEST_Types.Qureg,
                          rotQubit      ::Integer,
                          angle         ::Qreal,
                          axis          ::Union{NTuple{3,Qreal},Vector{Qreal}}) ::Nothing

    @assert length(axis) == 3

    q_axis = QuEST_Types.Vector(axis[1], axis[2], axis[3])

    ccall(:rotateAroundAxis,  Cvoid,
          (QuEST_Types.Qureg, Cint,     Qreal, QuEST_Types.Vector),
          qureg,              rotQubit, angle, q_axis)
end

function rotateX(qureg    ::QuEST_Types.Qureg,
                 rotQubit ::Integer,
                 angle    ::Qreal)           ::Nothing

    ccall(:rotateX,  Cvoid,
          (QuEST_Types.Qureg, Cint,      Qreal),
          qureg,              rotQubit,  angle)
end

function rotateY(qureg    ::QuEST_Types.Qureg,
                 rotQubit ::Integer,
                 angle    ::Qreal)           ::Nothing

    ccall(:rotateY, Cvoid,
          (QuEST_Types.Qureg, Cint,      Qreal),
          qureg,              rotQubit,  angle)
end

function rotateZ(qureg    ::QuEST_Types.Qureg,
                 rotQubit ::Integer,
                 angle    ::Qreal)           ::Nothing

    ccall(:rotateZ, Cvoid,
          (QuEST_Types.Qureg, Cint,     Qreal),
          qureg,              rotQubit, angle)
end

function sGate(qureg       ::QuEST_Types.Qureg,
               targetQubit ::Integer)         ::Nothing

    ccall(:sGate, Cvoid,
          (QuEST_Types.Qureg, Cint),
          qureg,              targetQubit)
end

function sqrtSwapGate(qureg  ::QuEST_Types.Qureg,
                      qubit1 ::Integer,
                      qubit2 ::Integer)         ::Nothing

    ccall(:sqrtSwapGate, Cvoid,
          (QuEST_Types.Qureg, Cint,   Cint),
          qureg,              qubit1, qubit2)
end

function swapGate(qureg  ::QuEST_Types.Qureg,
                  qubit1 ::Integer,
                  qubit2 ::Integer)         ::Nothing

    ccall(:swapGate,  Cvoid,
          (QuEST_Types.Qureg, Cint,   Cint),
          qureg,              qubit1, qubit2)
end

function tGate(qureg       ::QuEST_Types.Qureg,
               targetQubit ::Integer)         ::Nothing

    ccall(:tGate, Cvoid,
          (QuEST_Types.Qureg, Cint),
          qureg,              targetQubit)
end

function twoQubitUnitary(qureg           ::QuEST_Types.Qureg,
                         targetQubit1    ::Integer,
                         targetQubit2    ::Integer,
                         U               ::Matrix{Complex{Qreal}}) ::Nothing

    u = _quest_mtx_4(U)

    ccall(:twoQubitUnitary, Cvoid,
          (QuEST_Types.Qureg, Cint,          Cint,          QuEST_Types.ComplexMatrix4),
          qureg,              targetQubit1,  targetQubit2,  u)
end

function unitary(qureg           ::QuEST_Types.Qureg,
                 targetQubit     ::Integer,
                 U               ::Matrix{Complex{Qreal}}) ::Nothing

    u = _quest_mtx_2(U)

    ccall(:unitary, Cvoid,
          (QuEST_Types.Qureg, Cint,           QuEST_Types.ComplexMatrix2),
          qureg,              targetQubit,    u)
end

#EOF
