# QuEST_jl/src/base/unitaries.jl
#

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


function controlledNot(qureg         ::QuEST_Types.Qureg,
                       controlQubit  ::Integer,
                       targetQubit   ::Integer)         ::Nothing

    @assert 0 ≤ controlQubit < qureg.numQubitsRepresented
    @assert 0 ≤ targetQubit  < qureg.numQubitsRepresented
    @assert controlQubit != targetQubit

    ccall(:controlledNot, Cvoid,
          (QuEST_Types.Qureg, QubitIdx,      QubitIdx),
          qureg,              controlQubit,  targetQubit)
end


function controlledPauliY(qureg         ::QuEST_Types.Qureg,
                          controlQubit  ::Integer,
                          targetQubit   ::Integer)         ::Nothing

    @assert 0 ≤ controlQubit < qureg.numQubitsRepresented
    @assert 0 ≤ targetQubit  < qureg.numQubitsRepresented
    @assert controlQubit != targetQubit

    ccall(:controlledPauliY,  Cvoid,
          (QuEST_Types.Qureg, Cint,          Cint),
          qureg,              controlQubit,  targetQubit)
end

function controlledPhaseFlip(qureg         ::QuEST_Types.Qureg,
                             controlQubit  ::Integer,
                             targetQubit   ::Integer)   ::Nothing

    @assert 0 ≤ controlQubit < qureg.numQubitsRepresented
    @assert 0 ≤ targetQubit  < qureg.numQubitsRepresented
    @assert controlQubit != targetQubit

    ccall(:controlledPhaseFlip, Cvoid,
          (QuEST_Types.Qureg, Cint,         Cint),
          qureg,              controlQubit, targetQubit)
end


function controlledPhaseShift(qureg    ::QuEST_Types.Qureg,
                              idQubit1 ::Integer,
                              idQubit2 ::Integer,
                              angle    ::Qreal)   ::Nothing

    @assert 0 ≤ idQubit1 < qureg.numQubitsRepresented
    @assert 0 ≤ idQubit2  < qureg.numQubitsRepresented
    @assert idQubit1 != idQubit2

    ccall(:controlledPhaseShift,  Cvoid,
          (QuEST_Types.Qureg, Cint,     Cint,      Qreal),
          qureg,              idQubit1, idQubit2,  angle)
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

function controlledRotateX(qureg         ::QuEST_Types.Qureg,
                           controlQubit  ::Integer,
                           targetQubit   ::Integer,
                           angle         ::Qreal)  ::Nothing

    @assert 0 ≤ controlQubit < qureg.numQubitsRepresented
    @assert 0 ≤ targetQubit  < qureg.numQubitsRepresented
    @assert controlQubit != targetQubit

    ccall(:controlledRotateX, Cvoid,
          (QuEST_Types.Qureg, Cint,          Cint,          Qreal),
          qureg,              controlQubit,  targetQubit,   angle)
end


function controlledRotateY(qureg         ::QuEST_Types.Qureg,
                           controlQubit  ::Integer,
                           targetQubit   ::Integer,
                           angle         ::Qreal)  ::Nothing

    @assert 0 ≤ controlQubit < qureg.numQubitsRepresented
    @assert 0 ≤ targetQubit  < qureg.numQubitsRepresented
    @assert controlQubit != targetQubit

    ccall(:controlledRotateY,    Cvoid,
          (QuEST_Types.Qureg, Cint,          Cint,        Qreal),
          qureg,              controlQubit,  targetQubit, angle)
end

function controlledRotateZ(qureg         ::QuEST_Types.Qureg,
    controlQubit  ::Integer,
    targetQubit   ::Integer,
    angle         ::Qreal)  ::Nothing

    ccall(:controlledRotateZ,
          Cvoid,
          (QuEST_Types.Qureg, Cint, Cint, Qreal),
          qureg,
          Cint(controlQubit),
          Cint(targetQubit),
          angle)
return nothing
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

function hadamard(qureg       ::QuEST_Types.Qureg,
                  targetQubit ::Integer) ::Nothing

    ccall(:hadamard, Cvoid,
          (QuEST_Types.Qureg, Cint),
          qureg,              targetQubit)
end

function multiControlledMultiQubitUnitary(qureg  ::QuEST_Types.Qureg,
                                          ctrls  ::Vector{QubitIdx},
                                          targs  ::Vector{QubitIdx},
                                          u      ::QuEST_Types.ComplexMatrixN) ::Nothing

    ccall(:multiControlledMultiQubitUnitary, Cvoid,
          (QuEST_Types.Qureg, Ptr{Cint}, Cint,          Ptr{Cint}, Cint,          QuEST_Types.ComplexMatrixN),
          qureg,              ctrls,     length(ctrls), targs,     length(targs), u)
end

function multiControlledPhaseFlip(qureg         ::QuEST_Types.Qureg,
                                  controlQubits ::Vector{QubitIdx}) ::Nothing

    ccall(:multiControlledPhaseFlip, Cvoid,
          (QuEST_Types.Qureg,     Ptr{Cint},     Cint),
          qureg,                  controlQubits, length(controlQubits))
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
