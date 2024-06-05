struct ComplexMatrixN
    numQubits::Cint
    real::Ptr{Ptr{Cdouble}}
    imag::Ptr{Ptr{Cdouble}}
end

function bindArraysToStackComplexMatrixN(numQubits, re, im, reStorage, imStorage)
    @ccall libquest.bindArraysToStackComplexMatrixN(numQubits::Cint, re::Ptr{Ptr{Cdouble}}, im::Ptr{Ptr{Cdouble}}, reStorage::Ptr{Ptr{Cdouble}}, imStorage::Ptr{Ptr{Cdouble}})::ComplexMatrixN
end

@enum phaseGateType::UInt32 begin
    SIGMA_Z = 0
    S_GATE = 1
    T_GATE = 2
end

struct QASMLogger
    buffer::Cstring
    bufferSize::Cint
    bufferFill::Cint
    isLogging::Cint
end

struct ComplexArray
    real::Ptr{Cdouble}
    imag::Ptr{Cdouble}
end




@enum pauliOpType::UInt32 begin
    PAULI_I = 1
    PAULI_X = 2
    PAULI_Y = 3
    PAULI_Z = 4
end

@enum cPauliOpType::UInt32 begin
    cPAULI_I = 0
    cPAULI_X = 1
    cPAULI_Y = 2
    cPAULI_Z = 3
end

function convert_to_cPauliOpType(pauliOp::pauliOpType)
    pauli_int = Int(pauliOp)
    if pauli_int < 1
        error("Minimum pauliaOpType is 1")
    elseif pauli_int > 4
        error("Maximum pauliOpType is 4")
    else
        return cPauliOpType(pauli_int-1)
    end
end

function convert_to_cPauliOpType(pauliInt::Int)
    pauli_int = pauliInt
    if pauli_int < 1
        error("Minimum pauliaOpType is 1")
    elseif pauli_int > 4
        error("Maximum pauliOpType is 4")
    else
        return cPauliOpType(pauli_int-1)
    end
end



struct QComplex
    real::Cdouble
    imag::Cdouble
end


struct ComplexMatrix2
    real::NTuple{2, NTuple{2, Cdouble}}
    imag::NTuple{2, NTuple{2, Cdouble}}
end

struct ComplexMatrix4
    real::NTuple{4, NTuple{4, Cdouble}}
    imag::NTuple{4, NTuple{4, Cdouble}}
end

struct QVector
    x::Cdouble
    y::Cdouble
    z::Cdouble
end

@enum phaseFunc::UInt32 begin
    NORM = 0
    SCALED_NORM = 1
    INVERSE_NORM = 2
    SCALED_INVERSE_NORM = 3
    SCALED_INVERSE_SHIFTED_NORM = 4
    PRODUCT = 5
    SCALED_PRODUCT = 6
    INVERSE_PRODUCT = 7
    SCALED_INVERSE_PRODUCT = 8
    DISTANCE = 9
    SCALED_DISTANCE = 10
    INVERSE_DISTANCE = 11
    SCALED_INVERSE_DISTANCE = 12
    SCALED_INVERSE_SHIFTED_DISTANCE = 13
    SCALED_INVERSE_SHIFTED_WEIGHTED_DISTANCE = 14
end

@enum bitEncoding::UInt32 begin
    UNSIGNED = 0
    TWOS_COMPLEMENT = 1
end

struct PauliHamil
    pauliCodes::Ptr{pauliOpType}
    termCoeffs::Ptr{Cdouble}
    numSumTerms::Cint
    numQubits::Cint
end

struct DiagonalOp
    numQubits::Cint
    numElemsPerChunk::Clonglong
    numChunks::Cint
    chunkId::Cint
    real::Ptr{Cdouble}
    imag::Ptr{Cdouble}
    deviceOperator::ComplexArray
end

struct SubDiagonalOp
    numQubits::Cint
    numElems::Clonglong
    real::Ptr{Cdouble}
    imag::Ptr{Cdouble}
end

struct Qureg
    isDensityMatrix::Cint
    numQubitsRepresented::Cint
    numQubitsInStateVec::Cint
    numAmpsPerChunk::Clonglong
    numAmpsTotal::Clonglong
    chunkId::Cint
    numChunks::Cint
    stateVec::ComplexArray
    pairStateVec::ComplexArray
    deviceStateVec::ComplexArray
    firstLevelReduction::Ptr{Cdouble}
    secondLevelReduction::Ptr{Cdouble}
    cuStateVec::Ptr{Cvoid}
    deviceCuStateVec::Ptr{Cvoid}
    cuConfig::Ptr{Ptr{Cvoid}}
    qasmLog::Ptr{QASMLogger}
end

struct QuESTEnv
    rank::Cint
    numRanks::Cint
    seeds::Ptr{Culong}
    numSeeds::Cint
    cuConfig::Ptr{Ptr{Cvoid}}
end






function createQureg(numQubits, env)
    @ccall libquest.createQureg(numQubits::Cint, env::QuESTEnv)::Qureg
end

function createDensityQureg(numQubits, env)
    @ccall libquest.createDensityQureg(numQubits::Cint, env::QuESTEnv)::Qureg
end

function createCloneQureg(qureg, env)
    @ccall libquest.createCloneQureg(qureg::Qureg, env::QuESTEnv)::Qureg
end

function destroyQureg(qureg, env)
    @ccall libquest.destroyQureg(qureg::Qureg, env::QuESTEnv)::Cvoid
end

function createComplexMatrixN(numQubits)
    @ccall libquest.createComplexMatrixN(numQubits::Cint)::ComplexMatrixN
end

function destroyComplexMatrixN(matr)
    @ccall libquest.destroyComplexMatrixN(matr::ComplexMatrixN)::Cvoid
end

function initComplexMatrixN(m, real, imag)
    @ccall libquest.initComplexMatrixN(m::ComplexMatrixN, real::Ptr{Ptr{Cdouble}}, imag::Ptr{Ptr{Cdouble}})::Cvoid
end

function createPauliHamil(numQubits, numSumTerms)
    @ccall libquest.createPauliHamil(numQubits::Cint, numSumTerms::Cint)::PauliHamil
end

function destroyPauliHamil(hamil)
    @ccall libquest.destroyPauliHamil(hamil::PauliHamil)::Cvoid
end

function createPauliHamilFromFile(fn)
    @ccall libquest.createPauliHamilFromFile(fn::Cstring)::PauliHamil
end

function initPauliHamil(hamil, coeffs, codes)
    @ccall libquest.initPauliHamil(hamil::PauliHamil, coeffs::Ptr{Cdouble}, codes::Ptr{pauliOpType})::Cvoid
end

function createDiagonalOp(numQubits, env)
    @ccall libquest.createDiagonalOp(numQubits::Cint, env::QuESTEnv)::DiagonalOp
end

function destroyDiagonalOp(op, env)
    @ccall libquest.destroyDiagonalOp(op::DiagonalOp, env::QuESTEnv)::Cvoid
end

function syncDiagonalOp(op)
    @ccall libquest.syncDiagonalOp(op::DiagonalOp)::Cvoid
end

function initDiagonalOp(op, real, imag)
    @ccall libquest.initDiagonalOp(op::DiagonalOp, real::Ptr{Cdouble}, imag::Ptr{Cdouble})::Cvoid
end

function initDiagonalOpFromPauliHamil(op, hamil)
    @ccall libquest.initDiagonalOpFromPauliHamil(op::DiagonalOp, hamil::PauliHamil)::Cvoid
end

function createDiagonalOpFromPauliHamilFile(fn, env)
    @ccall libquest.createDiagonalOpFromPauliHamilFile(fn::Cstring, env::QuESTEnv)::DiagonalOp
end

function setDiagonalOpElems(op, startInd, real, imag, numElems)
    test_qubit_present(qureg,startInd)
    startInd = c_shift_index(startInd)
    @ccall libquest.setDiagonalOpElems(op::DiagonalOp, startInd::Clonglong, real::Ptr{Cdouble}, imag::Ptr{Cdouble}, numElems::Clonglong)::Cvoid
end

function applyDiagonalOp(qureg, op)
    @ccall libquest.applyDiagonalOp(qureg::Qureg, op::DiagonalOp)::Cvoid
end

function calcExpecDiagonalOp(qureg, op)
    @ccall libquest.calcExpecDiagonalOp(qureg::Qureg, op::DiagonalOp)::QComplex
end

function createSubDiagonalOp(numQubits)
    @ccall libquest.createSubDiagonalOp(numQubits::Cint)::SubDiagonalOp
end

function destroySubDiagonalOp(op)
    @ccall libquest.destroySubDiagonalOp(op::SubDiagonalOp)::Cvoid
end

function diagonalUnitary(qureg, targets,op)
    test_qubit_present(qureg,targets)
    targets = [c_shift_index(tq) for tq in targets]
    numTargets = length(targets)
    @ccall libquest.diagonalUnitary(qureg::Qureg, targets::Ptr{Cint}, numTargets::Cint, op::SubDiagonalOp)::Cvoid
end

function controlledMultiQubitUnitary(qureg, ctrl, targets, u)
    test_qubit_present(qureg,ctrl)
    test_qubit_present(qureg,targets)
    ctrl = c_shift_index(ctrl)
    targets = [c_shift_index(tq) for tq in targets]
    numTargets = length(targets)
    @ccall libquest.controlledMultiQubitUnitary(qureg::Qureg, ctrl::Cint, targets::Ptr{Cint}, numTargets::Cint, u::ComplexMatrixN)::Cvoid
end

function multiControlledMultiQubitUnitary(qureg, ctrls, targets, u)
    test_qubit_present(qureg,ctrl)
    test_qubit_present(qureg,targets)
    ctrls = Cint.([c_shift_index(cq) for cq in ctrls])
    targets = Cint.([c_shift_index(tq) for tq in targets])
    numCtrls = Cint(length(ctrls))
    numTargets = Cint(length(targets))
    @ccall libquest.multiControlledMultiQubitUnitary(qureg::Qureg, pointer(ctrls)::Ptr{Cint}, numCtrls::Cint, pointer(targets)::Ptr{Cint}, numTargets::Cint, make_QuEST_matrix(u)::ComplexMatrixN)::Cvoid
end

function applyGateSubDiagonalOp(qureg, targets, op)
    test_qubit_present(qureg,targets)
    targets = [c_shift_index(t) for t in targets]
    numTargets = length(targets)
    @ccall libquest.applyGateSubDiagonalOp(qureg::Qureg, targets::Ptr{Cint}, numTargets::Cint, op::SubDiagonalOp)::Cvoid
end

function applySubDiagonalOp(qureg, targets, op)
    test_qubit_present(qureg,targets)
    targets = [c_shift_index(t) for t in targets]
    numTargets = length(targets)
    @ccall libquest.applySubDiagonalOp(qureg::Qureg, targets::Ptr{Cint}, numTargets::Cint, op::SubDiagonalOp)::Cvoid
end

function reportState(qureg)
    @ccall libquest.reportState(qureg::Qureg)::Cvoid
end

function reportStateToScreen(qureg, env, reportRank)
    @ccall libquest.reportStateToScreen(qureg::Qureg, env::QuESTEnv, reportRank::Cint)::Cvoid
end

function reportQuregParams(qureg)
    @ccall libquest.reportQuregParams(qureg::Qureg)::Cvoid
end

function reportPauliHamil(hamil)
    @ccall libquest.reportPauliHamil(hamil::PauliHamil)::Cvoid
end

function getNumQubits(qureg)
    @ccall libquest.getNumQubits(qureg::Qureg)::Cint
end

function getNumAmps(qureg)
    @ccall libquest.getNumAmps(qureg::Qureg)::Clonglong
end

function initBlankState(qureg)
    @ccall libquest.initBlankState(qureg::Qureg)::Cvoid
end

function initZeroState(qureg)
    @ccall libquest.initZeroState(qureg::Qureg)::Cvoid
end

function initPlusState(qureg)
    @ccall libquest.initPlusState(qureg::Qureg)::Cvoid
end

function initClassicalState(qureg, stateInd)
    test_qubit_present(qureg, stateInd)
    stateInd = c_shift_index(stateInd)
    @ccall libquest.initClassicalState(qureg::Qureg, stateInd::Clonglong)::Cvoid
end

function initPureState(qureg, pure)
    @ccall libquest.initPureState(qureg::Qureg, pure::Qureg)::Cvoid
end

function initDebugState(qureg)
    @ccall libquest.initDebugState(qureg::Qureg)::Cvoid
end

function initStateFromAmps(qureg, reals, imags)
    @ccall libquest.initStateFromAmps(qureg::Qureg, reals::Ptr{Cdouble}, imags::Ptr{Cdouble})::Cvoid
end

function setAmps(qureg, startInd, reals, imags, numAmps)
    test_qubit_present(qureg, startInd)
    startInd = c_shift_index(startInd)
    @ccall libquest.setAmps(qureg::Qureg, startInd::Clonglong, reals::Ptr{Cdouble}, imags::Ptr{Cdouble}, numAmps::Clonglong)::Cvoid
end

function setDensityAmps(qureg, startRow, startCol, reals, imags, numAmps)
    test_row_col_in_size(qureg,startRow)
    test_row_col_in_size(qureg, startCol)
    startRow = c_shift_index(startRow)
    startCol = c_shift_index(startCol)
    @ccall libquest.setDensityAmps(qureg::Qureg, startRow::Clonglong, startCol::Clonglong, reals::Ptr{Cdouble}, imags::Ptr{Cdouble}, numAmps::Clonglong)::Cvoid
end

function setQuregToPauliHamil(qureg, hamil)
    @ccall libquest.setQuregToPauliHamil(qureg::Qureg, hamil::PauliHamil)::Cvoid
end

function cloneQureg(targetQureg, copyQureg)
    @ccall libquest.cloneQureg(targetQureg::Qureg, copyQureg::Qureg)::Cvoid
end

function phaseShift(qureg, targetQubit, angle)
    test_qubit_present(qureg, targetQubit)
    targetQubit = c_shift_index(targetQubit)
    @ccall libquest.phaseShift(qureg::Qureg, targetQubit::Cint, angle::Cdouble)::Cvoid
end

function controlledPhaseShift(qureg, qubit1, qubit2, angle)
    test_qubit_present(qureg, qubit1)
    test_qubit_present(qureg, qubit2)
    qubit1 = c_shift_index(qubit1)
    qubit2 = c_shift_index(qubit2)
    @ccall libquest.controlledPhaseShift(qureg::Qureg, qubit1::Cint, qubit2::Cint, angle::Cdouble)::Cvoid
end

function multiControlledPhaseShift(qureg, controlQubits, angle)
    test_qubit_present(qureg, controlQubits)
    controlQubits = Cint.([c_shift_index(cq) for cq in controlQubits])
    @ccall libquest.multiControlledPhaseShift(qureg::Qureg, pointer(controlQubits)::Ptr{Cint}, length(controlQubits)::Cint, angle::Cdouble)::Cvoid
end




function controlledPhaseFlip(qureg, qubit1, qubit2)
    test_qubit_present(qureg, qubit1)
    test_qubit_present(qureg, qubit2)
    qubit1 = c_shift_index(qubit1)
    qubit2 = c_shift_index(qubit2)
    @ccall libquest.controlledPhaseFlip(qureg::Qureg, qubit1::Cint, qubit2::Cint)::Cvoid
end

function multiControlledPhaseFlip(qureg, controlQubits) 
    test_qubit_present(qureg, controlQubits)
    controlQubits = Cint.([c_shift_index(cq) for cq in controlQubits])
    @ccall libquest.multiControlledPhaseFlip(qureg::Qureg, pointer(controlQubits)::Ptr{Cint}, length(controlQubits)::Cint)::Cvoid
end




function sGate(qureg, targetQubit)
    test_qubit_present(qureg, targetQubit)
    targetQubit = c_shift_index(targetQubit)
    @ccall libquest.sGate(qureg::Qureg, targetQubit::Cint)::Cvoid
end

function tGate(qureg, targetQubit)
    test_qubit_present(qureg, targetQubit)
    targetQubit = c_shift_index(targetQubit)
    @ccall libquest.tGate(qureg::Qureg, targetQubit::Cint)::Cvoid
end

function createQuESTEnv()
    @ccall libquest.createQuESTEnv()::QuESTEnv
end

function destroyQuESTEnv(env)
    @ccall libquest.destroyQuESTEnv(env::QuESTEnv)::Cvoid
end

function syncQuESTEnv(env)
    @ccall libquest.syncQuESTEnv(env::QuESTEnv)::Cvoid
end

function syncQuESTSuccess(successCode)
    @ccall libquest.syncQuESTSuccess(successCode::Cint)::Cint
end

function reportQuESTEnv(env)
    @ccall libquest.reportQuESTEnv(env::QuESTEnv)::Cvoid
end

function getEnvironmentString(env, str)
    @ccall libquest.getEnvironmentString(env::QuESTEnv, str::Ptr{Cchar})::Cvoid
end

function copyStateToGPU(qureg)
    @ccall libquest.copyStateToGPU(qureg::Qureg)::Cvoid
end

function copyStateFromGPU(qureg)
    @ccall libquest.copyStateFromGPU(qureg::Qureg)::Cvoid
end

function copySubstateToGPU(qureg, startInd, numAmps)
    test_qubit_present(qureg, startInd)
    startInd = c_shift_index(startInd)
    @ccall libquest.copySubstateToGPU(qureg::Qureg, startInd::Clonglong, numAmps::Clonglong)::Cvoid
end

function copySubstateFromGPU(qureg, startInd, numAmps)
    test_qubit_present(qureg, startInd)
    startInd = c_shift_index(startInd)
    @ccall libquest.copySubstateFromGPU(qureg::Qureg, startInd::Clonglong, numAmps::Clonglong)::Cvoid
end

function getAmp(qureg, index)
    test_qubit_present(qureg, index)
    index = c_shift_index(index)
    @ccall libquest.getAmp(qureg::Qureg, index::Clonglong)::QComplex
end

function getRealAmp(qureg, index)
    test_qubit_present(qureg, index)
    index = c_shift_index(index)
    @ccall libquest.getRealAmp(qureg::Qureg, index::Clonglong)::Cdouble
end

function getImagAmp(qureg, index)
    test_qubit_present(qureg, index)
    index = c_shift_index(index)
    @ccall libquest.getImagAmp(qureg::Qureg, index::Clonglong)::Cdouble
end

function getProbAmp(qureg, index)
    test_qubit_present(qureg, index)
    index = c_shift_index(index)
    @ccall libquest.getProbAmp(qureg::Qureg, index::Clonglong)::Cdouble
end

function getDensityAmp(qureg, row, col)
    test_row_col_in_size(qureg,row)
    test_row_col_in_size(qureg, col)
    row = c_shift_index(row)
    col = c_shift_index(col)
    @ccall libquest.getDensityAmp(qureg::Qureg, row::Clonglong, col::Clonglong)::QComplex
end

function calcTotalProb(qureg)
    @ccall libquest.calcTotalProb(qureg::Qureg)::Cdouble
end

function compactUnitary(qureg, targetQubit, alpha, beta)
    test_qubit_present(qureg,targetQubit)
    targetQubit = c_shift_index(targetQubit)
    @ccall libquest.compactUnitary(qureg::Qureg, targetQubit::Cint, alpha::QComplex, beta::QComplex)::Cvoid
end

function unitary(qureg, targetQubit, u)
    test_qubit_present(qureg,targetQubit)
    targetQubit = c_shift_index(targetQubit)
    @ccall libquest.unitary(qureg::Qureg, targetQubit::Cint, u::ComplexMatrix2)::Cvoid
end

function rotateX(qureg, rotQubit, angle)
    test_qubit_present(qureg,rotQubit)
    rotQubit = c_shift_index(rotQubit)
    @ccall libquest.rotateX(qureg::Qureg, rotQubit::Cint, angle::Cdouble)::Cvoid
end

function rotateY(qureg, rotQubit, angle)
    test_qubit_present(qureg,rotQubit)
    rotQubit = c_shift_index(rotQubit)
    @ccall libquest.rotateY(qureg::Qureg, rotQubit::Cint, angle::Cdouble)::Cvoid
end

function rotateZ(qureg, rotQubit, angle)
    test_qubit_present(qureg,rotQubit)
    rotQubit = c_shift_index(rotQubit)
    @ccall libquest.rotateZ(qureg::Qureg, rotQubit::Cint, angle::Cdouble)::Cvoid
end

function rotateAroundAxis(qureg, rotQubit, angle, axis)
    test_qubit_present(qureg,rotQubit)
    rotQubit = c_shift_index(rotQubit)
    @ccall libquest.rotateAroundAxis(qureg::Qureg, rotQubit::Cint, angle::Cdouble, axis::QVector)::Cvoid
end

function controlledRotateX(qureg, controlQubit, targetQubit, angle)
    test_qubit_present(qureg,controlQubit)
    test_qubit_present(qureg,targetQubit)
    controlQubit = c_shift_index(controlQubit)
    targetQubit = c_shift_index(targetQubit)
    @ccall libquest.controlledRotateX(qureg::Qureg, controlQubit::Cint, targetQubit::Cint, angle::Cdouble)::Cvoid
end

function controlledRotateY(qureg, controlQubit, targetQubit, angle)
    test_qubit_present(qureg,controlQubit)
    test_qubit_present(qureg,targetQubit)
    controlQubit = c_shift_index(controlQubit)
    targetQubit = c_shift_index(targetQubit)
    @ccall libquest.controlledRotateY(qureg::Qureg, controlQubit::Cint, targetQubit::Cint, angle::Cdouble)::Cvoid
end

function controlledRotateZ(qureg, controlQubit, targetQubit, angle)
    test_qubit_present(qureg,controlQubit)
    test_qubit_present(qureg,targetQubit)
    controlQubit = c_shift_index(controlQubit)
    targetQubit = c_shift_index(targetQubit)
    @ccall libquest.controlledRotateZ(qureg::Qureg, controlQubit::Cint, targetQubit::Cint, angle::Cdouble)::Cvoid
end

function controlledRotateAroundAxis(qureg, controlQubit, targetQubit, angle, axis)
    test_qubit_present(qureg,controlQubit)
    test_qubit_present(qureg,targetQubit)
    controlQubit = c_shift_index(controlQubit)
    targetQubit = c_shift_index(targetQubit)
    @ccall libquest.controlledRotateAroundAxis(qureg::Qureg, controlQubit::Cint, targetQubit::Cint, angle::Cdouble, axis::QVector)::Cvoid
end

function controlledCompactUnitary(qureg, controlQubit, targetQubit, alpha, beta)
    test_qubit_present(qureg,controlQubit)
    test_qubit_present(qureg,targetQubit)
    controlQubit = c_shift_index(controlQubit)
    targetQubit = c_shift_index(targetQubit)
    @ccall libquest.controlledCompactUnitary(qureg::Qureg, controlQubit::Cint, targetQubit::Cint, alpha::QComplex, beta::QComplex)::Cvoid
end

function controlledUnitary(qureg, controlQubit, targetQubit, u)
    test_qubit_present(qureg,controlQubit)
    test_qubit_present(qureg,targetQubit)
    test_qubits_different(controlQubit,targetQubit)
    controlQubit = c_shift_index(controlQubit)
    targetQubit = c_shift_index(targetQubit)    
    u = make_QuEST_matrix_2x2(u)
    @ccall libquest.controlledUnitary(qureg::Qureg, controlQubit::Cint, targetQubit::Cint, u::ComplexMatrix2)::Cvoid
end

function multiControlledUnitary(qureg, controlQubits,targetQubit, u)
    test_qubit_present(qureg,controlQubits)
    test_qubit_present(qureg,targetQubit)
    controlQubits = [c_shift_index(cq) for cq in controlQubits]
    targetQubit = c_shift_index(targetQubit)
    numControlQubits = length(numControlQubits)
    @ccall libquest.multiControlledUnitary(qureg::Qureg, controlQubits::Ptr{Cint}, numControlQubits::Cint, targetQubit::Cint, u::ComplexMatrix2)::Cvoid
end

function pauliX(qureg, targetQubit)
    test_qubit_present(qureg,targetQubit)
    targetQubit = c_shift_index(targetQubit)
    @ccall libquest.pauliX(qureg::Qureg, targetQubit::Cint)::Cvoid
end

function pauliY(qureg, targetQubit)
    test_qubit_present(qureg,targetQubit)
    targetQubit = c_shift_index(targetQubit)
    @ccall libquest.pauliY(qureg::Qureg, targetQubit::Cint)::Cvoid
end

function pauliZ(qureg, targetQubit)
    test_qubit_present(qureg,targetQubit)
    targetQubit = c_shift_index(targetQubit)
    @ccall libquest.pauliZ(qureg::Qureg, targetQubit::Cint)::Cvoid
end

function hadamard(qureg, targetQubit)
    test_qubit_present(qureg,targetQubit)
    targetQubit = c_shift_index(targetQubit)
    @ccall libquest.hadamard(qureg::Qureg, targetQubit::Cint)::Cvoid
end

function controlledNot(qureg, controlQubit, targetQubit)
    test_qubit_present(qureg,targetQubit)
    test_qubit_present(qureg,controlQubit)
    controlQubit = c_shift_index(controlQubit)
    targetQubit = c_shift_index(targetQubit)
    @ccall libquest.controlledNot(qureg::Qureg, controlQubit::Cint, targetQubit::Cint)::Cvoid
end

function multiControlledMultiQubitNot(qureg, ctrls, targets)
    test_qubit_present(qureg,ctrls)
    test_qubit_present(qureg,targets)
    ctrls = [c_shift_index(cq) for cq in ctrls] 
    targets = [c_shift_index(tq) for tq in targets]
    numCtrls = length(numCtrls)
    numTargets = length(numTargets)
    @ccall libquest.multiControlledMultiQubitNot(qureg::Qureg, ctrls::Ptr{Cint}, numCtrls::Cint, targets::Ptr{Cint}, numTargets::Cint)::Cvoid
end

function multiQubitNot(qureg, targets)
    test_qubit_present(qureg,targets)
    targets = [c_shift_index(tq) for tq in targets]
    numTargets = length(numTargets)
    @ccall libquest.multiQubitNot(qureg::Qureg, targets::Ptr{Cint}, numTargets::Cint)::Cvoid
end

function controlledPauliY(qureg, controlQubit, targetQubit)
    test_qubit_present(qureg,controlQubit)
    test_qubit_present(qureg,targetQubit)
    controlQubit = c_shift_index(controlQubit)
    targetQubit = c_shift_index(targetQubit)
    @ccall libquest.controlledPauliY(qureg::Qureg, controlQubit::Cint, targetQubit::Cint)::Cvoid
end

function calcProbOfOutcome(qureg, measureQubit, outcome)
    test_qubit_present(qureg,measureQubit)
    measureQubit = c_shift_index(measureQubit)
    @ccall libquest.calcProbOfOutcome(qureg::Qureg, measureQubit::Cint, outcome::Cint)::Cdouble
end

function calcProbOfAllOutcomes(outcomeProbs, qureg, qubits)
    test_qubit_present(qureg,qubits)
    qubits = [c_shift_index(q) for q in qubits]
    numQubits = length(numQubits)
    @ccall libquest.calcProbOfAllOutcomes(outcomeProbs::Ptr{Cdouble}, qureg::Qureg, qubits::Ptr{Cint}, numQubits::Cint)::Cvoid
end

function collapseToOutcome(qureg, measureQubit, outcome)
    test_qubit_present(qureg,measureQubit)
    measureQubit = c_shift_index(measureQubit)
    @ccall libquest.collapseToOutcome(qureg::Qureg, measureQubit::Cint, outcome::Cint)::Cdouble
end

function measure(qureg, measureQubit)
    test_qubit_present(qureg,measureQubit)
    measureQubit = c_shift_index(measureQubit)
    @ccall libquest.measure(qureg::Qureg, measureQubit::Cint)::Cint
end

function measureWithStats(qureg, measureQubit, outcomeProb)
    test_qubit_present(qureg,measureQubit)
    measureQubit = c_shift_index(measureQubit)
    @ccall libquest.measureWithStats(qureg::Qureg, measureQubit::Cint, outcomeProb::Ptr{Cdouble})::Cint
end

function calcInnerProduct(bra, ket)
    @ccall libquest.calcInnerProduct(bra::Qureg, ket::Qureg)::QComplex
end

function calcDensityInnerProduct(rho1, rho2)
    test_if_qureg_is_density_matrix(rho1)
    test_if_qureg_is_density_matrix(rho2)
    test_if_qureg_is_density_matrix(rho1,rho2)
    @ccall libquest.calcDensityInnerProduct(rho1::Qureg, rho2::Qureg)::Cdouble
end

function calcPurity(qureg)
    @ccall libquest.calcPurity(qureg::Qureg)::Cdouble
end

function calcFidelity(qureg, pureState)
    test_if_at_least_one_qureg_is_state_vector(pureState)
    @ccall libquest.calcFidelity(qureg::Qureg, pureState::Qureg)::Cdouble
end

function calcExpecPauliProd(qureg, targetQubits, pauliCodes, workspace)
    test_qubit_present(qureg,targetQubits)
    test_length_pauli_prod(targetQubits,pauliCodes)
    targetQubits = pointer(Cint.([c_shift_index(tq) for tq in targetQubits]))
    pauliCodes = [get_pauli_code(str) for str in pauliCodes] |> pointer
    numTargets = length(targetQubits)
    @ccall libquest.calcExpecPauliProd(qureg::Qureg, targetQubits::Ptr{Cint}, pauliCodes::Ptr{pauliOpType}, numTargets::Cint, workspace::Qureg)::Cdouble
end

function calcExpecPauliSum(qureg, allPauliCodes, termCoeffs, workspace)
    test_length_pauli_sum(qureg,numSumTerms,allPauliCodes)
    numSumTerms = length(termCoeffs)
    allPauliCodes = [get_pauli_code(str) for str in allPauliCodes] |> pointer
    termCoeffs = pointer(termCoeffs)
    @ccall libquest.calcExpecPauliSum(qureg::Qureg, allPauliCodes::Ptr{pauliOpType}, termCoeffs::Ptr{Cdouble}, numSumTerms::Cint, workspace::Qureg)::Cdouble
end

function calcExpecPauliHamil(qureg, hamil, workspace)
    @ccall libquest.calcExpecPauliHamil(qureg::Qureg, hamil::PauliHamil, workspace::Qureg)::Cdouble
end

function calcHilbertSchmidtDistance(a, b)
    @ccall libquest.calcHilbertSchmidtDistance(a::Qureg, b::Qureg)::Cdouble
end

function seedQuESTDefault(env)
    @ccall libquest.seedQuESTDefault(env::Ptr{QuESTEnv})::Cvoid
end

function seedQuEST(env, seedArray, numSeeds)
    @ccall libquest.seedQuEST(env::Ptr{QuESTEnv}, seedArray::Ptr{Culong}, numSeeds::Cint)::Cvoid
end

function getQuESTSeeds(env, seeds, numSeeds)
    @ccall libquest.getQuESTSeeds(env::QuESTEnv, seeds::Ptr{Ptr{Culong}}, numSeeds::Ptr{Cint})::Cvoid
end

function startRecordingQASM(qureg)
    @ccall libquest.startRecordingQASM(qureg::Qureg)::Cvoid
end

function stopRecordingQASM(qureg)
    @ccall libquest.stopRecordingQASM(qureg::Qureg)::Cvoid
end

function clearRecordedQASM(qureg)
    @ccall libquest.clearRecordedQASM(qureg::Qureg)::Cvoid
end

function printRecordedQASM(qureg)
    @ccall libquest.printRecordedQASM(qureg::Qureg)::Cvoid
end

function writeRecordedQASMToFile(qureg, filename)
    @ccall libquest.writeRecordedQASMToFile(qureg::Qureg, filename::Cstring)::Cvoid
end

function mixDephasing(qureg, targetQubit, prob)
    test_qubit_present(qureg,targetQubit)
    test_is_density_matrix(qureg)
    test_is_valid_probability(prob,0.5)
    targetQubit = c_shift_index(targetQubit)
    @ccall libquest.mixDephasing(qureg::Qureg, targetQubit::Cint, prob::Cdouble)::Cvoid
end

function mixTwoQubitDephasing(qureg, qubit1, qubit2, prob)
    test_qubit_present(qureg,qubit1)
    test_qubit_present(qureg,qubit2)
    test_is_density_matrix(qureg)
    test_is_valid_probability(prob,0.75)
    qubit1 = c_shift_index(qubit1)
    qubit2 = c_shift_index(qubit2)
    @ccall libquest.mixTwoQubitDephasing(qureg::Qureg, qubit1::Cint, qubit2::Cint, prob::Cdouble)::Cvoid
end

function mixDepolarising(qureg, targetQubit, prob)
    test_qubit_present(qureg,targetQubit)
    test_is_density_matrix(qureg)
    test_is_valid_probability(prob,0.75)
    targetQubit = c_shift_index(targetQubit)
    @ccall libquest.mixDepolarising(qureg::Qureg, targetQubit::Cint, prob::Cdouble)::Cvoid
end

function mixDamping(qureg, targetQubit, prob)
    test_qubit_present(qureg,targetQubit)
    test_is_density_matrix(qureg)
    test_is_probability(prob)
    targetQubit = c_shift_index(targetQubit)
    @ccall libquest.mixDamping(qureg::Qureg, targetQubit::Cint, prob::Cdouble)::Cvoid
end

function mixTwoQubitDepolarising(qureg, qubit1, qubit2, prob)
    test_qubit_present(qureg,qubit1)
    test_qubit_present(qureg,qubit2)
    test_is_density_matrix(qureg)
    test_is_valid_probability(prob, 15/16)
    qubit1 = c_shift_index(qubit1)
    qubit2 = c_shift_index(qubit2)
    @ccall libquest.mixTwoQubitDepolarising(qureg::Qureg, qubit1::Cint, qubit2::Cint, prob::Cdouble)::Cvoid
end

function mixPauli(qureg, targetQubit, prob)
    test_qubit_present(qureg,targetQubit)
    probX, probY, probZ = prob
    test_is_density_matrix(qureg)
    test_is_valid_probability(probX, 1 - probX - probY - probZ)
    test_is_valid_probability(probY, 1 - probX - probY - probZ)
    test_is_valid_probability(probZ, 1 - probX - probY - probZ)
    targetQubit = c_shift_index(targetQubit)
    @ccall libquest.mixPauli(qureg::Qureg, targetQubit::Cint, probX::Cdouble, probY::Cdouble, probZ::Cdouble)::Cvoid
end

function mixDensityMatrix(combineQureg, prob, otherQureg)
    test_is_density_matrix(combineQureg)
    test_is_density_matrix(otherQureg)
    test_is_probability(prob)
    @ccall libquest.mixDensityMatrix(combineQureg::Qureg, prob::Cdouble, otherQureg::Qureg)::Cvoid
end

function mixKrausMap(qureg, target, ops)
    test_qubit_present(qureg,target)
    test_kraus_sum_is_identity(ops)
    test_max_kraus_operators(ops,4)
    test_is_density_matrix(qureg)
    ops = [make_QuEST_matrix_2x2(o) for o in ops] 
    numOps = Cint(length(ops))
    target = c_shift_index(target)
    @ccall libquest.mixKrausMap(qureg::Qureg, target::Cint, ops::Ptr{ComplexMatrix2}, numOps::Cint)::Cvoid
end

function mixTwoQubitKrausMap(qureg, target1, target2, ops)
    test_qubit_present(qureg, target1)
    test_qubit_present(qureg, target2)
    test_kraus_sum_is_identity(ops)
    test_max_kraus_operators(ops,16)
    test_is_density_matrix(qureg)
    ops = [make_QuEST_matrix_4x4(o) for o in ops] 
    numOps = Cint(length(ops))
    target1 = c_shift_index(target1)
    target2 = c_shift_index(target2)
    @ccall libquest.mixTwoQubitKrausMap(qureg::Qureg, target1::Cint, target2::Cint, ops::Ptr{ComplexMatrix4}, numOps::Cint)::Cvoid
end

function mixMultiQubitKrausMap(qureg, targets, ops)
    test_qubit_present(qureg, targets)
    test_kraus_sum_is_identity(ops)
    rows,cols = length(ops) != 1 ? size(ops[1]) : size(ops)
    test_kraus_operator_dimension_square(rows,cols)
    num_qubits_in_op = Int(log2(rows))
    test_max_kraus_operators(ops,(2*num_qubits_in_op)^2)
    test_is_density_matrix(qureg)
    ops = [make_QuEST_matrix_NxN(o) for o in ops]
    targets = [c_shift_index(tq) for tq in targets]
    numTargets = Cint(length(targets))
    numOps = Cint(length(ops))
    @ccall libquest.mixMultiQubitKrausMap(qureg::Qureg, targets::Ptr{Cint}, numTargets::Cint, ops::Ptr{ComplexMatrixN}, numOps::Cint)::Cvoid
end

function mixNonTPKrausMap(qureg, target, ops)
    test_qubit_present(qureg, target)
    test_is_density_matrix(qureg)
    test_max_kraus_operators(ops,4)
    test_is_in_qubit_range(qureg, target)
    ops = [make_QuEST_matrix_2x2(o) for o in ops] 
    target = c_shift_index(target)
    numOps = Cint(length(ops))
    @ccall libquest.mixNonTPKrausMap(qureg::Qureg, target::Cint, ops::Ptr{ComplexMatrix2}, numOps::Cint)::Cvoid
end

function mixNonTPTwoQubitKrausMap(qureg, target1, target2, ops)
    test_qubit_present(qureg, target1)
    test_qubit_present(qureg, target2)
    test_is_density_matrix(qureg)
    test_max_kraus_operators(ops,16)
    ops = [make_QuEST_matrix_4x4(o) for o in ops]
    target1 = c_shift_index(target1)
    target2 = c_shift_index(target2)
    numOps = Cint(length(ops))
    @ccall libquest.mixNonTPTwoQubitKrausMap(qureg::Qureg, target1::Cint, target2::Cint, ops::Ptr{ComplexMatrix4}, numOps::Cint)::Cvoid
end

function mixNonTPMultiQubitKrausMap(qureg, targets, ops)
    test_qubit_present(qureg, targets)
    rows,cols = length(ops) != 1 ? size(ops[1]) : size(ops)
    test_kraus_operator_dimension_square(rows,cols)
    num_qubits_in_op = Int(log2(rows))
    test_max_kraus_operators(ops,(2*num_qubits_in_op)^2)
    test_is_density_matrix(qureg)
    ops = [make_QuEST_matrix_NxN(o) for o in ops]
    targets = [c_shift_index(tq) for tq in targets]
    numTargets = Cint(length(targets))
    numOps = Cint(length(ops))
    @ccall libquest.mixNonTPMultiQubitKrausMap(qureg::Qureg, targets::Ptr{Cint}, numTargets::Cint, ops::Ptr{ComplexMatrixN}, numOps::Cint)::Cvoid
end


function multiStateControlledUnitary(qureg, controlQubits, controlState,targetQubit, u)
    test_qubit_present(qureg, controlQubits)
    test_qubit_present(qureg, targetQubit)
    test_qubit_present(qureg, controlState)
    controlQubits = [c_shift_index(cq) for cq in controlQubits]
    targetQubit = c_shift_index(targetQubit)
    targetQubit = c_shift_index(controlState)
    numControlQubits = length(numControlQubits)
    @ccall libquest.multiStateControlledUnitary(qureg::Qureg, controlQubits::Ptr{Cint}, controlState::Ptr{Cint}, numControlQubits::Cint, targetQubit::Cint, u::ComplexMatrix2)::Cvoid
end

function multiRotateZ(qureg, qubits, angle)
    test_qubit_present(qureg, qubits)
    qubits = [c_shift_index(q) for q in qubits]
    numQubits = length(qubits)
    @ccall libquest.multiRotateZ(qureg::Qureg, qubits::Ptr{Cint}, numQubits::Cint, angle::Cdouble)::Cvoid
end

function multiRotatePauli(qureg::Qureg, targetQubits::Vector{Int}, targetPaulis::Union{Vector{pauliOpType},Vector{Int}}, angle)
    test_qubit_present(qureg, targetQubits)
    targetQubits = [c_shift_index(q) for q in targetQubits] |> pointer
    targetPaulis = [convert_to_cPauliOpType(tps) for tps in targetPaulis] |> pointer
    numTargets = length(targetQubits)
    @ccall libquest.multiRotatePauli(qureg::Qureg, targetQubits::Ptr{Cint}, targetPaulis::Ptr{cPauliOpType}, numTargets::Cint, angle::Cdouble)::Cvoid
end

function multiControlledMultiRotateZ(qureg, controlQubits,targetQubits,angle)
    test_qubit_present(qureg, controlQubits)
    test_qubit_present(qureg, targetQubits)
    controlQubits = [c_shift_index(cq) for cq in controlQubits]
    targetQubits = [c_shift_index(tq) for tq in targetQubits]
    numControls = length(numControls)
    numTargets = length(numTargets)
    @ccall libquest.multiControlledMultiRotateZ(qureg::Qureg, controlQubits::Ptr{Cint}, numControls::Cint, targetQubits::Ptr{Cint}, numTargets::Cint, angle::Cdouble)::Cvoid
end

function multiControlledMultiRotatePauli(qureg, controlQubits, targetQubits, targetPaulis::Union{Vector{pauliOpType},Vector{Int}},angle)
    test_qubit_present(qureg, controlQubits)
    test_qubit_present(qureg, targetQubits)
    controlQubits = [c_shift_index(cq) for cq in controlQubits]
    targetQubits = [c_shift_index(tq) for tq in targetQubits]
    targetPaulis = [convert_to_cPauliOpType(tps) for tps in targetPaulis] |> pointer
    numControls = length(numControls)
    numTargets = length(numTargets)
    @ccall libquest.multiControlledMultiRotatePauli(qureg::Qureg, controlQubits::Ptr{Cint}, numControls::Cint, targetQubits::Ptr{Cint}, targetPaulis::Ptr{pauliOpType}, numTargets::Cint, angle::Cdouble)::Cvoid
end



function twoQubitUnitary(qureg, targetQubit1, targetQubit2, u)
    test_qubit_present(qureg, targetQubit1)
    test_qubit_present(qureg, targetQubit2)
    targetQubit1 = c_shift_index(targetQubit1)
    targetQubit2 = c_shift_index(targetQubit2)
    @ccall libquest.twoQubitUnitary(qureg::Qureg, targetQubit1::Cint, targetQubit2::Cint, u::ComplexMatrix4)::Cvoid
end



function controlledTwoQubitUnitary(qureg, controlQubit, targetQubit1, targetQubit2, u)
    test_qubit_present_multi_control(qureg, controlQubit)
    test_qubit_present(qureg, targetQubit1)
    test_qubit_present(qureg, targetQubit2)
    test_qubits_different(targetQubit1, targetQubit2)
    test_qubits_different(controlQubit, targetQubit1)
    test_qubits_different(controlQubit, targetQubit2)
    controlQubit = c_shift_index(controlQubit)
    targetQubit1 = c_shift_index(targetQubit1)
    targetQubit2 = c_shift_index(targetQubit2)
    u = make_QuEST_matrix_4x4(u)
    @ccall libquest.controlledTwoQubitUnitary(qureg::Qureg, controlQubit::Cint, targetQubit1::Cint, targetQubit2::Cint, u::ComplexMatrix4)::Cvoid
end

function multiControlledTwoQubitUnitary(qureg, controlQubits, targetQubit1, targetQubit2, u)
    test_qubit_present_multi_control(qureg, controlQubits)
    test_qubit_present(qureg, targetQubit1)
    test_qubit_present(qureg, targetQubit2)
    test_qubits_different(targetQubit1, targetQubit2)
    test_qubits_different(controlQubit, targetQubit1)
    test_qubits_different(controlQubit, targetQubit2)
    controlQubits = Cint.([c_shift_index(cq) for cq in controlQubits])
    targetQubit1 = c_shift_index(targetQubit1)
    targetQubit2 = c_shift_index(targetQubit2)
    @ccall libquest.multiControlledTwoQubitUnitary(qureg::Qureg, pointer(controlQubits)::Ptr{Cint}, length(controlQubits)::Cint, targetQubit1::Cint, targetQubit2::Cint, make_QuEST_matrix_4x4(u)::ComplexMatrix4)::Cvoid
end


function multiQubitUnitary(qureg, targets, u)
    test_qubit_present(qureg, targets)
    targets = [c_shift_index(tq) for tq in targets]
    numTargets = lenght(numTargets)
    @ccall libquest.multiQubitUnitary(qureg::Qureg, targets::Ptr{Cint}, numTargets::Cint, u::ComplexMatrixN)::Cvoid
end



function setWeightedQureg(fac1, qureg1, fac2, qureg2, facOut, out)
    @ccall libquest.setWeightedQureg(fac1::QComplex, qureg1::Qureg, fac2::QComplex, qureg2::Qureg, facOut::QComplex, out::Qureg)::Cvoid
end

function swapGate(qureg, qubit1, qubit2)
    test_qubit_present(qureg, qubit1)
    test_qubit_present(qureg, qubit2)
    qubit1 = c_shift_index(qubit1)
    qubit2 = c_shift_index(qubit2)
    @ccall libquest.swapGate(qureg::Qureg, qubit1::Cint, qubit2::Cint)::Cvoid
end

function sqrtSwapGate(qureg, qubit1, qubit2)
    test_qubit_present(qureg, qubit1)
    test_qubit_present(qureg, qubit2)
    qubit1 = c_shift_index(qubit1)
    qubit2 = c_shift_index(qubit2)
    @ccall libquest.sqrtSwapGate(qureg::Qureg, qubit1::Cint, qubit2::Cint)::Cvoid
end

function applyPauliSum(inQureg, allPauliCodes, termCoeffs, numSumTerms, outQureg)
    @ccall libquest.applyPauliSum(inQureg::Qureg, allPauliCodes::Ptr{pauliOpType}, termCoeffs::Ptr{Cdouble}, numSumTerms::Cint, outQureg::Qureg)::Cvoid
end

function applyPauliHamil(inQureg, hamil, outQureg)
    @ccall libquest.applyPauliHamil(inQureg::Qureg, hamil::PauliHamil, outQureg::Qureg)::Cvoid
end

function applyTrotterCircuit(qureg, hamil, time, order, reps)
    @ccall libquest.applyTrotterCircuit(qureg::Qureg, hamil::PauliHamil, time::Cdouble, order::Cint, reps::Cint)::Cvoid
end

function applyMatrix2(qureg, targetQubit, u)
    test_qubit_present(qureg, targetQubit)
    targetQubit = c_shift_index(targetQubit)
    @ccall libquest.applyMatrix2(qureg::Qureg, targetQubit::Cint, u::ComplexMatrix2)::Cvoid
end

function applyMatrix4(qureg, targetQubit1, targetQubit2, u)
    test_qubit_present(qureg, targetQubit1)
    test_qubit_present(qureg, targetQubit2)
    targetQubit1 = c_shift_index(targetQubit1)
    targetQubit2 = c_shift_index(targetQubit2)
    @ccall libquest.applyMatrix4(qureg::Qureg, targetQubit1::Cint, targetQubit2::Cint, u::ComplexMatrix4)::Cvoid
end

function applyMatrixN(qureg, targets, u)
    test_qubit_present(qureg, targets)
    targets = [c_shift_index(tq) for tq in targets]
    numTargets = length(numTargets)
    @ccall libquest.applyMatrixN(qureg::Qureg, targets::Ptr{Cint}, numTargets::Cint, u::ComplexMatrixN)::Cvoid
end

function applyGateMatrixN(qureg, targets, u)
    test_qubit_present(qureg, targets)
    targets = [c_shift_index(tq) for tq in targets]
    numTargets = length(numTargets)
    @ccall libquest.applyGateMatrixN(qureg::Qureg, targets::Ptr{Cint}, numTargets::Cint, u::ComplexMatrixN)::Cvoid
end

function applyMultiControlledGateMatrixN(qureg, ctrls, targets, m)
    test_qubit_present(qureg, ctrls)
    test_qubit_present(qureg, targets)
    test_if_vector_has_repetitions(ctrls)
    test_if_vector_has_repetitions(targets)
    test_if_vec1_is_in_vec2_vise_versa(ctrls,targets)
    test_if_vec_has_length_at_least_one(ctrls)
    test_if_vec_has_length_at_least_one(targets)
    ctrls = [c_shift_index(cq) for cq in ctrls]
    targets = [c_shift_index(tq) for tq in targets]
    numCtrls = length(numCtrls)
    numTargets = length(numTargets)
    @ccall libquest.applyMultiControlledGateMatrixN(qureg::Qureg, ctrls::Ptr{Cint}, numCtrls::Cint, targets::Ptr{Cint}, numTargets::Cint, m::ComplexMatrixN)::Cvoid
end



function applyMultiControlledMatrixN(qureg, ctrls, targets, u)
    test_qubit_present(qureg, ctrls)
    test_qubit_present(qureg, targets)
    ctrls = [c_shift_index(cq) for cq in ctrls]
    targets = [c_shift_index(tq) for tq in targets]
    numCtrls = length(numCtrls)
    numTargets = length(numTargets)
    @ccall libquest.applyMultiControlledMatrixN(qureg::Qureg, ctrls::Ptr{Cint}, numCtrls::Cint, targets::Ptr{Cint}, numTargets::Cint, u::ComplexMatrixN)::Cvoid
end

function invalidQuESTInputError(errMsg, errFunc)
    @ccall libquest.invalidQuESTInputError(errMsg::Cstring, errFunc::Cstring)::Cvoid
end

function applyPhaseFunc(qureg, qubits, encoding, coeffs, exponents, numTerms)
    test_qubit_present(qureg, qubits)
    qubits = [c_shift_index(q) for q in qubits]
    numQubits = length(numQubits)
    @ccall libquest.applyPhaseFunc(qureg::Qureg, qubits::Ptr{Cint}, numQubits::Cint, encoding::bitEncoding, coeffs::Ptr{Cdouble}, exponents::Ptr{Cdouble}, numTerms::Cint)::Cvoid
end

function applyPhaseFuncOverrides(qureg, qubits, encoding, coeffs, exponents, numTerms, overrideInds, overridePhases, numOverrides)
    test_qubit_present(qureg, qubits)
    test_qubit_present(qureg, overrideInds)
    qubits = [c_shift_index(q) for q in qubits]
    overrideInds = [c_shift_index(oi) for oi in overrideInds]
    numQubits = length(numQubits)
    @ccall libquest.applyPhaseFuncOverrides(qureg::Qureg, qubits::Ptr{Cint}, numQubits::Cint, encoding::bitEncoding, coeffs::Ptr{Cdouble}, exponents::Ptr{Cdouble}, numTerms::Cint, overrideInds::Ptr{Clonglong}, overridePhases::Ptr{Cdouble}, numOverrides::Cint)::Cvoid
end

function applyMultiVarPhaseFunc(qureg, qubits, numQubitsPerReg, numRegs, encoding, coeffs, exponents, numTermsPerReg) 
    qubits = [c_shift_index(q) for q in qubits]
    @ccall libquest.applyMultiVarPhaseFunc(qureg::Qureg, qubits::Ptr{Cint}, numQubitsPerReg::Ptr{Cint}, numRegs::Cint, encoding::bitEncoding, coeffs::Ptr{Cdouble}, exponents::Ptr{Cdouble}, numTermsPerReg::Ptr{Cint})::Cvoid
end

function applyMultiVarPhaseFuncOverrides(qureg, qubits, numQubitsPerReg, numRegs, encoding, coeffs, exponents, numTermsPerReg, overrideInds, overridePhases, numOverrides) 
    test_qubit_present(qureg, qubits) 
    test_qubit_present(qureg, overrideInds)
    qubits = [c_shift_index(q) for q in qubits]
    overrideInds = [c_shift_index(oi) for oi in overrideInds]
    @ccall libquest.applyMultiVarPhaseFuncOverrides(qureg::Qureg, qubits::Ptr{Cint}, numQubitsPerReg::Ptr{Cint}, numRegs::Cint, encoding::bitEncoding, coeffs::Ptr{Cdouble}, exponents::Ptr{Cdouble}, numTermsPerReg::Ptr{Cint}, overrideInds::Ptr{Clonglong}, overridePhases::Ptr{Cdouble}, numOverrides::Cint)::Cvoid
end

function applyNamedPhaseFunc(qureg, qubits, numQubitsPerReg, numRegs, encoding, functionNameCode)
    test_qubit_present(qureg, qubits)
    qubits = [c_shift_index(q) for q in qubits]
    @ccall libquest.applyNamedPhaseFunc(qureg::Qureg, qubits::Ptr{Cint}, numQubitsPerReg::Ptr{Cint}, numRegs::Cint, encoding::bitEncoding, functionNameCode::phaseFunc)::Cvoid
end

function applyNamedPhaseFuncOverrides(qureg, qubits, numQubitsPerReg, numRegs, encoding, functionNameCode, overrideInds, overridePhases, numOverrides)
    test_qubit_present(qureg, qubits)
    test_qubit_present(qureg, overrideInds)
    qubits = [c_shift_index(q) for q in qubits]
    overrideInds = [c_shift_index(oi) for oi in overrideInds]
    @ccall libquest.applyNamedPhaseFuncOverrides(qureg::Qureg, qubits::Ptr{Cint}, numQubitsPerReg::Ptr{Cint}, numRegs::Cint, encoding::bitEncoding, functionNameCode::phaseFunc, overrideInds::Ptr{Clonglong}, overridePhases::Ptr{Cdouble}, numOverrides::Cint)::Cvoid
end

function applyParamNamedPhaseFunc(qureg, qubits, numQubitsPerReg, numRegs, encoding, functionNameCode, params, numParams)
    test_qubit_present(qureg, qubits)
    qubits = [c_shift_index(q) for q in qubits]
    @ccall libquest.applyParamNamedPhaseFunc(qureg::Qureg, qubits::Ptr{Cint}, numQubitsPerReg::Ptr{Cint}, numRegs::Cint, encoding::bitEncoding, functionNameCode::phaseFunc, params::Ptr{Cdouble}, numParams::Cint)::Cvoid
end

function applyParamNamedPhaseFuncOverrides(qureg, qubits, numQubitsPerReg, numRegs, encoding, functionNameCode, params, numParams, overrideInds, overridePhases, numOverrides)
    test_qubit_present(qureg, qubits)
    test_qubit_present(qureg, overrideInds)
    qubits = [c_shift_index(q) for q in qubits]
    overrideInds = [c_shift_index(oi) for oi in overrideInds]
    @ccall libquest.applyParamNamedPhaseFuncOverrides(qureg::Qureg, qubits::Ptr{Cint}, numQubitsPerReg::Ptr{Cint}, numRegs::Cint, encoding::bitEncoding, functionNameCode::phaseFunc, params::Ptr{Cdouble}, numParams::Cint, overrideInds::Ptr{Clonglong}, overridePhases::Ptr{Cdouble}, numOverrides::Cint)::Cvoid
end

function applyFullQFT(qureg)
    @ccall libquest.applyFullQFT(qureg::Qureg)::Cvoid
end

function applyQFT(qureg, qubits)
    test_qubit_present(qureg, qubits)
    qubits = [c_shift_index(q) for q in qubits]
    numQubits = length(numQubits)
    @ccall libquest.applyQFT(qureg::Qureg, qubits::Ptr{Cint}, numQubits::Cint)::Cvoid
end

function applyProjector(qureg, qubit, outcome)
    test_qubit_present(qureg, qubit)
    qubit = c_shift_index(qubit)
    @ccall libquest.applyProjector(qureg::Qureg, qubit::Cint, outcome::Cint)::Cvoid
end

