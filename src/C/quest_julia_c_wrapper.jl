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
    PAULI_I = 0
    PAULI_X = 1
    PAULI_Y = 2
    PAULI_Z = 3
end

struct Complex
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

struct Vector
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
    @ccall libquest.setDiagonalOpElems(op::DiagonalOp, startInd::Clonglong, real::Ptr{Cdouble}, imag::Ptr{Cdouble}, numElems::Clonglong)::Cvoid
end

function applyDiagonalOp(qureg, op)
    @ccall libquest.applyDiagonalOp(qureg::Qureg, op::DiagonalOp)::Cvoid
end

function calcExpecDiagonalOp(qureg, op)
    @ccall libquest.calcExpecDiagonalOp(qureg::Qureg, op::DiagonalOp)::Complex
end

function createSubDiagonalOp(numQubits)
    @ccall libquest.createSubDiagonalOp(numQubits::Cint)::SubDiagonalOp
end

function destroySubDiagonalOp(op)
    @ccall libquest.destroySubDiagonalOp(op::SubDiagonalOp)::Cvoid
end

function diagonalUnitary(qureg, targets, numTargets, op)
    @ccall libquest.diagonalUnitary(qureg::Qureg, targets::Ptr{Cint}, numTargets::Cint, op::SubDiagonalOp)::Cvoid
end

function applyGateSubDiagonalOp(qureg, targets, numTargets, op)
    @ccall libquest.applyGateSubDiagonalOp(qureg::Qureg, targets::Ptr{Cint}, numTargets::Cint, op::SubDiagonalOp)::Cvoid
end

function applySubDiagonalOp(qureg, targets, numTargets, op)
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
    @ccall libquest.setAmps(qureg::Qureg, startInd::Clonglong, reals::Ptr{Cdouble}, imags::Ptr{Cdouble}, numAmps::Clonglong)::Cvoid
end

function setDensityAmps(qureg, startRow, startCol, reals, imags, numAmps)
    @ccall libquest.setDensityAmps(qureg::Qureg, startRow::Clonglong, startCol::Clonglong, reals::Ptr{Cdouble}, imags::Ptr{Cdouble}, numAmps::Clonglong)::Cvoid
end

function setQuregToPauliHamil(qureg, hamil)
    @ccall libquest.setQuregToPauliHamil(qureg::Qureg, hamil::PauliHamil)::Cvoid
end

function cloneQureg(targetQureg, copyQureg)
    @ccall libquest.cloneQureg(targetQureg::Qureg, copyQureg::Qureg)::Cvoid
end

function phaseShift(qureg, targetQubit, angle)
    @ccall libquest.phaseShift(qureg::Qureg, targetQubit::Cint, angle::Cdouble)::Cvoid
end

function controlledPhaseShift(qureg, idQubit1, idQubit2, angle)
    @ccall libquest.controlledPhaseShift(qureg::Qureg, idQubit1::Cint, idQubit2::Cint, angle::Cdouble)::Cvoid
end

function multiControlledPhaseShift(qureg, controlQubits, numControlQubits, angle)
    @ccall libquest.multiControlledPhaseShift(qureg::Qureg, controlQubits::Ptr{Cint}, numControlQubits::Cint, angle::Cdouble)::Cvoid
end

function controlledPhaseFlip(qureg, idQubit1, idQubit2)
    @ccall libquest.controlledPhaseFlip(qureg::Qureg, idQubit1::Cint, idQubit2::Cint)::Cvoid
end

function multiControlledPhaseFlip(qureg, controlQubits, numControlQubits)
    @ccall libquest.multiControlledPhaseFlip(qureg::Qureg, controlQubits::Ptr{Cint}, numControlQubits::Cint)::Cvoid
end

function sGate(qureg, targetQubit)
    @ccall libquest.sGate(qureg::Qureg, targetQubit::Cint)::Cvoid
end

function tGate(qureg, targetQubit)
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
    @ccall libquest.copySubstateToGPU(qureg::Qureg, startInd::Clonglong, numAmps::Clonglong)::Cvoid
end

function copySubstateFromGPU(qureg, startInd, numAmps)
    @ccall libquest.copySubstateFromGPU(qureg::Qureg, startInd::Clonglong, numAmps::Clonglong)::Cvoid
end

function getAmp(qureg, index)
    @ccall libquest.getAmp(qureg::Qureg, index::Clonglong)::Complex
end

function getRealAmp(qureg, index)
    @ccall libquest.getRealAmp(qureg::Qureg, index::Clonglong)::Cdouble
end

function getImagAmp(qureg, index)
    @ccall libquest.getImagAmp(qureg::Qureg, index::Clonglong)::Cdouble
end

function getProbAmp(qureg, index)
    @ccall libquest.getProbAmp(qureg::Qureg, index::Clonglong)::Cdouble
end

function getDensityAmp(qureg, row, col)
    @ccall libquest.getDensityAmp(qureg::Qureg, row::Clonglong, col::Clonglong)::Complex
end

function calcTotalProb(qureg)
    @ccall libquest.calcTotalProb(qureg::Qureg)::Cdouble
end

function compactUnitary(qureg, targetQubit, alpha, beta)
    @ccall libquest.compactUnitary(qureg::Qureg, targetQubit::Cint, alpha::Complex, beta::Complex)::Cvoid
end

function unitary(qureg, targetQubit, u)
    @ccall libquest.unitary(qureg::Qureg, targetQubit::Cint, u::ComplexMatrix2)::Cvoid
end

function rotateX(qureg, rotQubit, angle)
    @ccall libquest.rotateX(qureg::Qureg, rotQubit::Cint, angle::Cdouble)::Cvoid
end

function rotateY(qureg, rotQubit, angle)
    @ccall libquest.rotateY(qureg::Qureg, rotQubit::Cint, angle::Cdouble)::Cvoid
end

function rotateZ(qureg, rotQubit, angle)
    @ccall libquest.rotateZ(qureg::Qureg, rotQubit::Cint, angle::Cdouble)::Cvoid
end

function rotateAroundAxis(qureg, rotQubit, angle, axis)
    @ccall libquest.rotateAroundAxis(qureg::Qureg, rotQubit::Cint, angle::Cdouble, axis::Vector)::Cvoid
end

function controlledRotateX(qureg, controlQubit, targetQubit, angle)
    @ccall libquest.controlledRotateX(qureg::Qureg, controlQubit::Cint, targetQubit::Cint, angle::Cdouble)::Cvoid
end

function controlledRotateY(qureg, controlQubit, targetQubit, angle)
    @ccall libquest.controlledRotateY(qureg::Qureg, controlQubit::Cint, targetQubit::Cint, angle::Cdouble)::Cvoid
end

function controlledRotateZ(qureg, controlQubit, targetQubit, angle)
    @ccall libquest.controlledRotateZ(qureg::Qureg, controlQubit::Cint, targetQubit::Cint, angle::Cdouble)::Cvoid
end

function controlledRotateAroundAxis(qureg, controlQubit, targetQubit, angle, axis)
    @ccall libquest.controlledRotateAroundAxis(qureg::Qureg, controlQubit::Cint, targetQubit::Cint, angle::Cdouble, axis::Vector)::Cvoid
end

function controlledCompactUnitary(qureg, controlQubit, targetQubit, alpha, beta)
    @ccall libquest.controlledCompactUnitary(qureg::Qureg, controlQubit::Cint, targetQubit::Cint, alpha::Complex, beta::Complex)::Cvoid
end

function controlledUnitary(qureg, controlQubit, targetQubit, u)
    @ccall libquest.controlledUnitary(qureg::Qureg, controlQubit::Cint, targetQubit::Cint, u::ComplexMatrix2)::Cvoid
end

function multiControlledUnitary(qureg, controlQubits, numControlQubits, targetQubit, u)
    @ccall libquest.multiControlledUnitary(qureg::Qureg, controlQubits::Ptr{Cint}, numControlQubits::Cint, targetQubit::Cint, u::ComplexMatrix2)::Cvoid
end

function pauliX(qureg, targetQubit)
    @ccall libquest.pauliX(qureg::Qureg, targetQubit::Cint)::Cvoid
end

function pauliY(qureg, targetQubit)
    @ccall libquest.pauliY(qureg::Qureg, targetQubit::Cint)::Cvoid
end

function pauliZ(qureg, targetQubit)
    @ccall libquest.pauliZ(qureg::Qureg, targetQubit::Cint)::Cvoid
end

function hadamard(qureg, targetQubit)
    @ccall libquest.hadamard(qureg::Qureg, targetQubit::Cint)::Cvoid
end

function controlledNot(qureg, controlQubit, targetQubit)
    @ccall libquest.controlledNot(qureg::Qureg, controlQubit::Cint, targetQubit::Cint)::Cvoid
end

function multiControlledMultiQubitNot(qureg, ctrls, numCtrls, targs, numTargs)
    @ccall libquest.multiControlledMultiQubitNot(qureg::Qureg, ctrls::Ptr{Cint}, numCtrls::Cint, targs::Ptr{Cint}, numTargs::Cint)::Cvoid
end

function multiQubitNot(qureg, targs, numTargs)
    @ccall libquest.multiQubitNot(qureg::Qureg, targs::Ptr{Cint}, numTargs::Cint)::Cvoid
end

function controlledPauliY(qureg, controlQubit, targetQubit)
    @ccall libquest.controlledPauliY(qureg::Qureg, controlQubit::Cint, targetQubit::Cint)::Cvoid
end

function calcProbOfOutcome(qureg, measureQubit, outcome)
    @ccall libquest.calcProbOfOutcome(qureg::Qureg, measureQubit::Cint, outcome::Cint)::Cdouble
end

function calcProbOfAllOutcomes(outcomeProbs, qureg, qubits, numQubits)
    @ccall libquest.calcProbOfAllOutcomes(outcomeProbs::Ptr{Cdouble}, qureg::Qureg, qubits::Ptr{Cint}, numQubits::Cint)::Cvoid
end

function collapseToOutcome(qureg, measureQubit, outcome)
    @ccall libquest.collapseToOutcome(qureg::Qureg, measureQubit::Cint, outcome::Cint)::Cdouble
end

function measure(qureg, measureQubit)
    @ccall libquest.measure(qureg::Qureg, measureQubit::Cint)::Cint
end

function measureWithStats(qureg, measureQubit, outcomeProb)
    @ccall libquest.measureWithStats(qureg::Qureg, measureQubit::Cint, outcomeProb::Ptr{Cdouble})::Cint
end

function calcInnerProduct(bra, ket)
    @ccall libquest.calcInnerProduct(bra::Qureg, ket::Qureg)::Complex
end

function calcDensityInnerProduct(rho1, rho2)
    @ccall libquest.calcDensityInnerProduct(rho1::Qureg, rho2::Qureg)::Cdouble
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
    @ccall libquest.mixDephasing(qureg::Qureg, targetQubit::Cint, prob::Cdouble)::Cvoid
end

function mixTwoQubitDephasing(qureg, qubit1, qubit2, prob)
    @ccall libquest.mixTwoQubitDephasing(qureg::Qureg, qubit1::Cint, qubit2::Cint, prob::Cdouble)::Cvoid
end

function mixDepolarising(qureg, targetQubit, prob)
    @ccall libquest.mixDepolarising(qureg::Qureg, targetQubit::Cint, prob::Cdouble)::Cvoid
end

function mixDamping(qureg, targetQubit, prob)
    @ccall libquest.mixDamping(qureg::Qureg, targetQubit::Cint, prob::Cdouble)::Cvoid
end

function mixTwoQubitDepolarising(qureg, qubit1, qubit2, prob)
    @ccall libquest.mixTwoQubitDepolarising(qureg::Qureg, qubit1::Cint, qubit2::Cint, prob::Cdouble)::Cvoid
end

function mixPauli(qureg, targetQubit, probX, probY, probZ)
    @ccall libquest.mixPauli(qureg::Qureg, targetQubit::Cint, probX::Cdouble, probY::Cdouble, probZ::Cdouble)::Cvoid
end

function mixDensityMatrix(combineQureg, prob, otherQureg)
    @ccall libquest.mixDensityMatrix(combineQureg::Qureg, prob::Cdouble, otherQureg::Qureg)::Cvoid
end

function calcPurity(qureg)
    @ccall libquest.calcPurity(qureg::Qureg)::Cdouble
end

function calcFidelity(qureg, pureState)
    @ccall libquest.calcFidelity(qureg::Qureg, pureState::Qureg)::Cdouble
end

function swapGate(qureg, qubit1, qubit2)
    @ccall libquest.swapGate(qureg::Qureg, qubit1::Cint, qubit2::Cint)::Cvoid
end

function sqrtSwapGate(qureg, qb1, qb2)
    @ccall libquest.sqrtSwapGate(qureg::Qureg, qb1::Cint, qb2::Cint)::Cvoid
end

function multiStateControlledUnitary(qureg, controlQubits, controlState, numControlQubits, targetQubit, u)
    @ccall libquest.multiStateControlledUnitary(qureg::Qureg, controlQubits::Ptr{Cint}, controlState::Ptr{Cint}, numControlQubits::Cint, targetQubit::Cint, u::ComplexMatrix2)::Cvoid
end

function multiRotateZ(qureg, qubits, numQubits, angle)
    @ccall libquest.multiRotateZ(qureg::Qureg, qubits::Ptr{Cint}, numQubits::Cint, angle::Cdouble)::Cvoid
end

function multiRotatePauli(qureg, targetQubits, targetPaulis, numTargets, angle)
    @ccall libquest.multiRotatePauli(qureg::Qureg, targetQubits::Ptr{Cint}, targetPaulis::Ptr{pauliOpType}, numTargets::Cint, angle::Cdouble)::Cvoid
end

function multiControlledMultiRotateZ(qureg, controlQubits, numControls, targetQubits, numTargets, angle)
    @ccall libquest.multiControlledMultiRotateZ(qureg::Qureg, controlQubits::Ptr{Cint}, numControls::Cint, targetQubits::Ptr{Cint}, numTargets::Cint, angle::Cdouble)::Cvoid
end

function multiControlledMultiRotatePauli(qureg, controlQubits, numControls, targetQubits, targetPaulis, numTargets, angle)
    @ccall libquest.multiControlledMultiRotatePauli(qureg::Qureg, controlQubits::Ptr{Cint}, numControls::Cint, targetQubits::Ptr{Cint}, targetPaulis::Ptr{pauliOpType}, numTargets::Cint, angle::Cdouble)::Cvoid
end

function calcExpecPauliProd(qureg, targetQubits, pauliCodes, numTargets, workspace)
    @ccall libquest.calcExpecPauliProd(qureg::Qureg, targetQubits::Ptr{Cint}, pauliCodes::Ptr{pauliOpType}, numTargets::Cint, workspace::Qureg)::Cdouble
end

function calcExpecPauliSum(qureg, allPauliCodes, termCoeffs, numSumTerms, workspace)
    @ccall libquest.calcExpecPauliSum(qureg::Qureg, allPauliCodes::Ptr{pauliOpType}, termCoeffs::Ptr{Cdouble}, numSumTerms::Cint, workspace::Qureg)::Cdouble
end

function calcExpecPauliHamil(qureg, hamil, workspace)
    @ccall libquest.calcExpecPauliHamil(qureg::Qureg, hamil::PauliHamil, workspace::Qureg)::Cdouble
end

function twoQubitUnitary(qureg, targetQubit1, targetQubit2, u)
    @ccall libquest.twoQubitUnitary(qureg::Qureg, targetQubit1::Cint, targetQubit2::Cint, u::ComplexMatrix4)::Cvoid
end

function controlledTwoQubitUnitary(qureg, controlQubit, targetQubit1, targetQubit2, u)
    @ccall libquest.controlledTwoQubitUnitary(qureg::Qureg, controlQubit::Cint, targetQubit1::Cint, targetQubit2::Cint, u::ComplexMatrix4)::Cvoid
end

function multiControlledTwoQubitUnitary(qureg, controlQubits, numControlQubits, targetQubit1, targetQubit2, u)
    @ccall libquest.multiControlledTwoQubitUnitary(qureg::Qureg, controlQubits::Ptr{Cint}, numControlQubits::Cint, targetQubit1::Cint, targetQubit2::Cint, u::ComplexMatrix4)::Cvoid
end

function multiQubitUnitary(qureg, targs, numTargs, u)
    @ccall libquest.multiQubitUnitary(qureg::Qureg, targs::Ptr{Cint}, numTargs::Cint, u::ComplexMatrixN)::Cvoid
end

function controlledMultiQubitUnitary(qureg, ctrl, targs, numTargs, u)
    @ccall libquest.controlledMultiQubitUnitary(qureg::Qureg, ctrl::Cint, targs::Ptr{Cint}, numTargs::Cint, u::ComplexMatrixN)::Cvoid
end

function multiControlledMultiQubitUnitary(qureg, ctrls, numCtrls, targs, numTargs, u)
    @ccall libquest.multiControlledMultiQubitUnitary(qureg::Qureg, ctrls::Ptr{Cint}, numCtrls::Cint, targs::Ptr{Cint}, numTargs::Cint, u::ComplexMatrixN)::Cvoid
end

function mixKrausMap(qureg, target, ops, numOps)
    @ccall libquest.mixKrausMap(qureg::Qureg, target::Cint, ops::Ptr{ComplexMatrix2}, numOps::Cint)::Cvoid
end

function mixTwoQubitKrausMap(qureg, target1, target2, ops, numOps)
    @ccall libquest.mixTwoQubitKrausMap(qureg::Qureg, target1::Cint, target2::Cint, ops::Ptr{ComplexMatrix4}, numOps::Cint)::Cvoid
end

function mixMultiQubitKrausMap(qureg, targets, numTargets, ops, numOps)
    @ccall libquest.mixMultiQubitKrausMap(qureg::Qureg, targets::Ptr{Cint}, numTargets::Cint, ops::Ptr{ComplexMatrixN}, numOps::Cint)::Cvoid
end

function mixNonTPKrausMap(qureg, target, ops, numOps)
    @ccall libquest.mixNonTPKrausMap(qureg::Qureg, target::Cint, ops::Ptr{ComplexMatrix2}, numOps::Cint)::Cvoid
end

function mixNonTPTwoQubitKrausMap(qureg, target1, target2, ops, numOps)
    @ccall libquest.mixNonTPTwoQubitKrausMap(qureg::Qureg, target1::Cint, target2::Cint, ops::Ptr{ComplexMatrix4}, numOps::Cint)::Cvoid
end

function mixNonTPMultiQubitKrausMap(qureg, targets, numTargets, ops, numOps)
    @ccall libquest.mixNonTPMultiQubitKrausMap(qureg::Qureg, targets::Ptr{Cint}, numTargets::Cint, ops::Ptr{ComplexMatrixN}, numOps::Cint)::Cvoid
end

function calcHilbertSchmidtDistance(a, b)
    @ccall libquest.calcHilbertSchmidtDistance(a::Qureg, b::Qureg)::Cdouble
end

function setWeightedQureg(fac1, qureg1, fac2, qureg2, facOut, out)
    @ccall libquest.setWeightedQureg(fac1::Complex, qureg1::Qureg, fac2::Complex, qureg2::Qureg, facOut::Complex, out::Qureg)::Cvoid
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
    @ccall libquest.applyMatrix2(qureg::Qureg, targetQubit::Cint, u::ComplexMatrix2)::Cvoid
end

function applyMatrix4(qureg, targetQubit1, targetQubit2, u)
    @ccall libquest.applyMatrix4(qureg::Qureg, targetQubit1::Cint, targetQubit2::Cint, u::ComplexMatrix4)::Cvoid
end

function applyMatrixN(qureg, targs, numTargs, u)
    @ccall libquest.applyMatrixN(qureg::Qureg, targs::Ptr{Cint}, numTargs::Cint, u::ComplexMatrixN)::Cvoid
end

function applyGateMatrixN(qureg, targs, numTargs, u)
    @ccall libquest.applyGateMatrixN(qureg::Qureg, targs::Ptr{Cint}, numTargs::Cint, u::ComplexMatrixN)::Cvoid
end

function applyMultiControlledGateMatrixN(qureg, ctrls, numCtrls, targs, numTargs, m)
    @ccall libquest.applyMultiControlledGateMatrixN(qureg::Qureg, ctrls::Ptr{Cint}, numCtrls::Cint, targs::Ptr{Cint}, numTargs::Cint, m::ComplexMatrixN)::Cvoid
end

function applyMultiControlledMatrixN(qureg, ctrls, numCtrls, targs, numTargs, u)
    @ccall libquest.applyMultiControlledMatrixN(qureg::Qureg, ctrls::Ptr{Cint}, numCtrls::Cint, targs::Ptr{Cint}, numTargs::Cint, u::ComplexMatrixN)::Cvoid
end

function invalidQuESTInputError(errMsg, errFunc)
    @ccall libquest.invalidQuESTInputError(errMsg::Cstring, errFunc::Cstring)::Cvoid
end

function applyPhaseFunc(qureg, qubits, numQubits, encoding, coeffs, exponents, numTerms)
    @ccall libquest.applyPhaseFunc(qureg::Qureg, qubits::Ptr{Cint}, numQubits::Cint, encoding::bitEncoding, coeffs::Ptr{Cdouble}, exponents::Ptr{Cdouble}, numTerms::Cint)::Cvoid
end

function applyPhaseFuncOverrides(qureg, qubits, numQubits, encoding, coeffs, exponents, numTerms, overrideInds, overridePhases, numOverrides)
    @ccall libquest.applyPhaseFuncOverrides(qureg::Qureg, qubits::Ptr{Cint}, numQubits::Cint, encoding::bitEncoding, coeffs::Ptr{Cdouble}, exponents::Ptr{Cdouble}, numTerms::Cint, overrideInds::Ptr{Clonglong}, overridePhases::Ptr{Cdouble}, numOverrides::Cint)::Cvoid
end

function applyMultiVarPhaseFunc(qureg, qubits, numQubitsPerReg, numRegs, encoding, coeffs, exponents, numTermsPerReg)
    @ccall libquest.applyMultiVarPhaseFunc(qureg::Qureg, qubits::Ptr{Cint}, numQubitsPerReg::Ptr{Cint}, numRegs::Cint, encoding::bitEncoding, coeffs::Ptr{Cdouble}, exponents::Ptr{Cdouble}, numTermsPerReg::Ptr{Cint})::Cvoid
end

function applyMultiVarPhaseFuncOverrides(qureg, qubits, numQubitsPerReg, numRegs, encoding, coeffs, exponents, numTermsPerReg, overrideInds, overridePhases, numOverrides)
    @ccall libquest.applyMultiVarPhaseFuncOverrides(qureg::Qureg, qubits::Ptr{Cint}, numQubitsPerReg::Ptr{Cint}, numRegs::Cint, encoding::bitEncoding, coeffs::Ptr{Cdouble}, exponents::Ptr{Cdouble}, numTermsPerReg::Ptr{Cint}, overrideInds::Ptr{Clonglong}, overridePhases::Ptr{Cdouble}, numOverrides::Cint)::Cvoid
end

function applyNamedPhaseFunc(qureg, qubits, numQubitsPerReg, numRegs, encoding, functionNameCode)
    @ccall libquest.applyNamedPhaseFunc(qureg::Qureg, qubits::Ptr{Cint}, numQubitsPerReg::Ptr{Cint}, numRegs::Cint, encoding::bitEncoding, functionNameCode::phaseFunc)::Cvoid
end

function applyNamedPhaseFuncOverrides(qureg, qubits, numQubitsPerReg, numRegs, encoding, functionNameCode, overrideInds, overridePhases, numOverrides)
    @ccall libquest.applyNamedPhaseFuncOverrides(qureg::Qureg, qubits::Ptr{Cint}, numQubitsPerReg::Ptr{Cint}, numRegs::Cint, encoding::bitEncoding, functionNameCode::phaseFunc, overrideInds::Ptr{Clonglong}, overridePhases::Ptr{Cdouble}, numOverrides::Cint)::Cvoid
end

function applyParamNamedPhaseFunc(qureg, qubits, numQubitsPerReg, numRegs, encoding, functionNameCode, params, numParams)
    @ccall libquest.applyParamNamedPhaseFunc(qureg::Qureg, qubits::Ptr{Cint}, numQubitsPerReg::Ptr{Cint}, numRegs::Cint, encoding::bitEncoding, functionNameCode::phaseFunc, params::Ptr{Cdouble}, numParams::Cint)::Cvoid
end

function applyParamNamedPhaseFuncOverrides(qureg, qubits, numQubitsPerReg, numRegs, encoding, functionNameCode, params, numParams, overrideInds, overridePhases, numOverrides)
    @ccall libquest.applyParamNamedPhaseFuncOverrides(qureg::Qureg, qubits::Ptr{Cint}, numQubitsPerReg::Ptr{Cint}, numRegs::Cint, encoding::bitEncoding, functionNameCode::phaseFunc, params::Ptr{Cdouble}, numParams::Cint, overrideInds::Ptr{Clonglong}, overridePhases::Ptr{Cdouble}, numOverrides::Cint)::Cvoid
end

function applyFullQFT(qureg)
    @ccall libquest.applyFullQFT(qureg::Qureg)::Cvoid
end

function applyQFT(qureg, qubits, numQubits)
    @ccall libquest.applyQFT(qureg::Qureg, qubits::Ptr{Cint}, numQubits::Cint)::Cvoid
end

function applyProjector(qureg, qubit, outcome)
    @ccall libquest.applyProjector(qureg::Qureg, qubit::Cint, outcome::Cint)::Cvoid
end

