
@testset "test_env" begin
    test_env()
end

@testset "test_qureg" begin
    test_qureg()
end

@testset "test_qureg_density" begin
    test_qureg_density()
end

@testset "check_createCloneQureg" begin
    check_createCloneQureg()
end

@testset "check_ComplexMatrixN" begin
    check_ComplexMatrixN()
end

@testset "test_DiagonalOp" begin
    test_DiagonalOp()
end

@testset begin
    test_PauliHaiml()
end

@testset "test_PauliHaiml" begin
    test_compactUnitary()
end

@testset "test_compactUnitary" begin
    test_compactUnitary()
end

@testset "test_controlledCompactUnitary" begin
    test_controlledCompactUnitary()
end

@testset "test_controlledMultiQubitUnitary" begin
    test_controlledMultiQubitUnitary()
end

@testset "test_controlledNot" begin    
    test_controlledNot()
end

@testset begin    
    test_controlledPauliY()
end

@testset "test_controlledPauliY" begin    
    test_controlledPhaseFlip()
end

@testset "test_controlledPhaseFlip" begin
    test_controlledPhaseFlip()
end

@testset "test_controlledPhaseShift" begin    
    test_controlledPhaseShift()
end

@testset "test_controlledRotateAroundAxis" begin    
    test_controlledRotateAroundAxis()
end

@testset "test_controlledRotateX" begin    
    test_controlledRotateX()
end

@testset "test_controlledRotateY" begin    
    test_controlledRotateY()
end

@testset "test_controlledRotateZ" begin    
    test_controlledRotateZ()
end

@testset "test_controlledTwoQubitUnitary" begin    
    test_controlledTwoQubitUnitary()
end

@testset "test_controlledUnitary" begin    
    test_controlledUnitary()
end

@testset "test_hadamard" begin    
    test_hadamard()
end

@testset "test_multiControlledMultiQubitUnitary" begin    
    test_multiControlledMultiQubitUnitary()
end

@testset "test_multiControlledPhaseFlip" begin
    test_multiControlledPhaseFlip()
end

@testset "test_multiControlledPhaseShift" begin
    test_multiControlledPhaseShift()
end

@testset "test_multiControlledTwoQubitUnitary" begin
    test_multiControlledTwoQubitUnitary()
end

@testset "test_multiControlledUnitary" begin
    test_multiControlledUnitary()
end

@testset "test_multiQubitUnitary" begin
    test_multiQubitUnitary()
end

@testset "test_multiRotatePauli" begin
    test_multiRotatePauli()
end

@testset "test_multiRotateZ" begin
    test_multiRotateZ()
end

@testset "test_multiStateControlledUnitary" begin
    test_multiStateControlledUnitary()
end

@testset "test_pauliX" begin
    test_pauliX()
end

@testset "test_pauliY" begin
    test_pauliY()
end

@testset "test_pauliZ" begin
    test_pauliZ()
end

@testset "test_phaseShift" begin
    test_phaseShift()
end

@testset begin
    test_rotateAroundAxis()
end

@testset "test_rotateAroundAxis" begin
    test_rotateX()
end

@testset "test_rotateY" begin
    test_rotateY()
end

@testset "test_rotateZ" begin
    test_rotateZ()
end

@testset "test_sGate" begin
    test_sGate()
end

@testset "test_sqrtSwapGate" begin
    test_sqrtSwapGate()
end

@testset "test_swapGate" begin
    test_swapGate()
end

@testset "test_tGate" begin
    test_tGate()
end

@testset "test_twoQubitUnitary" begin
    test_twoQubitUnitary()
end

@testset "test_unitary" begin
    test_unitary()
end

    ##############
    ##############

@testset "test_collapseToOutcome" begin
    test_collapseToOutcome()
end

@testset "test_measure" begin
    test_measure()
end

@testset "test_measureWithStats" begin
    test_measureWithStats()
end

@testset "test_DiagonalOp" begin
    test_DiagonalOp()
end

@testset "test_applyMatrix2" begin
    test_applyMatrix2()
end

@testset "test_applyMatrix4" begin
    test_applyMatrix4()
end

@testset "test_applyMatrixN" begin
    test_applyMatrixN()
end

@testset "test_applyMultiControlledMatrixN" begin
    test_applyMultiControlledMatrixN()
end

@testset "test_applyPauliHamil" begin
    test_applyPauliHamil()
end

@testset "test_applyPauliSum" begin
    test_applyPauliSum()
end

@testset "test_applyTrotterCircuit" begin
    test_applyTrotterCircuit()
end

    #####################
    #####################
@testset "test_calcDensityInnerProduct" begin
    test_calcDensityInnerProduct()
end
@testset "test_calcExpecDiagonalOp" begin
    test_calcExpecDiagonalOp()
end
@testset "test_calcExpecPauliHamil" begin
    test_calcExpecPauliHamil()
end
@testset "test_calcExpecPauliProd" begin
    test_calcExpecPauliProd()
end
@testset "test_calcExpecPauliSum" begin
    test_calcExpecPauliSum()
end
@testset "test_calcFidelity" begin
    test_calcFidelity()
end
@testset "test_calcHilbertSchmidtDistance" begin
    test_calcHilbertSchmidtDistance()
end
@testset "test_calcInnerProduct" begin
    test_calcInnerProduct()
end
@testset "test_calcProbOfOutcome" begin
    test_calcProbOfOutcome()
end
@testset "test_calcPurity" begin
    test_calcPurity()
end
@testset "test_calcTotalProb" begin
    test_calcTotalProb()
end
@testset "test_amp_funcs" begin
    test_amp_funcs()
end
@testset "test_getDensityAmp" begin
    test_getDensityAmp()
end

    ####################
    ####################
@testset "test_cloneQureg" begin
    test_cloneQureg()
end
@testset "test_initBlankState" begin
    test_initBlankState()
end
@testset "test_initClassicalState" begin
    test_initClassicalState()
end
@testset "test_initPlusState" begin
    test_initPlusState()
end
@testset "test_initPureState" begin
    test_initPureState()
end
@testset "test_initStateFromAmps" begin
    test_initStateFromAmps()
end
@testset "test_initZeroState" begin
    test_initZeroState()
end
@testset "test_setAmps" begin
    test_setAmps()
end
@testset "test_setWeightedQureg" begin
    test_setWeightedQureg()
end

    #################
    #################

@testset "test_QASM" begin
    test_QASM()
end

    ###############
    ###############

@testset "test_GPU" begin
    test_GPU()
end
@testset "test_initDebugState" begin
    test_initDebugState()
end
@testset "test_reportPauliHamil" begin
    test_reportPauliHamil()
end

@testset "test_reportQuregParams" begin
    test_reportQuregParams()
end

@testset "test_reporQuESTEnv" begin
    test_reporQuESTEnv()
end

@testset "test_reportState" begin
    test_reportState()
end
@testset "test_seed" begin
    test_seed()
end
@testset "test_sync" begin
    test_sync()
end