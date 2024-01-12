##################################################################
# Filename  : test_sets.jl
# Author    : Jonathan Miller
# Date      : 2024-01-09
# Aim       : aim_script
#           : The aim of this script is to provide the test sets
#           : for the QuEST.jl module.
##################################################################

println("Running tests on QuEST.jl")
println("Starting with structs")
@testset "Test: QuEST Environmnet struct" begin
    test_env()
end  

@testset "Test: QuEST Qureg struct (GPU disabled)" begin
    test_state_vec_qureg_no_gpu()
    test_density_matrix_qureg_no_gpu()
end

@testset "Test: Clone Qureg" begin
    test_createCloneQureg()
end

#= Tests for ComplexMatrixN and DiagonalOp are not implemented
@testset "Test: Complex MatrixN" begin
    test_ComplexMatrixN()
end

@testset "Test: Diagonal Op" begin
    test_DiagonalOp()
end
=#
@testset "Test: Pauli Hamil" begin
    test_PauliHamil()
end

@testset "Test: Compact Unitary" begin
    test_compactUnitary()
end

@testset "Test: Controlled Compact Unitary" begin
    test_controlledCompactUnitary()
end

@testset "Test: Controlled not gate" begin
    test_controlledNot()
end

@testset "Test: Controlled Pauli Y gate" begin
    test_controlledPauliY()
end

@testset "Test: Controlled phase flip gate" begin
    test_controlledPhaseFlip()
end

@testset "Test: Controlled phase shift gate" begin
    test_controlledPhaseShift()

end

@testset "Test: Controlled rotate around axis gate" begin
    test_controlledRotateAroundAxis()
end


@testset "Test: Controlled rotate X gate" begin
    test_controlledRotateX()
end

@testset "Test: Controlled rotate Y gate" begin
    test_controlledRotateY()
end

@testset "Test: Controlled rotate Z gate" begin
    test_controlledRotateZ()
end

@testset "Test: Controlled two qubit unitary gate" begin
    test_controlledTwoQubitUnitary()
end

@testset "Test: Controlled Unitary" begin
    test_controlledUnitary()
end

@testset "Test: Hadamard gate" begin
    test_hadamard()
end

@testset "Test: Multi Controlled Multi Qubit Unitary" begin
    test_multiControlledMultiQubitUnitary()
end


@testset "Test: Multi Controlled Phase Flip" begin
    test_multiControlledPhaseFlip()
end


@testset "Test: Multi Controlled Phase Shift" begin
    test_multiControlledPhaseShift()
end

@testset "Test: Multi Controlled Two Qubit Unitary" begin
    test_multiControlledTwoQubitUnitary()
end





@testset "Test: Sync" begin
    test_sync()
end