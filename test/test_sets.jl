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
