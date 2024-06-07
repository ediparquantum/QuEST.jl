using Pkg
Pkg.activate(".")
using Test

include("src/QuEST.jl")
using ..QuEST




function test_PauliHamil()
    env = createQuESTEnv()
    paulies = QuEST.jpauliOpType
    numQubits = rand(1:12)
    numSumTerms = rand(1:20)
    hamil = createPauliHamil(numQubits, numSumTerms)
    @test hamil.numSumTerms == numSumTerms
    @test hamil.numQubits   == numQubits
    destroyPauliHamil(hamil)
    hamil = createPauliHamil(numQubits, numSumTerms)
    codes =  rand(instances(paulies), numQubits*numSumTerms)
    coeffs = rand(Float64, numSumTerms)
    initPauliHamil(hamil, coeffs, codes)
    hamil_codes = unsafe_wrap(Vector{pauliOpType}, hamil.pauliCodes, numQubits*numSumTerms)
    hamil_coeffs = unsafe_wrap(Vector{Float64}, hamil.termCoeffs, numSumTerms)
    for ind =1:numSumTerms
        @test coeffs[ind] ≈ hamil_coeffs[ind]  atol = tolerance
        for qubit = 1:numQubits
            @test codes[(ind-1)*numQubits+qubit] == hamil_codes[(ind-1)*numQubits+qubit]
        end
    end
    destroyPauliHamil(hamil)

    

    f = open("hamil","w")
    
    for row=1:numSumTerms
        write(f, string(coeffs[row]))
        write(f, " ")
        for qubit =1:numQubits
            write(f, string(Int(codes[(row-1)*numQubits+qubit])))
            write(f, " ")
        end
        write(f, "\n")
    end
    close(f)

    hamil = createPauliHamilFromFile("hamil")
    hamil_codes = unsafe_wrap(Vector{pauliOpType}, hamil.pauliCodes, numQubits*numSumTerms)
    hamil_coeffs = unsafe_wrap(Vector{Float64}, hamil.termCoeffs, numSumTerms)
    for ind=1:numSumTerms
        @test coeffs[ind] ≈ hamil_coeffs[ind] atol = tolerance
        for qubit = 1:numQubits
            @test codes[(ind-1)*numQubits+qubit] == hamil_codes[(ind-1)*numQubits+qubit]
        end
    end

    rm("hamil")

    destroyQuESTEnv(env)
end