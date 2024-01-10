
function test_env()
    env = createQuESTEnv()
    @test env |> typeof == QuESTEnv
    destroyQuESTEnv(env)
end

function test_state_vec_qureg_no_gpu()
    num_qubits = 1
    num_elements = 2^num_qubits
    env = createQuESTEnv()
    ψ = createQureg(num_qubits, env)
    @test env |> typeof == QuESTEnv
    @test ψ |> typeof == Qureg
    @test ψ.isDensityMatrix == 0 # 0 for pure state, 1 for density matrix
    @test ψ.numQubitsRepresented == num_qubits
    @test ψ.numQubitsInStateVec == num_qubits
    @test ψ.numAmpsPerChunk == num_elements   
    @test ψ.numAmpsTotal == num_elements 
    @test ψ.chunkId == 0 # Assume 0 as 1 chunk used
    @test ψ.numChunks == 1 # Assume 1 as single CPU running test chunk used
    @test unsafe_load_state_vec(ψ.stateVec,num_elements)[1].re == 1.0

    # For the field ψ.pairStateVec, it is temporary storage
    # for a chunk of the state vector received from another 
    # process in the MPI version. No test at this time.
    # Date 2024-01-09 - Jonathan Miller

    # For the field ψ.deviceStateVec, it is storage for 
    # wavefunction amplitudes in the GPU version.
    # This package did not enable GPU support at compilation
    # No test at this time.
    # Date 2024-01-09 - Jonathan Miller

    # For the field ψ.firstLevelReduction, it is storage for 
    # wavefunction amplitudes in the GPU version.
    # This package did not enable GPU support at compilation
    # No test at this time.
    # Date 2024-01-09 - Jonathan Miller

    # For the field ψ.secondLevelReduction, it is storage for
    # wavefunction amplitudes in the GPU version.
    # This package did not enable GPU support at compilation
    # No test at this time.
    # Date 2024-01-09 - Jonathan Miller

    # For the field ψ.cuStateVec, ψ.deviceCuStateVec and ψ.cuConfig, it is 
    # storage for wavefunction amplitues and config (copy of 
    # QuESTEnv's handle) in cuQuantum deployment.
    # This package did not enable GPU support at compilation
    # No test at this time.
    # Date 2024-01-09 - Jonathan Miller

    qasm_logger = unsafe_load(ψ.qasmLog)
    @test unsafe_string(qasm_logger.buffer) isa String
    @test qasm_logger.bufferSize isa Integer
    @test qasm_logger.bufferFill isa Integer
    @test qasm_logger.isLogging isa Integer

    destroyQureg(ψ, env)
    destroyQuESTEnv(env)
end


function test_density_matrix_qureg_no_gpu()
    num_qubits = 1
    num_elements = 2^num_qubits
    env = createQuESTEnv()
    ψ = createDensityQureg(num_qubits, env)
    @test env |> typeof == QuESTEnv
    @test ψ |> typeof == Qureg
    @test ψ.isDensityMatrix == 1 # 0 for pure state, 1 for density matrix
    @test ψ.numQubitsRepresented == num_qubits
    @test ψ.numQubitsInStateVec == 2*num_qubits
    @test ψ.numAmpsPerChunk == 2*num_elements   
    @test ψ.numAmpsTotal == 2*num_elements 
    @test ψ.chunkId == 0 # Assume 0 as 1 chunk used
    @test ψ.numChunks == 1 # Assume 1 as single CPU running test chunk used
    @test unsafe_load_state_vec(ψ.stateVec,num_elements)[1].re == 1.0

    # For the field ψ.pairStateVec, it is temporary storage
    # for a chunk of the state vector received from another 
    # process in the MPI version. No test at this time.
    # Date 2024-01-09 - Jonathan Miller

    # For the field ψ.deviceStateVec, it is storage for 
    # wavefunction amplitudes in the GPU version.
    # This package did not enable GPU support at compilation
    # No test at this time.
    # Date 2024-01-09 - Jonathan Miller

    # For the field ψ.firstLevelReduction, it is storage for 
    # wavefunction amplitudes in the GPU version.
    # This package did not enable GPU support at compilation
    # No test at this time.
    # Date 2024-01-09 - Jonathan Miller

    # For the field ψ.secondLevelReduction, it is storage for
    # wavefunction amplitudes in the GPU version.
    # This package did not enable GPU support at compilation
    # No test at this time.
    # Date 2024-01-09 - Jonathan Miller

    # For the field ψ.cuStateVec, ψ.deviceCuStateVec and ψ.cuConfig, it is 
    # storage for wavefunction amplitues and config (copy of 
    # QuESTEnv's handle) in cuQuantum deployment.
    # This package did not enable GPU support at compilation
    # No test at this time.
    # Date 2024-01-09 - Jonathan Miller

    qasm_logger = unsafe_load(ψ.qasmLog)
    @test unsafe_string(qasm_logger.buffer) isa String
    @test qasm_logger.bufferSize isa Integer
    @test qasm_logger.bufferFill isa Integer
    @test qasm_logger.isLogging isa Integer

    destroyQureg(ψ, env)
    destroyQuESTEnv(env)
end

function test_createCloneQureg()
    env= createQuESTEnv()
    num_qubits = 2
    qureg1 = createQureg(num_qubits, env)
    rotateX(qureg1, 1, rand())
    controlledNot(qureg1, 1, 2)
    qureg2 = createCloneQureg(qureg1, env)
    @test qureg1.numQubitsRepresented == qureg2.numQubitsRepresented
    @test qureg1.numAmpsTotal == qureg2.numAmpsTotal
    @test getNumAmps(qureg1) == getNumAmps(qureg2)
    @test getNumQubits(qureg1) == getNumQubits(qureg2)
    for i in Base.OneTo(qureg1.numAmpsTotal)
        @test getProbAmp(qureg1, i) == getProbAmp(qureg2, i)
        @test getRealAmp(qureg1, i) == getRealAmp(qureg2, i)
        @test getImagAmp(qureg1, i) == getImagAmp(qureg2, i)
    end
  destroyQuESTEnv(env)
end



function test_PauliHamil()
    env = createQuESTEnv()
    paulies = pauliOpType
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
    hamil_codes = unsafe_wrap(Base.Vector{pauliOpType}, hamil.pauliCodes, numQubits*numSumTerms)
    hamil_coeffs = unsafe_wrap(Base.Vector{Float64}, hamil.termCoeffs, numSumTerms)
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
    hamil_codes = unsafe_wrap(Base.Vector{pauliOpType}, hamil.pauliCodes, numQubits*numSumTerms)
    hamil_coeffs = unsafe_wrap(Base.Vector{Float64}, hamil.termCoeffs, numSumTerms)
    for ind=1:numSumTerms
        @test coeffs[ind] ≈ hamil_coeffs[ind] atol = tolerance
        for qubit = 1:numQubits
            @test codes[(ind-1)*numQubits+qubit] == hamil_codes[(ind-1)*numQubits+qubit]
        end
    end

    rm("hamil")

    destroyQuESTEnv(env)
end


function test_compactUnitary()
    env= createQuESTEnv()
    for t=1:10
        α=Base.Complex(rand(Cdouble),rand(Cdouble))
        β=Base.Complex(rand(Cdouble),rand(Cdouble))
        norm = sqrt(α*conj(α) + β*conj(β))
        α /= norm
        β /= norm
        numQubits = rand(1:12)
        target = rand(1:numQubits)
        qureg = createQureg(numQubits, env)
        α = QuEST.Complex(α.re, α.im)
        β = QuEST.Complex(β.re, β.im)
        compactUnitary(qureg, target, α, β)
        reals = unsafe_wrap(Base.Vector{Float64}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Base.Vector{Float64}, qureg.stateVec.imag, 2^numQubits)
        @test reals[1] ≈ α.real atol = tolerance
        @test imags[1] ≈ α.imag atol = tolerance
        @test reals[2^(target-1)+1] ≈ β.real atol = tolerance
        @test imags[2^(target-1)+1] ≈ β.imag atol = tolerance
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end