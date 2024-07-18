
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


function test_compactUnitary()
    env= createQuESTEnv()
    for t=1:10
        α=Complex(rand(Cdouble),rand(Cdouble))
        β=Complex(rand(Cdouble),rand(Cdouble))
        norm = sqrt(α*conj(α) + β*conj(β))
        α /= norm
        β /= norm
        numQubits = rand(1:12)
        target = rand(1:numQubits)
        qureg = createQureg(numQubits, env)
        α = QComplex(α.re, α.im)
        β = QComplex(β.re, β.im)
        compactUnitary(qureg, target, α, β)
        reals = unsafe_wrap(Vector{Float64}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{Float64}, qureg.stateVec.imag, 2^numQubits)
        @test reals[1] ≈ α.real atol = tolerance
        @test imags[1] ≈ α.imag atol = tolerance
        @test reals[2^(target-1)+1] ≈ β.real atol = tolerance
        @test imags[2^(target-1)+1] ≈ β.imag atol = tolerance
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end


function test_controlledCompactUnitary()
    env = createQuESTEnv()
    for t=1:10
        α=rand(Complex{Float64})
        β=rand(Complex{Float64})
        norm = sqrt(α*conj(α) + β*conj(β))
        α /= norm
        β /= norm
        @test abs(α)^2 + abs(β)^2 ≈ 1.0 atol = tolerance
        numQubits = 4
        target = 4
        control = 3
        qureg = createQureg(numQubits, env)
        pauliX(qureg, control)
        α = qComplex(α)
        β = qComplex(β)
        controlledCompactUnitary(qureg, control, target, α, β)
        reals = unsafe_wrap(Vector{Float64}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{Float64}, qureg.stateVec.imag, 2^numQubits)
        @test reals[2^(control-1)+1] ≈ α.real atol = tolerance
        @test imags[2^(control-1)+1] ≈ α.imag atol = tolerance
        @test reals[2^(control-1)+2^(target-1)+1] ≈ β.real atol = tolerance
        @test imags[2^(control-1)+2^(target-1)+1] ≈ β.imag atol = tolerance
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_controlledNot()
    env= createQuESTEnv()
    for t=1:10    
        numQubits = rand(2:12)
        target = rand(1:numQubits)
        control = rand(filter(x->x≠target, 1:numQubits))

        qureg = createQureg(numQubits, env)
        pauliX(qureg, control)
        controlledNot(qureg, control, target)
        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        one_ind=2^(control-1)+2^(target-1)+1
        for ind =1:2^numQubits
            if ind == one_ind
                @test state_vec[ind] ≈ Complex(1.0) atol = tolerance
            else
                @test state_vec[ind] ≈ Complex(0.0) atol = tolerance
            end
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_controlledPauliY()
    env= createQuESTEnv()
    for t=1:10
        
        numQubits = rand(2:12)
        target = rand(1:numQubits)
        control = rand(filter(x->x≠target, 1:numQubits))
        qureg = createQureg(numQubits, env)
        pauliX(qureg, control)
        controlledPauliY(qureg, control, target)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        one_ind=2^(control-1)+2^(target-1)+1
        for ind =1:2^numQubits
            if ind == one_ind
                @test state_vec[ind] ≈ Complex(1.0im) atol = tolerance
            else
                @test state_vec[ind] ≈ Complex(0.0) atol = tolerance
            end
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_controlledPhaseFlip()
    env= createQuESTEnv()
    for t=1:10
        
        numQubits = rand(2:12)
        target = rand(1:numQubits)
        control = rand(1:numQubits)
        control = rand(filter(x->x≠target, 1:numQubits))

        qureg = createQureg(numQubits, env)

        pauliX(qureg, control)
        pauliX(qureg, target)
        controlledPhaseFlip(qureg, control, target)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        one_ind=2^(control-1)+2^(target-1)+1
        for ind =1:2^numQubits
            if ind == one_ind
                @test state_vec[ind] ≈ Complex(-1.0) atol = tolerance
            else
                @test state_vec[ind] ≈ Complex(0.0) atol = tolerance
            end
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end


function test_controlledPhaseShift()
    env= createQuESTEnv()
    for t=1:10
        
        numQubits = rand(2:12)
        target = rand(1:numQubits)
        control = rand(filter(x->x≠target, 1:numQubits))
        θ = rand(Float64)
        qureg = createQureg(numQubits, env)

        pauliX(qureg, control)
        pauliX(qureg, target)
        controlledPhaseShift(qureg, control, target, θ)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        one_ind=2^(control-1)+2^(target-1)+1
        for ind =1:2^numQubits
            if ind == one_ind
                @test state_vec[ind] ≈ ℯ^(im*θ) atol = tolerance
            else
                @test state_vec[ind] ≈ 0 atol = tolerance
            end
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end


function test_controlledRotateAroundAxis()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(2:12)
        target = rand(1:numQubits)
        control = rand(filter(x->x≠target, 1:numQubits))
        θ = rand(Float64)
        n̂ = rand(Float64, 3) 
        n̂ /= norm(n̂)
        U = ℯ^(-im*θ*(n̂[1]*X+n̂[2]*Y+n̂[3]*Z)/2)
        qureg = createQureg(numQubits, env)
        v = [1.0, 0.0]
        v = U*v
        pauliX(qureg, control)
        
        controlledRotateAroundAxis(qureg, control, target, θ, qVector(n̂))

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        one_ind=2^(control-1)+2^(target-1)+1
        zero_ind=2^(control-1)+1
        for ind =1:2^numQubits
            if ind == one_ind
                @test state_vec[ind] ≈ v[2] atol = tolerance
            elseif ind == zero_ind
                @test state_vec[ind] ≈ v[1] atol = tolerance
            else
                @test state_vec[ind] ≈ Complex(0,0) atol = tolerance
            end
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_controlledRotateAroundAxis()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(2:12)
        target = rand(1:numQubits)
        control = rand(filter(x->x≠target, 1:numQubits))
        θ = rand(Float64)
        n̂ = rand(Float64, 3) 
        n̂ /= norm(n̂)
        U = ℯ^(-im*θ*(n̂[1]*X+n̂[2]*Y+n̂[3]*Z)/2)
        qureg = createQureg(numQubits, env)
        v = [1.0, 0.0]
        v = U*v
        pauliX(qureg, control)
        
        controlledRotateAroundAxis(qureg, control, target, θ, qVector(n̂))

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        one_ind=2^(control-1)+2^(target-1)+1
        zero_ind=2^(control-1)+1
        for ind =1:2^numQubits
            if ind == one_ind
                @test state_vec[ind] ≈ v[2] atol = tolerance
            elseif ind == zero_ind
                @test state_vec[ind] ≈ v[1] atol = tolerance
            else
                @test state_vec[ind] ≈ Complex(0,0) atol = tolerance
            end
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_controlledRotateAroundAxis()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(2:12)
        target = rand(1:numQubits)
        control = rand(filter(x->x≠target, 1:numQubits))
        θ = rand(Float64)
        n̂ = rand(Float64, 3) 
        n̂ /= norm(n̂)
        U = ℯ^(-im*θ*(n̂[1]*X+n̂[2]*Y+n̂[3]*Z)/2)
        qureg = createQureg(numQubits, env)
        v = [1.0, 0.0]
        v = U*v
        pauliX(qureg, control)
        
        controlledRotateAroundAxis(qureg, control, target, θ, qVector(n̂))

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        one_ind=2^(control-1)+2^(target-1)+1
        zero_ind=2^(control-1)+1
        for ind =1:2^numQubits
            if ind == one_ind
                @test state_vec[ind] ≈ v[2] atol = tolerance
            elseif ind == zero_ind
                @test state_vec[ind] ≈ v[1] atol = tolerance
            else
                @test state_vec[ind] ≈ Complex(0,0) atol = tolerance
            end
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end


function test_controlledRotateX()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(2:12)
        target = rand(1:numQubits)
        control = rand(filter(x->x≠target, 1:numQubits))
        θ = rand(Float64)
        
        U = ℯ^(-im*θ*X/2)
        qureg = createQureg(numQubits, env)
        v = [1.0, 0.0]
        v = U*v
        pauliX(qureg, control)
        
        controlledRotateX(qureg, control, target, θ)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        one_ind=2^(control-1)+2^(target-1)+1
        zero_ind=2^(control-1)+1
        for ind =1:2^numQubits
            if ind == one_ind
                @test state_vec[ind] ≈ v[2] atol = tolerance
            elseif ind == zero_ind
                @test state_vec[ind] ≈ v[1] atol = tolerance
            else
                @test state_vec[ind] ≈ Complex(0,0) atol = tolerance
            end
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_controlledRotateY()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(2:12)
        target = rand(1:numQubits)
        control = rand(filter(x->x≠target, 1:numQubits))
        θ = rand(Float64)
        
        U = ℯ^(-im*θ*Y/2)
        qureg = createQureg(numQubits, env)
        v = [1.0, 0.0]
        v = U*v
        pauliX(qureg, control)
        
        controlledRotateY(qureg, control, target, θ)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        one_ind=2^(control-1)+2^(target-1)+1
        zero_ind=2^(control-1)+1
        for ind =1:2^numQubits
            if ind == one_ind
                @test state_vec[ind] ≈ v[2] atol = tolerance
            elseif ind == zero_ind
                @test state_vec[ind] ≈ v[1] atol = tolerance
            else
                @test state_vec[ind] ≈ Complex(0,0) atol = tolerance
            end
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end



function test_controlledRotateZ()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(2:12)
        target = rand(1:numQubits)
        control = rand(filter(x->x≠target, 1:numQubits))
        θ = rand(Float64)
        
        U = ℯ^(-im*θ*Z/2)
        qureg = createQureg(numQubits, env)
        v = [1.0, 0.0]
        v = U*v
        pauliX(qureg, control)
        
        controlledRotateZ(qureg, control, target, θ)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        one_ind=2^(control-1)+2^(target-1)+1
        zero_ind=2^(control-1)+1
        for ind =1:2^numQubits
            if ind == one_ind
                @test state_vec[ind] ≈ v[2] atol = tolerance
            elseif ind == zero_ind
                @test state_vec[ind] ≈ v[1] atol = tolerance
            else
                @test state_vec[ind] ≈ Complex(0,0) atol = tolerance
            end
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end


function test_controlledTwoQubitUnitary()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(3:12)
    
        qureg = createQureg(numQubits, env)
        v = zeros(Float64, 2^numQubits)
        v[1]=1.0
        M_j = rand(Haar(2), 4)
        res = numQubits-3
        Idm = 1.0*Matrix(I,2^res, 2^res)
        U=kron(Idm, M_j, X)
        v=U*v
        
        pauliX(qureg, 1)
        controlledTwoQubitUnitary(qureg, 1, 2, 3, M_j)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        
        for ind =1:2^numQubits
            @test state_vec[ind] ≈ v[ind] atol = tolerance
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_controlledUnitary()
    env= createQuESTEnv()
    for t=1:10
        numQubits = rand(2:12)
        qureg = createQureg(numQubits, env)
        target = rand(1:numQubits)
        control = rand(filter(x->x≠target, 1:numQubits))
        M_j = rand(Haar(2), 2)
        v = [1.0, 0.0]
        v = M_j*v
        pauliX(qureg, control)
        controlledUnitary(qureg, control, target, M_j)
        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        one_ind=2^(control-1)+2^(target-1)+1
        zero_ind=2^(control-1)+1
        for ind =1:2^numQubits
            if ind == one_ind
                @test state_vec[ind] ≈ v[2] atol = tolerance
            elseif ind == zero_ind
                @test state_vec[ind] ≈ v[1] atol = tolerance
            else
                @test state_vec[ind] ≈ Complex(0,0) atol = tolerance
            end
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end


function test_hadamard()
    env= createQuESTEnv()
    for t=1:10
        numQubits = rand(1:12)
        qureg = createQureg(numQubits, env)
        target = rand(1:numQubits)
        state_vec = create_state_vec_zero_state(numQubits)
        state_vec=apply_hadamard_to_state_vector_at_target_index(state_vec, target)
        hadamard(qureg, target)
        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        for ind =1:2^numQubits
            @test state_vec[ind] ≈ state_vec[ind] atol = tolerance
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_multiControlledMultiQubitUnitary()
   
    env= createQuESTEnv()
    for t=1:10
        numQubits = rand(4:12)
        num_ctrls = rand(2:numQubits-2)
        num_targs = rand(2:numQubits-num_ctrls)

        ctrls = [x for x in 1:num_ctrls]
        targs = [x for x in num_ctrls+1:num_ctrls+num_targs]
        #targs = Cint.([x for x in num_ctrls:num_ctrls+num_targs-1])
        #print(ctrls)
        #print(targs)
        
        M = rand(Haar(2), 2^num_targs)
        
    
        qureg = createQureg(numQubits, env)
        
        for ind in ctrls
            pauliX(qureg, ind)
        end
        
        multiControlledMultiQubitUnitary(qureg, ctrls, targs, M)
        
        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)

        X=Matrix([0 1.0;1 0])
        big_X=Matrix([0 1.0;1 0])
        for iter = 1:num_ctrls-1
            big_X=kron(big_X, X)
        end
        res = numQubits-num_targs-num_ctrls
        Idm = 1.0*Matrix(I,2^res, 2^res)
        U=kron(Idm, M, big_X)
        v = zeros(Float64, 2^numQubits)
        v[1]=1.0
        v=U*v
        for ind = 1:2^numQubits
            @test state_vec[ind] ≈ v[ind] atol = tolerance
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
    
end

function test_multiControlledPhaseFlip()
    env= createQuESTEnv()
    for t=1:10
        
        numQubits = rand(2:12)
        num_ctrls = rand(1:numQubits)

        ctrls = [x for x in 1:num_ctrls]

        qureg = createQureg(numQubits, env)

        for ind = 1:num_ctrls
            pauliX(qureg, ind)
        end

        multiControlledPhaseFlip(qureg, ctrls)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        one_ind=2^(num_ctrls)
        for ind =1:2^numQubits
            if ind == one_ind
                @test state_vec[ind] ≈ Complex(-1,0) atol = tolerance
            else
                @test state_vec[ind] ≈ Complex(0,0) atol = tolerance
            end
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_multiControlledPhaseShift()
    env= createQuESTEnv()
    for t=1:10
        
        numQubits = rand(2:12)
        num_ctrls = rand(1:numQubits)

        ctrls = [x for x in 1:num_ctrls]
        θ = rand(Float64)
        qureg = createQureg(numQubits, env)

        for ind in ctrls
            pauliX(qureg, ind)
        end
        
        multiControlledPhaseShift(qureg, ctrls, θ)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        one_ind=2^(num_ctrls)
        for ind =1:2^numQubits
            if ind == one_ind
                @test state_vec[ind] ≈ℯ^(im*θ) atol = tolerance
            else
                @test state_vec[ind] ≈ Complex(0,0) atol = tolerance
            end
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end


function test_multiControlledTwoQubitUnitary()
    env= createQuESTEnv()
    X=pauliX_matrix()
   
    numQubits = 5
    qureg = createQureg(numQubits, env)
    num_ctrls = 3
    targ1,targ2 = 4,5

    ctrls = [1,2,3]
    v = create_state_vec_zero_state(numQubits)
    M = rand(Haar(2), 4)
    
    big_X = reduce(kron,[X for x in Base.OneTo(num_ctrls)])
    
    res = numQubits-2-num_ctrls
    Idm = 1.0*Matrix(I,2^res, 2^res)
    U=kron(M,big_X)
    v=U*v
    
    [pauliX(qureg, ind) for ind in ctrls]
    multiControlledTwoQubitUnitary(qureg, ctrls, targ1,targ2, M)
    state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
    @test state_vec ≈ v atol = tolerance
    destroyQureg(qureg, env)

    destroyQuESTEnv(env)
end












function test_sync()
    env= createQuESTEnv()
    numQubits = rand(3:12)
    qureg = createQureg(numQubits, env)
    state = rand(Complex{Float64}, 2^numQubits)
    state /= norm(state)
    initStateFromAmps(qureg, real.(state), imag.(state))

    syncQuESTEnv(env)
    @test syncQuESTSuccess(1) == 1

    destroyQureg(qureg, env)
    destroyQuESTEnv(env)

end