using Pkg
Pkg.activate(".")
include("src/QuEST.jl")
using .QuEST
using Test
using RandomMatrices: Haar
using LinearAlgebra
tolerance = 1e-10









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
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(2:12)
        qureg = createQureg(numQubits, env)
        target = rand(1:numQubits)
        control = rand(filter(x->x≠target, 1:numQubits))
        M_j = rand(Haar(2), 2)
        if Float64 == Float32
            M_j = Matrix{Complex{Float64}}(M_j)
        end
        v = [1.0, 0.0]
        v = M_j*v
        pauliX(qureg, control)
        
        controlledUnitary(qureg, control, target, M_j)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        one_ind=2^(control-1)+2^(target-1)+1
        zero_ind=2^(control-1)+1
        for ind =1:2^numQubits
            if ind == one_ind
                @test reals[ind] ≈ real(v[2]) atol = tolerance
                @test imags[ind] ≈ imag(v[2]) atol = tolerance
            elseif ind == zero_ind
                @test reals[ind] ≈ real(v[1]) atol = tolerance
                @test imags[ind] ≈ imag(v[1]) atol = tolerance
            else
                @test reals[ind] ≈ 0 atol = tolerance
                @test imags[ind] ≈ 0 atol = tolerance
            end
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_hadamard()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    H=Matrix([1.0 1.0;1.0 -1.0])/sqrt(2)
    for t=1:10
        
        numQubits = rand(1:12)
        qureg = createQureg(numQubits, env)
        target = rand(0:numQubits-1)

        v = zeros(Float64, 2^numQubits)
        v[1]=1
        i_right=1.0*Matrix(I, 2^target, 2^target)
        res = numQubits-target
        i_left =1.0*Matrix(I, 2^res, 2^res)
        U=kron(i_left, H, i_right)
        v=U*v

        hadamard(qureg, target)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)

        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_multiControlledMultiQubitUnitary()
    print("warning check test_multiControlledMultiQubitUnitary\n")
    env= createQuESTEnv()
    for t=1:10
        numQubits = rand(2:12)
        num_ctrls = rand(1:numQubits-1)
        num_targs = rand(1:numQubits-num_ctrls)

        ctrls = Cint.([x for x in 0:num_ctrls-1])
        targs = Cint.([x for x in num_targs+num_ctrls-1:-1:num_ctrls])
        targs = Cint.([x for x in num_ctrls:num_ctrls+num_targs-1])
        #print(ctrls)
        #print(targs)
        
        M_j = rand(Haar(2), 2^num_targs)
        if Float64 == Float32
            M_j = Matrix{Complex{Float64}}(M_j)
        end
        M = make_QuEST_matrix(M_j)
    
        qureg = createQureg(numQubits, env)
        
        for ind in ctrls
            pauliX(qureg, ind)
        end
        #pauliX(qureg, 0)
        multiControlledMultiQubitUnitary(qureg, ctrls, targs, M)
        
        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)

        X=Matrix([0 1.0;1 0])
        big_X=Matrix([0 1.0;1 0])
        for iter = 1:num_ctrls-1
            big_X=kron(big_X, X)
        end
        res = numQubits-num_targs-num_ctrls
        Idm = 1.0*Matrix(I,2^res, 2^res)
        U=kron(Idm, M_j, big_X)
        v = zeros(Float64, 2^numQubits)
        v[1]=1.0
        v=U*v
        #print(reals)
        #print(imags)
        #print(v)
        for ind = 1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end
        destroyComplexMatrixN(M)
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
    
end


function test_multiControlledPhaseFlip()
    env= createQuESTEnv()
    for t=1:10
        
        numQubits = rand(2:12)
        num_ctrls = rand(1:numQubits)

        ctrls = Cint.([x for x in 0:num_ctrls-1])

        qureg = createQureg(numQubits, env)

        for ind = 0:num_ctrls-1
            pauliX(qureg, ind)
        end

        multiControlledPhaseFlip(qureg, ctrls)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        one_ind=2^(num_ctrls)
        for ind =1:2^numQubits
            if ind == one_ind
                @test reals[ind] ≈ -1 atol = tolerance
                @test imags[ind] ≈ 0 atol = tolerance
            else
                @test reals[ind] ≈ 0 atol = tolerance
                @test imags[ind] ≈ 0 atol = tolerance
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

        ctrls = Cint.([x for x in 0:num_ctrls-1])
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
                @test reals[ind] ≈ real(ℯ^(im*θ)) atol = tolerance
                @test imags[ind] ≈ imag(ℯ^(im*θ)) atol = tolerance
            else
                @test reals[ind] ≈ 0 atol = tolerance
                @test imags[ind] ≈ 0 atol = tolerance
            end
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_multiControlledTwoQubitUnitary()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
    
        
        numQubits = rand(3:12)
        qureg = createQureg(numQubits, env)
        num_ctrls = rand(1:numQubits-2)

        ctrls = Cint.([x for x in 0:num_ctrls-1])
        v = zeros(Float64, 2^numQubits)
        v[1]=1.0
        M_j = rand(Haar(2), 4)
        if Float64 == Float32
            M_j = Matrix{Complex{Float64}}(M_j)
        end

        big_X=Matrix([0 1.0;1 0])
        for iter = 1:num_ctrls-1
            big_X=kron(big_X, X)
        end
        res = numQubits-2-num_ctrls
        Idm = 1.0*Matrix(I,2^res, 2^res)
        U=kron(Idm, M_j, big_X)

        U=kron(Idm, M_j, big_X)
        v=U*v
        
        for ind in ctrls
            pauliX(qureg, ind)
        end
        multiControlledTwoQubitUnitary(qureg, ctrls, num_ctrls, num_ctrls+1, M_j)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        
        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_multiControlledUnitary()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(3:12)
        qureg = createQureg(numQubits, env)
        num_ctrls = rand(1:numQubits-1)

        ctrls = Cint.([x for x in 0:num_ctrls-1])

        M_j = rand(Haar(2), 2)
        if Float64 == Float32
            M_j = Matrix{Complex{Float64}}(M_j)
        end
        
        v = zeros(Float64, 2^numQubits)
        v[1]=1.0

        big_X=Matrix([0 1.0;1 0])
        for iter = 1:num_ctrls-1
            big_X=kron(big_X, X)
        end
        res = numQubits-1-num_ctrls
        Idm = 1.0*Matrix(I,2^res, 2^res)
        U=kron(Idm, M_j, big_X)

        U=kron(Idm, M_j, big_X)
        v=U*v
        
        for ind in ctrls
            pauliX(qureg, ind)
        end
        
        multiControlledUnitary(qureg, ctrls, num_ctrls, M_j)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        

        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end


        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_multiQubitUnitary()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(3:12)
        qureg = createQureg(numQubits, env)
        num_targs = rand(1:numQubits-1)

        targs = Cint.([x for x in 0:num_targs-1])

        M_j = rand(Haar(2), 2^num_targs)
        if Float64 == Float32
            M_j = Matrix{Complex{Float64}}(M_j)
        end
        M = make_QuEST_matrix(M_j)

        v = zeros(Float64, 2^numQubits)
        v[1]=1.0

        res = numQubits-num_targs
        Idm = 1.0*Matrix(I,2^res, 2^res)
        U=kron(Idm, M_j)
        v=U*v
        
        
        multiQubitUnitary(qureg, targs, M)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        

        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end


        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_multiRotatePauli()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    Id=Matrix([1.0 0;0 1])
    paulies = [QuEST_Types.PAULI_I, 
               QuEST_Types.PAULI_X, 
               QuEST_Types.PAULI_Y, 
               QuEST_Types.PAULI_Z]

    pauli_dict=Dict()
    pauli_dict[QuEST_Types.PAULI_I]=Id
    pauli_dict[QuEST_Types.PAULI_X]=X
    pauli_dict[QuEST_Types.PAULI_Y]=Y
    pauli_dict[QuEST_Types.PAULI_Z]=Z
    for t=1:10
        
        numQubits = rand(3:12)
        qureg = createQureg(numQubits, env)
        
        codes = rand(paulies, numQubits)

        c_codes = Array{QuEST_Types.pauliOpType, 1}()
        c_targs = Array{Cint, 1}()
        for i =1:numQubits
            if codes[i] != QuEST_Types.PAULI_I
                push!(c_codes, codes[i])
                push!(c_targs, i-1)
            end
        end

        
        M = ones(1,1)
        for ind =1:numQubits
            M=kron(pauli_dict[codes[ind]], M)
        end

        v = zeros(Float64, 2^numQubits)
        v[1]=1.0

        θ = rand(Float64)
        U=ℯ^(-im*θ*M/2)
        v=U*v
        
        
        multiRotatePauli(qureg, c_targs, c_codes, θ)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        

        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end


        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_multiRotateZ()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    Id=Matrix([1.0 0;0 1])
    
    for t=1:10
        
        numQubits = rand(3:12)
        qureg = createQureg(numQubits, env)
        
        codes = rand(0:1, numQubits)
        θ = rand(Float64)
        
        c_targs = Array{Cint, 1}()
        push!(c_targs, 0)
        for ind =2:numQubits
            if codes[ind] == 1 
                push!(c_targs, ind-1)
            end
        end

        
        M = ones(1,1)
        for ind =1:numQubits
            if codes[ind] == 1
                M=kron(Z, M)
            else
                M=kron(Id, M)
            end
        end

        v = zeros(Float64, 2^numQubits)
        v[1]=1.0

        
        U=ℯ^(-im*θ*M/2)
        v=U*v
        
        
        multiRotateZ(qureg, c_targs, θ)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        

        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end


        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end


function test_multiStateControlledUnitary()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    Id=Matrix([1.0 0;0 1])
    for t=1:10
        
        numQubits = rand(3:12)
        qureg = createQureg(numQubits, env)
        
        codes = rand(0:1, numQubits)
        target = rand(1:numQubits) 
        codes[target] = 3

        c_qubits = Array{Cint, 1}()
        c_state  = Array{Cint, 1}()

        for ind =1:numQubits
            if codes[ind] in [0,1]
                push!(c_qubits, ind-1)
                push!(c_state, codes[ind])
            end
        end

        M_j = rand(Haar(2), 2)
        if Float64 == Float32
            M_j = Matrix{Complex{Float64}}(M_j)
        end
        
        v = zeros(Float64, 2^numQubits)
        v[1]=1.0

        U = ones(1, 1)
        for ind =1:numQubits
            if codes[ind] == 0
                U = kron(Id, U)
            end

            if codes[ind] == 2
                U = kron(Id, U)
            end

            if codes[ind] == 1
                U = kron(X, U)
            end

            if codes[ind] == 3
                U = kron(M_j, U)
            end
        end
        v=U*v
        
        for ind =1:numQubits
            if codes[ind]==1
                pauliX(qureg, ind-1)
            end
        end
        multiStateControlledUnitary(qureg, c_qubits, c_state, target, M_j)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        

        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end


        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_pauliX()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    H=Matrix([1.0 1.0;1.0 -1.0])/sqrt(2)
    for t=1:10
        
        numQubits = rand(1:12)
        qureg = createQureg(numQubits, env)
        target = rand(0:numQubits-1)

        v = zeros(Float64, 2^numQubits)
        v[1]=1
        i_right=1.0*Matrix(I, 2^target, 2^target)
        res = numQubits-target
        i_left =1.0*Matrix(I, 2^res, 2^res)
        U=kron(i_left, X, i_right)
        v=U*v

        pauliX(qureg, target)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)

        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_pauliY()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    H=Matrix([1.0 1.0;1.0 -1.0])/sqrt(2)
    for t=1:10
        
        numQubits = rand(1:12)
        qureg = createQureg(numQubits, env)
        target = rand(0:numQubits-1)

        v = zeros(Float64, 2^numQubits)
        v[1]=1
        i_right=1.0*Matrix(I, 2^target, 2^target)
        res = numQubits-target
        i_left =1.0*Matrix(I, 2^res, 2^res)
        U=kron(i_left, Y, i_right)
        v=U*v

        pauliY(qureg, target)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)

        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_pauliZ()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    H=Matrix([1.0 1.0;1.0 -1.0])/sqrt(2)
    for t=1:10
        
        numQubits = rand(1:12)
        qureg = createQureg(numQubits, env)
        target = rand(0:numQubits-1)

        v = zeros(Float64, 2^numQubits)
        v[1]=1
        i_right=1.0*Matrix(I, 2^target, 2^target)
        res = numQubits-target
        i_left =1.0*Matrix(I, 2^res, 2^res)
        U=kron(i_left, Z, i_right)
        v=U*v

        pauliZ(qureg, target)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)

        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_phaseShift()
    env= createQuESTEnv()
    for t=1:10
        
        numQubits = rand(2:12)
        target = rand(0:numQubits-1)
        
        θ = rand(Float64)
        qureg = createQureg(numQubits, env)

        
        pauliX(qureg, target)
        phaseShift(qureg, target, θ)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        one_ind=2^target+1
        for ind =1:2^numQubits
            if ind == one_ind
                @test reals[ind] ≈ real(ℯ^(im*θ)) atol = tolerance
                @test imags[ind] ≈ imag(ℯ^(im*θ)) atol = tolerance
            else
                @test reals[ind] ≈ 0 atol = tolerance
                @test imags[ind] ≈ 0 atol = tolerance
            end
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_rotateAroundAxis()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(2:12)
        target = rand(0:numQubits-1)
        
        θ = rand(Float64)
        n̂ = rand(Float64, 3)
        n̂ /= norm(n̂)
        U = ℯ^(-im*θ*(n̂[1]*X+n̂[2]*Y+n̂[3]*Z)/2)
        qureg = createQureg(numQubits, env)

        v = [1.0, 0.0]
        v = U*v
        
        rotateAroundAxis(qureg, target, θ, n̂)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        one_ind=2^target+1
        zero_ind=1
        for ind =1:2^numQubits
            if ind == one_ind
                @test reals[ind] ≈ real(v[2]) atol = tolerance
                @test imags[ind] ≈ imag(v[2]) atol = tolerance
            elseif ind == zero_ind
                @test reals[ind] ≈ real(v[1]) atol = tolerance
                @test imags[ind] ≈ imag(v[1]) atol = tolerance
            else
                @test reals[ind] ≈ 0 atol = tolerance
                @test imags[ind] ≈ 0 atol = tolerance
            end
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_rotateX()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(2:12)
        target = rand(0:numQubits-1)
        control = rand(1:numQubits)

        θ = rand(Float64)
        
        U = ℯ^(-im*θ*X/2)
        qureg = createQureg(numQubits, env)

        v = [1.0, 0.0]
        v = U*v
        
        
        rotateX(qureg, target, θ)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        one_ind=2^target+1
        zero_ind=+1
        for ind =1:2^numQubits
            if ind == one_ind
                @test reals[ind] ≈ real(v[2]) atol = tolerance
                @test imags[ind] ≈ imag(v[2]) atol = tolerance
            elseif ind == zero_ind
                @test reals[ind] ≈ real(v[1]) atol = tolerance
                @test imags[ind] ≈ imag(v[1]) atol = tolerance
            else
                @test reals[ind] ≈ 0 atol = tolerance
                @test imags[ind] ≈ 0 atol = tolerance
            end
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_rotateY()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(2:12)
        target = rand(0:numQubits-1)
        control = rand(1:numQubits)

        θ = rand(Float64)
        
        U = ℯ^(-im*θ*Y/2)
        qureg = createQureg(numQubits, env)

        v = [1.0, 0.0]
        v = U*v
        
        
        rotateY(qureg, target, θ)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        one_ind=2^target+1
        zero_ind=+1
        for ind =1:2^numQubits
            if ind == one_ind
                @test reals[ind] ≈ real(v[2]) atol = tolerance
                @test imags[ind] ≈ imag(v[2]) atol = tolerance
            elseif ind == zero_ind
                @test reals[ind] ≈ real(v[1]) atol = tolerance
                @test imags[ind] ≈ imag(v[1]) atol = tolerance
            else
                @test reals[ind] ≈ 0 atol = tolerance
                @test imags[ind] ≈ 0 atol = tolerance
            end
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_rotateZ()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(2:12)
        target = rand(0:numQubits-1)
        control = rand(1:numQubits)

        θ = rand(Float64)
        
        U = ℯ^(-im*θ*Z/2)
        qureg = createQureg(numQubits, env)

        v = [1.0, 0.0]
        v = U*v
        
        
        rotateZ(qureg, target, θ)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        one_ind=2^target+1
        zero_ind=+1
        for ind =1:2^numQubits
            if ind == one_ind
                @test reals[ind] ≈ real(v[2]) atol = tolerance
                @test imags[ind] ≈ imag(v[2]) atol = tolerance
            elseif ind == zero_ind
                @test reals[ind] ≈ real(v[1]) atol = tolerance
                @test imags[ind] ≈ imag(v[1]) atol = tolerance
            else
                @test reals[ind] ≈ 0 atol = tolerance
                @test imags[ind] ≈ 0 atol = tolerance
            end
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_sGate()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    H=Matrix([1.0 1.0;1.0 -1.0])/sqrt(2)
    S=Matrix([1.0 0;0 im])
    for t=1:10
        
        numQubits = rand(1:12)
        qureg = createQureg(numQubits, env)
        target = rand(0:numQubits-1)

        v = zeros(Float64, 2^numQubits)
        v[1]=1
        i_right=1.0*Matrix(I, 2^target, 2^target)
        res = numQubits-target
        i_left =1.0*Matrix(I, 2^res, 2^res)
        U=kron(i_left, S, i_right)
        v=U*v

        sGate(qureg, target)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)

        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_sqrtSwapGate()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    H=Matrix([1.0 1.0;1.0 -1.0])/sqrt(2)
    S=Matrix([1.0 0;0 im])
    Swap=Matrix([1 0 0 0;0 0 1 0;0 1 0 0;0 0 0 1])
    SqrtSwap=sqrt(Swap)
    for t=1:10
        
        numQubits = rand(3:12)
        qureg = createQureg(numQubits, env)
        qubit1 = rand(0:numQubits-2)

        v = zeros(Float64, 2^numQubits)
        v[1]=1
        i_right=1.0*Matrix(I, 2^qubit1, 2^qubit1)
        res = numQubits-qubit1-2
        i_left =1.0*Matrix(I, 2^res, 2^res)
        U=kron(i_left, SqrtSwap, i_right)
        v=U*v

        sqrtSwapGate(qureg, qubit1, qubit1+1)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)

        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_swapGate()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    H=Matrix([1.0 1.0;1.0 -1.0])/sqrt(2)
    S=Matrix([1.0 0;0 im])
    Swap=Matrix([1 0 0 0;0 0 1 0;0 1 0 0;0 0 0 1])
    
    for t=1:10
        
        numQubits = rand(3:12)
        qureg = createQureg(numQubits, env)
        qubit1 = rand(0:numQubits-2)

        v = zeros(Float64, 2^numQubits)
        v[1]=1
        i_right=1.0*Matrix(I, 2^qubit1, 2^qubit1)
        res = numQubits-qubit1-2
        i_left =1.0*Matrix(I, 2^res, 2^res)
        U=kron(i_left, Swap, i_right)
        v=U*v

        swapGate(qureg, qubit1, qubit1+1)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)

        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_tGate()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    H=Matrix([1.0 1.0;1.0 -1.0])/sqrt(2)
    T=Matrix([1.0 0;0 ℯ^(im*π/4)])
    for t=1:10
        
        numQubits = rand(1:12)
        qureg = createQureg(numQubits, env)
        target = rand(0:numQubits-1)

        v = zeros(Float64, 2^numQubits)
        v[1]=1
        i_right=1.0*Matrix(I, 2^target, 2^target)
        res = numQubits-target
        i_left =1.0*Matrix(I, 2^res, 2^res)
        U=kron(i_left, T, i_right)
        v=U*v

        tGate(qureg, target)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)

        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_twoQubitUnitary()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        
        M_j = rand(Haar(2), 4)
        if Float64 == Float32
            M_j = Matrix{Complex{Float64}}(M_j)
        end

        numQubits = rand(3:12)
        qureg = createQureg(numQubits, env)
        qubit1 = rand(0:numQubits-2)

        v = zeros(Float64, 2^numQubits)
        v[1]=1.0
        i_right=1.0*Matrix(I, 2^qubit1, 2^qubit1)
        res = numQubits-qubit1-2
        i_left =1.0*Matrix(I, 2^res, 2^res)
        U=kron(i_left, M_j, i_right)
        v=U*v
        
        
        twoQubitUnitary(qureg, qubit1, qubit1+1, M_j)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        
        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end


function test_unitary()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(2:12)
        qureg = createQureg(numQubits, env)
        target = rand(0:numQubits-1)
        
        M_j = rand(Haar(2), 2)
        if Float64 == Float32
            M_j = Matrix{Complex{Float64}}(M_j)
        end
        v = [1.0, 0.0]
        v = M_j*v
        
        
        unitary(qureg, target, M_j)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        one_ind=2^target+1
        zero_ind=1
        for ind =1:2^numQubits
            if ind == one_ind
                @test reals[ind] ≈ real(v[2]) atol = tolerance
                @test imags[ind] ≈ imag(v[2]) atol = tolerance
            elseif ind == zero_ind
                @test reals[ind] ≈ real(v[1]) atol = tolerance
                @test imags[ind] ≈ imag(v[1]) atol = tolerance
            else
                @test reals[ind] ≈ 0 atol = tolerance
                @test imags[ind] ≈ 0 atol = tolerance
            end
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end



function test_collapseToOutcome()
    env= createQuESTEnv()
    for t=1:100
        numQubits = rand(2:12)
        θ=rand(Float64)*π
        qureg = createQureg(numQubits, env)
        
        target=rand(0:numQubits-1)

        for ind =0:numQubits-1
            if ind!=target
                rotateX(qureg, ind, rand(Float64)*π)
            end
        end

        rotateX(qureg, target, θ)

        p = collapseToOutcome(qureg, target, 0)

        @test p ≈ (1 + cos(θ))/2 atol = tolerance

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)


    end
    
    destroyQuESTEnv(env)

end

function test_measure()
    env= createQuESTEnv()
    for t=1:100
        numQubits = rand(2:12)
        θ=rand(Float64)*π
        qureg = createQureg(numQubits, env)
        
        target=rand(0:numQubits-1)

        for ind =0:numQubits-1
            if ind!=target
                rotateX(qureg, ind, rand(Float64)*π)
            end
        end

        rotateX(qureg, target, θ)

        p = measure(qureg, target)

        @test p in [0, 1]


    end
    
    destroyQuESTEnv(env)

end

function test_measureWithStats()
    env= createQuESTEnv()
    for t=1:100
        numQubits = 1
        θ=rand(Float64)*π
        qureg = createQureg(numQubits, env)        

        rotateX(qureg, 0, θ)

        cbit, prob = measureWithStats(qureg, 0)

        @test cbit in [0, 1]
        if cbit == 0 
            @test prob ≈ (1 + cos(θ))/2 atol = tolerance
        else
            @test prob ≈ (1 - cos(θ))/2 atol = tolerance
        end

    end
    
    destroyQuESTEnv(env)

end

function test_DiagonalOp()
    env= createQuESTEnv()
    for i=1:100
        num_qubits = rand(1:12)
        qureg = createQureg(num_qubits, env)

        op = createDiagonalOp(num_qubits, env)
        op_j = rand(Complex{Float64}, 2^num_qubits)
        initDiagonalOp(op, real(op_j), imag(op_j))

        codes=rand(0:1, num_qubits)
        #codes=zeros(num_qubits)
        basis = 1
        for ind=0:num_qubits-1
            if codes[ind+1] == 1
                pauliX(qureg, ind)
                basis += 2^ind
            end
        end
        applyDiagonalOp(qureg, op)
        reals = unsafe_wrap(Vector{Float64}, qureg.stateVec.real, 2^num_qubits)
        imags = unsafe_wrap(Vector{Float64}, qureg.stateVec.imag, 2^num_qubits)
    
        @test real(op_j[basis]) ≈ reals[basis] atol = tolerance
        @test imag(op_j[basis]) ≈ imags[basis] atol = tolerance

        destroyDiagonalOp(op, env)
    end
    destroyQuESTEnv(env)
end

function test_applyMatrix2()
    env= createQuESTEnv()
    for t=1:10
        
        numQubits = rand(2:12)
        qureg = createQureg(numQubits, env)
        target = rand(0:numQubits-1)
        
        M_j = rand(Haar(2), 2)
        if Float64 == Float32
            M_j = Matrix{Complex{Float64}}(M_j)
        end
        v = [1.0, 0.0]
        v = M_j*v
        
        M = _quest_mtx_2(M_j)
        applyMatrix2(qureg, target, M)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        one_ind=2^target+1
        zero_ind=1
        for ind =1:2^numQubits
            if ind == one_ind
                @test reals[ind] ≈ real(v[2]) atol = tolerance
                @test imags[ind] ≈ imag(v[2]) atol = tolerance
            elseif ind == zero_ind
                @test reals[ind] ≈ real(v[1]) atol = tolerance
                @test imags[ind] ≈ imag(v[1]) atol = tolerance
            else
                @test reals[ind] ≈ 0 atol = tolerance
                @test imags[ind] ≈ 0 atol = tolerance
            end
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_applyMatrix4()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        
        M_j = rand(Haar(2), 4)
        if Float64 == Float32
            M_j = Matrix{Complex{Float64}}(M_j)
        end

        M = _quest_mtx_4(M_j)

        numQubits = rand(3:12)
        qureg = createQureg(numQubits, env)
        qubit1 = rand(0:numQubits-2)

        v = zeros(Float64, 2^numQubits)
        v[1]=1.0
        i_right=1.0*Matrix(I, 2^qubit1, 2^qubit1)
        res = numQubits-qubit1-2
        i_left =1.0*Matrix(I, 2^res, 2^res)
        U=kron(i_left, M_j, i_right)
        v=U*v
        
        
        applyMatrix4(qureg, qubit1, qubit1+1, M)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        
        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_applyMatrixN()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(3:12)
        qureg = createQureg(numQubits, env)
        num_targs = rand(1:numQubits-1)

        targs = Cint.([x for x in 0:num_targs-1])

        M_j = rand(Haar(2), 2^num_targs)
        if Float64 == Float32
            M_j = Matrix{Complex{Float64}}(M_j)
        end
        M = make_QuEST_matrix(M_j)

        v = zeros(Float64, 2^numQubits)
        v[1]=1.0

        res = numQubits-num_targs
        Idm = 1.0*Matrix(I,2^res, 2^res)
        U=kron(Idm, M_j)
        v=U*v
        
        
        applyMatrixN(qureg, targs, M)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        

        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end


        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end

function test_applyMultiControlledMatrixN()
    print("warning check test_multiControlledMultiQubitUnitary\n")
    env= createQuESTEnv()
    for t=1:10
        numQubits = rand(2:12)
        num_ctrls = rand(1:numQubits-1)
        num_targs = rand(1:numQubits-num_ctrls)

        ctrls = Cint.([x for x in 0:num_ctrls-1])
        targs = Cint.([x for x in num_targs+num_ctrls-1:-1:num_ctrls])
        targs = Cint.([x for x in num_ctrls:num_ctrls+num_targs-1])
        #print(ctrls)
        #print(targs)
        
        M_j = rand(Haar(2), 2^num_targs)
        if Float64 == Float32
            M_j = Matrix{Complex{Float64}}(M_j)
        end
        M = make_QuEST_matrix(M_j)
    
        qureg = createQureg(numQubits, env)
        
        for ind in ctrls
            pauliX(qureg, ind)
        end
        #pauliX(qureg, 0)
        applyMultiControlledMatrixN(qureg, ctrls, targs, M)
        
        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)

        X=Matrix([0 1.0;1 0])
        big_X=Matrix([0 1.0;1 0])
        for iter = 1:num_ctrls-1
            big_X=kron(big_X, X)
        end
        res = numQubits-num_targs-num_ctrls
        Idm = 1.0*Matrix(I,2^res, 2^res)
        U=kron(Idm, M_j, big_X)
        v = zeros(Float64, 2^numQubits)
        v[1]=1.0
        v=U*v
        #print(reals)
        #print(imags)
        #print(v)
        for ind = 1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end
        destroyComplexMatrixN(M)
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
    
end


function test_applyPauliHamil()
    env= createQuESTEnv()

    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    Id=Matrix([1.0 0;0 1])

    paulies = [QuEST_Types.PAULI_I, 
               QuEST_Types.PAULI_X, 
               QuEST_Types.PAULI_Y, 
               QuEST_Types.PAULI_Z]

    pauli_dict=Dict()
    pauli_dict[QuEST_Types.PAULI_I]=Id
    pauli_dict[QuEST_Types.PAULI_X]=X
    pauli_dict[QuEST_Types.PAULI_Y]=Y
    pauli_dict[QuEST_Types.PAULI_Z]=Z
    
    for i=1:10
        numQubits = rand(1:12)
        numSumTerms = rand(1:20)
        hamil = createPauliHamil(numQubits, numSumTerms)


        codes = rand(paulies, numQubits*numSumTerms)
        coeffs = rand(Float64, numSumTerms)
        initPauliHamil(hamil, coeffs, codes)
        
        H=zeros(2^numQubits, 2^numQubits)
        for term=0:numSumTerms-1
            h_tmp=ones(1, 1)
            for qubit =0:numQubits-1
                h_tmp = kron(pauli_dict[codes[term*numQubits+qubit+1]], h_tmp)
            end
            H += coeffs[term+1]*h_tmp
        end

        v = zeros(Float64, 2^numQubits)
        v[1]=1.0
        v=H*v

        in_qureg = createQureg(numQubits, env)
        out_qureg = createQureg(numQubits, env)

        applyPauliHamil(in_qureg, hamil, out_qureg)

        reals = unsafe_wrap(Vector{Float64}, out_qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{Float64}, out_qureg.stateVec.imag, 2^numQubits)

        for ind = 1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end

        destroyQureg(in_qureg, env)
        destroyQureg(out_qureg, env)
        destroyPauliHamil(hamil)


    end
    destroyQuESTEnv(env)
end

function test_applyPauliSum()
    env= createQuESTEnv()

    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    Id=Matrix([1.0 0;0 1])

    paulies = [QuEST_Types.PAULI_I, 
               QuEST_Types.PAULI_X, 
               QuEST_Types.PAULI_Y, 
               QuEST_Types.PAULI_Z]

    pauli_dict=Dict()
    pauli_dict[QuEST_Types.PAULI_I]=Id
    pauli_dict[QuEST_Types.PAULI_X]=X
    pauli_dict[QuEST_Types.PAULI_Y]=Y
    pauli_dict[QuEST_Types.PAULI_Z]=Z
    
    for i=1:10
        numQubits = rand(1:12)
        numSumTerms = rand(1:20)


        codes = rand(paulies, numQubits*numSumTerms)
        coeffs = rand(Float64, numSumTerms)
        
        H=zeros(2^numQubits, 2^numQubits)
        for term=0:numSumTerms-1
            h_tmp=ones(1, 1)
            for qubit =0:numQubits-1
                h_tmp = kron(pauli_dict[codes[term*numQubits+qubit+1]], h_tmp)
            end
            H += coeffs[term+1]*h_tmp
        end

        v = zeros(Float64, 2^numQubits)
        v[1]=1.0
        v=H*v

        in_qureg = createQureg(numQubits, env)
        out_qureg = createQureg(numQubits, env)

        applyPauliSum(in_qureg, codes, coeffs, numSumTerms, out_qureg)

        reals = unsafe_wrap(Vector{Float64}, out_qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{Float64}, out_qureg.stateVec.imag, 2^numQubits)

        for ind = 1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end

        destroyQureg(in_qureg, env)
        destroyQureg(out_qureg, env)


    end
    destroyQuESTEnv(env)
end

function test_applyTrotterCircuit()
    env= createQuESTEnv()

    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    Id=Matrix([1.0 0;0 1])

    paulies = [QuEST_Types.PAULI_I, 
               QuEST_Types.PAULI_X, 
               #QuEST_Types.PAULI_Y, 
               QuEST_Types.PAULI_Z]

    pauli_dict=Dict()
    pauli_dict[QuEST_Types.PAULI_I]=Id
    pauli_dict[QuEST_Types.PAULI_X]=X
    pauli_dict[QuEST_Types.PAULI_Y]=Y
    pauli_dict[QuEST_Types.PAULI_Z]=Z
    
    for i=1:10
        numQubits = rand(1:12)
        numSumTerms = rand(1:4)
        hamil = createPauliHamil(numQubits, numSumTerms)


        codes = rand(paulies, numQubits*numSumTerms)
        coeffs = rand(Float64, numSumTerms)
        reps = rand(1:10)
        t=rand(Float64)

        initPauliHamil(hamil, coeffs, codes)
        qureg = createQureg(numQubits, env)
        test_qureg = createQureg(numQubits, env)
        
        #startRecordingQASM(test_qureg)
        
        for trot=1:reps
            for term=0:numSumTerms-1
                z_list=Vector{Cint}()
                for qubit =0:numQubits-1
                    if codes[term*numQubits+qubit+1] == QuEST_Types.PAULI_X
                        hadamard(test_qureg, qubit)
                    end
                end
                for qubit =0:numQubits-1
                    if codes[term*numQubits+qubit+1] != QuEST_Types.PAULI_I
                        push!(z_list, qubit)
                    end
                end

                for ind=1:length(z_list)-1
                    controlledNot(test_qureg, z_list[ind], z_list[ind+1])
                end

                if length(z_list) !=0
                    rotateZ(test_qureg, z_list[length(z_list)], 2*t*coeffs[term+1]/reps)
                end

                for ind=length(z_list)-1:-1:1
                    controlledNot(test_qureg, z_list[ind], z_list[ind+1])
                end

                for qubit =0:numQubits-1
                    if codes[term*numQubits+qubit+1] == QuEST_Types.PAULI_X
                        hadamard(test_qureg, qubit)
                    end
                end
            end

        end
            
        #stopRecordingQASM(test_qureg)
        #printRecordedQASM(test_qureg)
        #print("##################################################")
        #startRecordingQASM(qureg)
        applyTrotterCircuit(qureg, hamil, t, 1, reps)
        #stopRecordingQASM(qureg)
        #printRecordedQASM(qureg)

        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)

        reals2 = unsafe_wrap(Vector{Float64}, test_qureg.stateVec.real, 2^numQubits)
        imags2 = unsafe_wrap(Vector{Float64}, test_qureg.stateVec.imag, 2^numQubits)

        # for ind = 1:2^numQubits
        #    @test isapprox(reals[ind], real(v[ind]);atol = 1e-2)
        #    @test isapprox(imags[ind], imag(v[ind]); atol = 1e-2)
        # end

        for ind = 1:2^numQubits
            @test reals[ind] ≈ reals2[ind] atol = tolerance
            @test imags[ind] ≈ imags2[ind] atol = tolerance
        end

        destroyQureg(qureg, env)
        destroyPauliHamil(hamil)


    end
    destroyQuESTEnv(env)
end


function test_calcDensityInnerProduct()
    env= createQuESTEnv()
    
    for i=1:20
        numQubits = rand(1:5)
        
        targs = Cint.([x for x in 0:numQubits-1])

        M1 = rand(Haar(2), 2^numQubits)
        if Float64 == Float32
            M1 = Matrix{Complex{Float64}}(M1)
        end
        U1 = make_QuEST_matrix(M1)

        M2 = rand(Haar(2), 2^numQubits)
        if Float64 == Float32
            M2 = Matrix{Complex{Float64}}(M2)
        end
        U2 = make_QuEST_matrix(M2)

        v = zeros(Float64, 2^numQubits)
        v[1]=1.0

        state1 = M1*v
        state2 = M2*v

        my_inner_prod = transpose(state1)*conj(state2)
        my_inner_prod = my_inner_prod*conj(my_inner_prod)

        qureg1 = createDensityQureg(numQubits, env)
        qureg2 = createDensityQureg(numQubits, env)
        
        multiQubitUnitary(qureg1, targs, U1)
        multiQubitUnitary(qureg2, targs, U2)

        @test my_inner_prod ≈ calcDensityInnerProduct(qureg1, qureg2) atol = tolerance

        destroyQureg(qureg1, env)
        destroyQureg(qureg2, env)

    end
    destroyQuESTEnv(env)
end


function test_calcExpecDiagonalOp()
    env= createQuESTEnv()
    
    for i=1:20
        numQubits = rand(1:12)
        
        qureg = createQureg(numQubits, env)

        state = rand(Complex{Float64}, 2^numQubits)
        state /= norm(state)

        initStateFromAmps(qureg, real.(state), imag.(state))
        

        op = createDiagonalOp(numQubits, env)
        op_j = rand(Complex{Float64}, 2^numQubits)
        initDiagonalOp(op, real(op_j), imag(op_j))
        
        rhs = calcExpecDiagonalOp(qureg, op)

        lhs = sum([state[x]*conj(state[x])*op_j[x] for x =1:2^numQubits])

        @test rhs.real ≈ real(lhs) atol = tolerance
        @test rhs.imag ≈ imag(lhs) atol = tolerance

        destroyDiagonalOp(op, env)

        
        destroyQureg(qureg, env)

    end
    destroyQuESTEnv(env)
end

function test_calcExpecPauliHamil()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    Id=Matrix([1.0 0;0 1])

    paulies = [QuEST_Types.PAULI_I, 
               QuEST_Types.PAULI_X, 
               QuEST_Types.PAULI_Y, 
               QuEST_Types.PAULI_Z]

    pauli_dict=Dict()
    pauli_dict[QuEST_Types.PAULI_I]=Id
    pauli_dict[QuEST_Types.PAULI_X]=X
    pauli_dict[QuEST_Types.PAULI_Y]=Y
    pauli_dict[QuEST_Types.PAULI_Z]=Z
    
    for i=1:20
        numQubits = rand(1:12)
        numSumTerms = rand(1:20)

        hamil = createPauliHamil(numQubits, numSumTerms)
        codes = rand(paulies, numQubits*numSumTerms)
        coeffs = rand(Float64, numSumTerms)
        initPauliHamil(hamil, coeffs, codes)

        qureg = createQureg(numQubits, env)
        work = createQureg(numQubits, env)

        H=zeros(2^numQubits, 2^numQubits)
        for term=0:numSumTerms-1
            h_tmp=ones(1, 1)
            for qubit =0:numQubits-1
                h_tmp = kron(pauli_dict[codes[term*numQubits+qubit+1]], h_tmp)
            end
            H += coeffs[term+1]*h_tmp
        end

        state = rand(Complex{Float64}, 2^numQubits)
        state /= norm(state)

        initStateFromAmps(qureg, real.(state), imag.(state))
        

        lhs = transpose(conj(state))*H*state

        rhs = calcExpecPauliHamil(qureg, hamil, work)

        @test rhs ≈ real(lhs) atol = tolerance
        


        destroyQureg(qureg, env)
        destroyQureg(work, env)
    end
    destroyQuESTEnv(env)
end

function test_calcExpecPauliProd()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    Id=Matrix([1.0 0;0 1])

    paulies = [QuEST_Types.PAULI_I, 
               QuEST_Types.PAULI_X, 
               QuEST_Types.PAULI_Y, 
               QuEST_Types.PAULI_Z]

    pauli_dict=Dict()
    pauli_dict[QuEST_Types.PAULI_I]=Id
    pauli_dict[QuEST_Types.PAULI_X]=X
    pauli_dict[QuEST_Types.PAULI_Y]=Y
    pauli_dict[QuEST_Types.PAULI_Z]=Z
    
    for i=1:20
        numQubits = rand(1:12)

        codes = rand(paulies, numQubits)

        qureg = createQureg(numQubits, env)
        work = createQureg(numQubits, env)

        H=zeros(2^numQubits, 2^numQubits)
        
        H=ones(1, 1)
        c_targs = Vector{Cint}([x for x in 0:numQubits-1])
        #c_codes = Vector{QuEST_Types.pauliOpType}()

        H=zeros(2^numQubits, 2^numQubits)
        H=ones(1, 1)
        h_tmp=ones(1, 1)
        for qubit =0:numQubits-1
            H = kron(pauli_dict[codes[qubit+1]], H)
        end


        state = rand(Complex{Float64}, 2^numQubits)
        state /= norm(state)

        initStateFromAmps(qureg, real.(state), imag.(state))
        

        lhs = transpose(conj(state))*H*state

        rhs = calcExpecPauliProd(qureg, c_targs, codes, work)

        @test rhs ≈ real(lhs) atol = tolerance
        


        destroyQureg(qureg, env)
        destroyQureg(work, env)
    end
    destroyQuESTEnv(env)
end

function test_calcExpecPauliSum()
    env= createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    Id=Matrix([1.0 0;0 1])

    paulies = [QuEST_Types.PAULI_I, 
               QuEST_Types.PAULI_X, 
               QuEST_Types.PAULI_Y, 
               QuEST_Types.PAULI_Z]

    pauli_dict=Dict()
    pauli_dict[QuEST_Types.PAULI_I]=Id
    pauli_dict[QuEST_Types.PAULI_X]=X
    pauli_dict[QuEST_Types.PAULI_Y]=Y
    pauli_dict[QuEST_Types.PAULI_Z]=Z
    
    for i=1:20
        numQubits = rand(1:12)
        numSumTerms = rand(1:20)

        codes = rand(paulies, numQubits*numSumTerms)
        coeffs = rand(Float64, numSumTerms)

        qureg = createQureg(numQubits, env)
        work = createQureg(numQubits, env)

        H=zeros(2^numQubits, 2^numQubits)
        for term=0:numSumTerms-1
            h_tmp=ones(1, 1)
            for qubit =0:numQubits-1
                h_tmp = kron(pauli_dict[codes[term*numQubits+qubit+1]], h_tmp)
            end
            H += coeffs[term+1]*h_tmp
        end

        state = rand(Complex{Float64}, 2^numQubits)
        state /= norm(state)

        initStateFromAmps(qureg, real.(state), imag.(state))
        

        lhs = transpose(conj(state))*H*state

        rhs = calcExpecPauliSum(qureg, codes, coeffs, work)

        @test isapprox(rhs, real(lhs); atol=1e-5)
        


        destroyQureg(qureg, env)
        destroyQureg(work, env)
    end
    destroyQuESTEnv(env)
end

function test_calcFidelity()
    env= createQuESTEnv()
    
    for i=1:20
        numQubits = rand(1:12)


        qureg1 = createQureg(numQubits, env)
        qureg2 = createQureg(numQubits, env)

        
        state1 = rand(Complex{Float64}, 2^numQubits)
        state1 /= norm(state1)

        state2 = rand(Complex{Float64}, 2^numQubits)
        state2 /= norm(state2)

        initStateFromAmps(qureg1, real.(state1), imag.(state1))
        initStateFromAmps(qureg2, real.(state2), imag.(state2))

        lhs = abs(transpose(conj(state1))*state2)^2

        rhs = calcFidelity(qureg1, qureg2)

        @test isapprox(rhs, real(lhs); atol=1e-5)
        


        destroyQureg(qureg1, env)
        destroyQureg(qureg2, env)
    end
    destroyQuESTEnv(env)
end

function test_calcHilbertSchmidtDistance()
    env= createQuESTEnv()
    
    for i=1:20
        numQubits = rand(1:5)
        
        targs = Cint.([x for x in 0:numQubits-1])

        M1 = rand(Haar(2), 2^numQubits)
        if Float64 == Float32
            M1 = Matrix{Complex{Float64}}(M1)
        end
        U1 = make_QuEST_matrix(M1)

        M2 = rand(Haar(2), 2^numQubits)
        if Float64 == Float32
            M2 = Matrix{Complex{Float64}}(M2)
        end
        U2 = make_QuEST_matrix(M2)

        v = zeros(Float64, 2^numQubits)
        v[1]=1.0

        state1 = M1*v
        state2 = M2*v

        my_inner_prod = sqrt(sum([abs( state1[i]*conj(state1[j]) - state2[i]*conj(state2[j]))^2 for i=1:2^numQubits for j=1:2^numQubits]))

        qureg1 = createDensityQureg(numQubits, env)
        qureg2 = createDensityQureg(numQubits, env)
        
        multiQubitUnitary(qureg1, targs, U1)
        multiQubitUnitary(qureg2, targs, U2)

        @test real(my_inner_prod) ≈ calcHilbertSchmidtDistance(qureg1, qureg2) atol = tolerance

        destroyQureg(qureg1, env)
        destroyQureg(qureg2, env)

    end
    destroyQuESTEnv(env)
end

function test_calcInnerProduct()
    env= createQuESTEnv()
    
    for i=1:20
        numQubits = rand(1:5)
        
        targs = Cint.([x for x in 0:numQubits-1])

        M1 = rand(Haar(2), 2^numQubits)
        if Float64 == Float32
            M1 = Matrix{Complex{Float64}}(M1)
        end
        U1 = make_QuEST_matrix(M1)

        M2 = rand(Haar(2), 2^numQubits)
        if Float64 == Float32
            M2 = Matrix{Complex{Float64}}(M2)
        end
        U2 = make_QuEST_matrix(M2)

        v = zeros(Float64, 2^numQubits)
        v[1]=1.0

        state1 = M1*v
        state2 = M2*v

        my_inner_prod = transpose(conj(state1))*state2

        qureg1 = createQureg(numQubits, env)
        qureg2 = createQureg(numQubits, env)
        
        multiQubitUnitary(qureg1, targs, U1)
        multiQubitUnitary(qureg2, targs, U2)

        @test my_inner_prod ≈ calcInnerProduct(qureg1, qureg2) atol = tolerance

        destroyQureg(qureg1, env)
        destroyQureg(qureg2, env)

    end
    destroyQuESTEnv(env)
end

function test_calcProbOfOutcome()
    env= createQuESTEnv()
    zero_zero = Matrix([1.0 0;0 0])
    one_one = Matrix([0 0;0 1.0])
    for i=1:20
        numQubits = rand(1:12)
        
        target = rand(0:numQubits-1)

        qureg = createQureg(numQubits, env)
        
        state = rand(Complex{Float64}, 2^numQubits)
        state /= norm(state)

        initStateFromAmps(qureg, real.(state), imag.(state))
        
        Id_right = 1.0*Matrix(I,2^target, 2^target)
        Id_left = 1.0*Matrix(I,2^(numQubits - target -1), 2^(numQubits - target -1))
        
        M_zero = kron(Id_left, zero_zero, Id_right)
        M_one = kron(Id_left, one_one, Id_right)

        prob_zero = transpose(conj(state))*M_zero*state
        prob_one = transpose(conj(state))*M_one*state

        @test real(prob_zero) ≈ calcProbOfOutcome(qureg, target, 0) atol = tolerance
        @test real(prob_one) ≈ calcProbOfOutcome(qureg, target, 1) atol = tolerance

        destroyQureg(qureg, env)

    end
    destroyQuESTEnv(env)
end

function test_calcPurity()
    env= createQuESTEnv()
    
    for i=1:20
        numQubits = rand(1:5)
        
        targs = Cint.([x for x in 0:numQubits-1])

        M1 = rand(Haar(2), 2^numQubits)
        if Float64 == Float32
            M1 = Matrix{Complex{Float64}}(M1)
        end
        U1 = make_QuEST_matrix(M1)

        qureg1 = createDensityQureg(numQubits, env)
        
        multiQubitUnitary(qureg1, targs, U1)

        @test 1 ≈ calcPurity(qureg1) atol = tolerance

        destroyQureg(qureg1, env)

    end
    destroyQuESTEnv(env)
end

function test_calcTotalProb()
    env= createQuESTEnv()
    
    for i=1:20
        numQubits = rand(1:5)
        
        targs = Cint.([x for x in 0:numQubits-1])

        M1 = rand(Haar(2), 2^numQubits)
        if Float64 == Float32
            M1 = Matrix{Complex{Float64}}(M1)
        end
        U1 = make_QuEST_matrix(M1)

        qureg1 = createQureg(numQubits, env)
        
        multiQubitUnitary(qureg1, targs, U1)

        @test 1 ≈ calcTotalProb(qureg1) atol = tolerance

        destroyQureg(qureg1, env)

    end
    destroyQuESTEnv(env)
end

function test_amp_funcs()
    env= createQuESTEnv()
    for i=1:20
        numQubits = rand(1:12)
        
        qureg = createQureg(numQubits, env)
        
        state = rand(Complex{Float64}, 2^numQubits)
        state /= norm(state)

        initStateFromAmps(qureg, real.(state), imag.(state))
        
        for ind =0:2^numQubits-1
            @test getAmp(qureg, ind) ≈ state[ind+1] atol = tolerance
            @test getImagAmp(qureg, ind) ≈ imag(state[ind+1]) atol = tolerance
            @test getRealAmp(qureg, ind) ≈ real(state[ind+1]) atol = tolerance
            @test getProbAmp(qureg, ind) ≈ abs(state[ind+1])^2 atol = tolerance
            @test getNumQubits(qureg) == numQubits
            @test getNumAmps(qureg) == 2^numQubits
        end
        destroyQureg(qureg, env)

    end
    destroyQuESTEnv(env)
end

function test_getDensityAmp()
    env= createQuESTEnv()
    
    for i=1:20
        numQubits = rand(1:5)
        
        targs = Cint.([x for x in 0:numQubits-1])

        M1 = rand(Haar(2), 2^numQubits)
        if Float64 == Float32
            M1 = Matrix{Complex{Float64}}(M1)
        end
        U1 = make_QuEST_matrix(M1)

        v = zeros(Float64, 2^numQubits)
        v[1]=1.0

        state1 = M1*v

        qureg1 = createDensityQureg(numQubits, env)
        
        multiQubitUnitary(qureg1, targs, U1)

        for i=0:2^numQubits-1
            for j=0:2^numQubits-1
                @test getDensityAmp(qureg1, i, j) ≈ state1[i+1]*conj(state1[j+1]) atol = tolerance
            end
        end
        destroyQureg(qureg1, env)
    end
    destroyQuESTEnv(env)
end

function test_cloneQureg()
    env= createQuESTEnv()
    for i=1:20
        numQubits = rand(1:12)
        
        qureg1 = createQureg(numQubits, env)
        qureg2 = createQureg(numQubits, env)
        
        state = rand(Complex{Float64}, 2^numQubits)
        state /= norm(state)

        initStateFromAmps(qureg1, real.(state), imag.(state))

        cloneQureg(qureg2, qureg1)
        
        
        for ind =0:2^numQubits-1
            @test getAmp(qureg1, ind) ≈ getAmp(qureg2, ind) atol = tolerance
        end
        destroyQureg(qureg1, env)
        destroyQureg(qureg2, env)

    end
    destroyQuESTEnv(env)
end

function test_initBlankState()
    env= createQuESTEnv()
    for i=1:20
        numQubits = rand(1:12)
        
        qureg = createQureg(numQubits, env)
        
        initBlankState(qureg)
        
        
        for ind =0:2^numQubits-1
            @test getAmp(qureg, ind) ≈ 0 atol = tolerance
        end
        destroyQureg(qureg, env)

    end
    destroyQuESTEnv(env)
end

function test_initClassicalState()
    env= createQuESTEnv()
    for i=1:20
        numQubits = rand(1:12)
        
        qureg = createQureg(numQubits, env)
        state = rand(0:2^numQubits-1)

        initClassicalState(qureg, state)
        
        
        for ind =0:2^numQubits-1
            if ind == state
                @test getAmp(qureg, ind) ≈ 1 atol = tolerance
            else
                @test getAmp(qureg, ind) ≈ 0 atol = tolerance
            end
        end
        destroyQureg(qureg, env)

    end
    destroyQuESTEnv(env)
end

function test_initPlusState()
    env= createQuESTEnv()
    for i=1:20
        numQubits = rand(1:12)
        
        qureg = createQureg(numQubits, env)

        initPlusState(qureg)
        
        
        for ind =0:2^numQubits-1
            @test getAmp(qureg, ind) ≈ 1/sqrt(2^numQubits) atol = tolerance
        end
        destroyQureg(qureg, env)

    end
    destroyQuESTEnv(env)
end

function test_initPureState()
    env= createQuESTEnv()
    for i=1:20
        numQubits = rand(1:12)
        
        qureg1 = createQureg(numQubits, env)
        qureg2 = createQureg(numQubits, env)
        
        state = rand(Complex{Float64}, 2^numQubits)
        state /= norm(state)

        initStateFromAmps(qureg1, real.(state), imag.(state))

        initPureState(qureg2, qureg1)
        
        
        for ind =0:2^numQubits-1
            @test getAmp(qureg1, ind) ≈ getAmp(qureg2, ind) atol = tolerance
        end
        destroyQureg(qureg1, env)
        destroyQureg(qureg2, env)

    end
    destroyQuESTEnv(env)
end

function test_initStateFromAmps()
    env= createQuESTEnv()
    for i=1:20
        numQubits = rand(1:12)
        
        qureg1 = createQureg(numQubits, env)
        
        state = rand(Complex{Float64}, 2^numQubits)
        state /= norm(state)

        initStateFromAmps(qureg1, real.(state), imag.(state))
        
        for ind =0:2^numQubits-1
            @test getAmp(qureg1, ind) ≈ state[ind+1] atol = tolerance
        end
        destroyQureg(qureg1, env)

    end
    destroyQuESTEnv(env)
end

function test_initZeroState()
    env= createQuESTEnv()
    for i=1:20
        numQubits = rand(1:12)
        
        qureg1 = createQureg(numQubits, env)
        
        

        initZeroState(qureg1)
        
        @test getAmp(qureg1, 0) ≈ 1 atol = tolerance

        for ind =1:2^numQubits-1
            @test getAmp(qureg1, ind) ≈ 0 atol = tolerance
        end
        destroyQureg(qureg1, env)

    end
    destroyQuESTEnv(env)
end

function test_setAmps()
    env= createQuESTEnv()
    for i=1:20
        numQubits = rand(1:12)
        
        qureg1 = createQureg(numQubits, env)
        
        state = rand(Complex{Float64}, 2^numQubits)
        state /= norm(state)

        initStateFromAmps(qureg1, real.(state), imag.(state))

        start_ind = rand(0:2^numQubits-2)
        finish_ind = rand(start_ind:2^numQubits-1)
        numAmps = finish_ind - start_ind +1

        new_amps = rand(Complex{Float64}, numAmps)

        setAmps(qureg1, start_ind, real.(new_amps), imag.(new_amps), numAmps)
        
        for ind =0:2^numQubits-1
            if ind < start_ind
                @test getAmp(qureg1, ind) ≈ state[ind+1] atol = tolerance
            elseif ind <= finish_ind
                @test getAmp(qureg1, ind) ≈ new_amps[ind - start_ind+1] atol = tolerance
            else
                @test getAmp(qureg1, ind) ≈ state[ind+1] atol = tolerance
            end
        end
        destroyQureg(qureg1, env)

    end
    destroyQuESTEnv(env)
end

function test_setWeightedQureg()
    env= createQuESTEnv()
    for i=1:10
        numQubits = rand(1:12)
        
        qureg1 = createQureg(numQubits, env)
        qureg2 = createQureg(numQubits, env)
        qureg3 = createQureg(numQubits, env)
        
        state1 = rand(Complex{Float64}, 2^numQubits)
        state1 /= norm(state1)

        state2 = rand(Complex{Float64}, 2^numQubits)
        state2 /= norm(state2)

        state3 = rand(Complex{Float64}, 2^numQubits)
        state3 /= norm(state3)

        initStateFromAmps(qureg1, real.(state1), imag.(state1))
        initStateFromAmps(qureg2, real.(state2), imag.(state2))
        initStateFromAmps(qureg3, real.(state3), imag.(state3))

        fac1 = rand(Complex{Float64})
        fac2 = rand(Complex{Float64})
        fac3 = rand(Complex{Float64})

        setWeightedQureg(fac1, qureg1, fac2, qureg2, fac3, qureg3)
        
        for ind =0:2^numQubits-1
            rhs  = fac1 * state1[ind+1]
            rhs += fac2 * state2[ind+1]
            rhs += fac3 * state3[ind+1]    
            @test getAmp(qureg3, ind) ≈ rhs atol = tolerance
        end
        
        destroyQureg(qureg1, env)
        destroyQureg(qureg2, env)
        destroyQureg(qureg3, env)

    end
    destroyQuESTEnv(env)
end

function test_QASM()
    env= createQuESTEnv()
    for i=1:10
        numQubits = rand(3:12)
        

        qureg = createQureg(numQubits, env)
        startRecordingQASM(qureg)
        qubits = rand(0:numQubits-1, 10)
        θ = Float64(rand(0.0:0.00001:1))
        α=rand(Complex{Float64})
        β=rand(Complex{Float64})
        norm = sqrt(α*conj(α) + β*conj(β))
        α /= norm
        β /= norm

        hadamard(qureg, qubits[1])
        tGate(qureg, qubits[2])
        sGate(qureg, qubits[3])
        pauliX(qureg, qubits[4])
        pauliY(qureg, qubits[5])
        pauliZ(qureg, qubits[6])
        rotateX(qureg, qubits[7], θ)
        rotateY(qureg, qubits[8], θ)
        rotateZ(qureg, qubits[9], θ)

        original_stdout = stdout
        read_pipe, write_pipe = redirect_stdout()
        lines = Vector{String}()
        printRecordedQASM(qureg)
        Base.Libc.flush_cstdio()
        redirect_stdout(original_stdout)
        close(write_pipe)
        for l =1:13
            push!(lines, readline(read_pipe))
        end
        
        
        @test lines[1] == "OPENQASM 2.0;"
        @test lines[2] == "qreg q[$numQubits];"
        @test lines[3] == "creg c[$numQubits];"
        @test lines[4] == "h q[$(qubits[1])];"
        @test lines[5] == "t q[$(qubits[2])];"
        @test lines[6] == "s q[$(qubits[3])];"
        @test lines[7] == "x q[$(qubits[4])];"
        @test lines[8] == "y q[$(qubits[5])];"
        @test lines[9] == "z q[$(qubits[6])];"
        @test lines[10] == "Rx($θ) q[$(qubits[7])];"
        @test lines[11] == "Ry($θ) q[$(qubits[8])];"
        @test lines[12] == "Rz($θ) q[$(qubits[9])];"
        @test lines[13] == ""
        
        stopRecordingQASM(qureg)
        pauliX(qureg, qubits[1])

        original_stdout = stdout
        read_pipe, write_pipe = redirect_stdout()
        lines = Vector{String}()
        printRecordedQASM(qureg)
        Base.Libc.flush_cstdio()
        redirect_stdout(original_stdout)
        close(write_pipe)
        for l =1:13
            push!(lines, readline(read_pipe))
        end
        

        @test lines[1] == "OPENQASM 2.0;"
        @test lines[2] == "qreg q[$numQubits];"
        @test lines[3] == "creg c[$numQubits];"
        @test lines[4] == "h q[$(qubits[1])];"
        @test lines[5] == "t q[$(qubits[2])];"
        @test lines[6] == "s q[$(qubits[3])];"
        @test lines[7] == "x q[$(qubits[4])];"
        @test lines[8] == "y q[$(qubits[5])];"
        @test lines[9] == "z q[$(qubits[6])];"
        @test lines[10] == "Rx($θ) q[$(qubits[7])];"
        @test lines[11] == "Ry($θ) q[$(qubits[8])];"
        @test lines[12] == "Rz($θ) q[$(qubits[9])];"
        @test lines[13] == ""

        startRecordingQASM(qureg)
        pauliZ(qureg, qubits[1])
        controlledNot(qureg, 0, 1)

        original_stdout = stdout
        read_pipe, write_pipe = redirect_stdout()
        lines = Vector{String}()
        printRecordedQASM(qureg)
        Base.Libc.flush_cstdio()
        redirect_stdout(original_stdout)
        close(write_pipe)
        for l =1:15
            push!(lines, readline(read_pipe))
        end

        @test lines[1] == "OPENQASM 2.0;"
        @test lines[2] == "qreg q[$numQubits];"
        @test lines[3] == "creg c[$numQubits];"
        @test lines[4] == "h q[$(qubits[1])];"
        @test lines[5] == "t q[$(qubits[2])];"
        @test lines[6] == "s q[$(qubits[3])];"
        @test lines[7] == "x q[$(qubits[4])];"
        @test lines[8] == "y q[$(qubits[5])];"
        @test lines[9] == "z q[$(qubits[6])];"
        @test lines[10] == "Rx($θ) q[$(qubits[7])];"
        @test lines[11] == "Ry($θ) q[$(qubits[8])];"
        @test lines[12] == "Rz($θ) q[$(qubits[9])];"
        @test lines[13] == "z q[$(qubits[1])];"
        @test lines[14] == "cx q[0],q[1];"
        @test lines[15] == ""
        
        stopRecordingQASM(qureg)
        clearRecordedQASM(qureg)
        startRecordingQASM(qureg)

        hadamard(qureg, qubits[1])
        tGate(qureg, qubits[2])
        sGate(qureg, qubits[3])

        original_stdout = stdout
        read_pipe, write_pipe = redirect_stdout()
        lines = Vector{String}()
        printRecordedQASM(qureg)
        Base.Libc.flush_cstdio()
        redirect_stdout(original_stdout)
        close(write_pipe)
        for l =1:7
            push!(lines, readline(read_pipe))
        end

        @test lines[1] == "h q[$(qubits[1])];"
        @test lines[2] == "t q[$(qubits[2])];"
        @test lines[3] == "s q[$(qubits[3])];"
        @test lines[4] == ""


        writeRecordedQASMToFile(qureg, "test_QASM_file.txt")

        io = open("test_QASM_file.txt", "r")

        @test readline(io) == "h q[$(qubits[1])];"
        @test readline(io) == "t q[$(qubits[2])];"
        @test readline(io) == "s q[$(qubits[3])];"
        @test readline(io) == ""

        close(io)

        rm("test_QASM_file.txt")

        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end
#

function test_GPU()
    env= createQuESTEnv()
    for i=1:1
        numQubits = rand(3:12)
        qureg = createQureg(numQubits, env)
        state = rand(Complex{Float64}, 2^numQubits)
        state /= norm(state)
        initStateFromAmps(qureg, real.(state), imag.(state))
        copyStateToGPU(qureg)
        copyStateFromGPU(qureg)
        #env_str = getEnvironmentString(env, qureg)
        #print(env_str)
        print("Add more stuff here to test GPU")

    end

    destroyQuESTEnv(env)
end

function test_initDebugState()
    env= createQuESTEnv()
    for i=1:10
        numQubits = rand(3:12)
        qureg = createQureg(numQubits, env)
        initDebugState(qureg)
        state_vec = unsafe_load_state_vec(qureg.stateVec, 2^numQubits)
        for i=0:2^numQubits-1
            @test reals[i+1] ≈ 2i/10 atol = tolerance
            @test imags[i+1] ≈ (2i+1)/10 atol = tolerance
        end

    end

    destroyQuESTEnv(env)
end

function test_reportPauliHamil()
    paulies = [QuEST_Types.PAULI_I, 
               QuEST_Types.PAULI_X, 
               QuEST_Types.PAULI_Y, 
               QuEST_Types.PAULI_Z]
    env= createQuESTEnv()
    for i=1:10

        numQubits = rand(3:12)
        numSumTerms = rand(1:20)
        
        
        hamil = createPauliHamil(numQubits, numSumTerms)
        codes = rand(paulies, numQubits*numSumTerms)
        coeffs = rand(Float64, numSumTerms)
        initPauliHamil(hamil, coeffs, codes)


        original_stdout = stdout
        read_pipe, write_pipe = redirect_stdout()
        lines = Vector{String}()
        reportPauliHamil(hamil)
        Base.Libc.flush_cstdio()
        redirect_stdout(original_stdout)
        close(write_pipe)
        for l =1:numSumTerms+1
            push!(lines, readline(read_pipe))
        end

        for row = 0:numSumTerms-1
            if occursin("\t", lines[row+1])
                coeff_str = split(lines[row+1], "\t")[1]
                codes_str = split(lines[row+1], "\t")[2]
                codes_list = split(codes_str, " ")
            else
                coeff_str = split(lines[row+1], " ")[1]
                codes_list = split(lines[row+1], " ")[2:end]
            end
            
            #print(codes_list)
            @test length(codes_list) == (numQubits+1)
            @test isapprox(parse(Float64, coeff_str), coeffs[row+1], atol=1e-5)
            for q = 1:numQubits
                @test paulies[parse(Cint, codes_list[q])+1] == codes[row*numQubits + q]
            end
        end

        @test lines[end] == ""

        

        destroyPauliHamil(hamil)
    end

    destroyQuESTEnv(env)
end

function test_reportQuregParams()

    env= createQuESTEnv()
    numQubits = rand(3:12)
    qureg = createQureg(numQubits, env)


    original_stdout = stdout
    read_pipe, write_pipe = redirect_stdout()
    lines = Vector{String}()

    reportQuregParams(qureg)
    Base.Libc.flush_cstdio()
    redirect_stdout(original_stdout)
    close(write_pipe)

    for l =1:5
        push!(lines, readline(read_pipe))
    end

    @test lines[1] == "QUBITS:"
    @test lines[2] == "Number of qubits is $numQubits."
    @test lines[3] == "Number of amps is $(2^numQubits)."
    @test lines[5] == ""

    destroyQureg(qureg, env)
    destroyQuESTEnv(env)
end

function test_reporQuESTEnv()
    
    env= createQuESTEnv()
    numQubits = rand(3:12)
    qureg = createQureg(numQubits, env)

    original_stdout = stdout
    read_pipe, write_pipe = redirect_stdout()
    lines = Vector{String}()

    reportQuESTEnv(env)
    Base.Libc.flush_cstdio()
    redirect_stdout(original_stdout)

    close(write_pipe)
    
    for l =1:7
        push!(lines, readline(read_pipe))
    end

    @test lines[1] == "EXECUTION ENVIRONMENT:"
    @test lines[end] == ""

    destroyQureg(qureg, env)
    destroyQuESTEnv(env)
end

function test_reportState()
    env= createQuESTEnv()
    for i=1:1
        numQubits = rand(3:5)
        qureg = createQureg(numQubits, env)
        state = rand(Complex{Float64}, 2^numQubits)
        state /= norm(state)
        initStateFromAmps(qureg, real.(state), imag.(state))

        original_stdout = stdout
        read_pipe, write_pipe = redirect_stdout()
        lines = Vector{String}()
        reportStateToScreen(qureg, env, 1)
        Base.Libc.flush_cstdio()
        redirect_stdout(original_stdout)
        close(write_pipe)
        for l =1:2^numQubits+4
            push!(lines, readline(read_pipe))
        end
        @test lines[2] == "real, imag"
        @test lines[end-1] == "]"
        @test lines[end] == ""
        for l =1:2^numQubits
            this_re = parse(Float64, split(lines[l+2], " ")[1][1:end-1])
            this_im = parse(Float64, split(lines[l+2], " ")[2])
            @test this_im ≈ imag(state[l]) atol = tolerance
            @test this_re ≈ real(state[l]) atol = tolerance
        end

        reportState(qureg)

        io = open("state_rank_0.csv", "r")

        @test readline(io) == "real, imag"

        for l=1:2^numQubits
            this_line = readline(io)
            this_re = parse(Float64, split(this_line, " ")[1][1:end-1])
            this_im = parse(Float64, split(this_line, " ")[2])
            @test this_im ≈ imag(state[l]) atol = tolerance
            @test this_re ≈ real(state[l]) atol = tolerance
        end

        @test readline(io) == ""
        
        close(io)

        rm("state_rank_0.csv")

    end

    destroyQuESTEnv(env)
    
end

function test_seed()
    env= createQuESTEnv()

    num_seeds = rand(1:30)
    seeds = rand(Culong, num_seeds)
    seedQuESTDefault()
    seedQuEST(seeds)

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
