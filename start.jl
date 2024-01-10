using Pkg
Pkg.activate(".")

include("src/QuEST.jl")
using .QuEST

using Test


tolerance = 1e-10



function test_controlledCompactUnitary()
    env = createQuESTEnv()
    for t=1:10
        α=Base.Complex(rand(Cdouble),rand(Cdouble))
        β=Base.Complex(rand(Cdouble),rand(Cdouble))
        norm = sqrt(α*conj(α) + β*conj(β))
        α /= norm
        β /= norm
        @test abs(α)^2 + abs(β)^2 ≈ 1.0 atol = tolerance
        numQubits = 4
        target = 4
        control = 3
        qureg = createQureg(numQubits, env)
        α = QuEST.Complex(α.re, α.im)
        β = QuEST.Complex(β.re, β.im)
        controlledCompactUnitary(qureg, control, target, α, β)
        unsafe_load_state_vec(qureg.stateVec,2^numQubits)
        reals = unsafe_wrap(Base.Vector{Float64}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Base.Vector{Float64}, qureg.stateVec.imag, 2^numQubits)
    
        @test reals[2^(control-1)+1] ≈ α.real atol = tolerance
        @test imags[2^(control-1)+1] ≈ α.imag atol = tolerance
        @test reals[2^(control-1)+2^(target-1)+1] ≈ β.real atol = tolerance
        @test imags[2^(control-1)+2^(target-1)+1] ≈ β.imag atol = tolerance
        destroyQureg(qureg, env)
    end
    destroyQuESTEnv(env)
end


function test_controlledCompactUnitary()
    env= QuEST.createQuESTEnv()
    for t=1:10
        α=rand(Complex{qreal})
        β=rand(Complex{qreal})
        norm = sqrt(α*conj(α) + β*conj(β))
        α /= norm
        β /= norm
        numQubits = rand(2:12)
        target = rand(1:numQubits)
        control = rand(1:numQubits)
        while target == control
            target = rand(1:numQubits)
        end
        qureg = QuEST.createQureg(numQubits, env)
        QuEST.pauliX(qureg, control-1)
        QuEST.controlledCompactUnitary(qureg, control-1, target-1, α, β)
        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
        #print(reals)
        #print(imags)
        @test reals[2^(control-1)+1] ≈ real(α) atol = tolerance
        @test imags[2^(control-1)+1] ≈ imag(α) atol = tolerance
        @test reals[2^(control-1)+2^(target-1)+1] ≈ real(β) atol = tolerance
        @test imags[2^(control-1)+2^(target-1)+1] ≈ imag(β) atol = tolerance
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end
