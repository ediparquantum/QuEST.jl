using Pkg
Pkg.activate(".")
include("src/QuEST.jl")
using .QuEST
using Test
tolerance = 1e-10

complex
qComplex

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


Complex(α)

function test_controlledCompactUnitary()
    env= createQuESTEnv()
    for t=1:10
        α=rand(Complex{Float64})
        β=rand(Complex{Float64})
        norm = sqrt(α*conj(α) + β*conj(β))
        α /= norm
        β /= norm
        numQubits = rand(2:12)
        target = rand(1:numQubits)
        control = rand(1:numQubits)
        while target == control
            target = rand(1:numQubits)
        end
        qureg = createQureg(numQubits, env)
        pauliX(qureg, control)
        α = QComplex(α.re, α.im)
        β = QComplex(β.re, β.im)
        controlledCompactUnitary(qureg, control, target, α, β)
        reals = unsafe_wrap(Vector{Float64}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{Float64}, qureg.stateVec.imag, 2^numQubits)
        #print(reals)
        #print(imags)
        @test reals[2^(control)+1] ≈ real(α) atol = tolerance
        @test imags[2^(control-1)+1] ≈ imag(α) atol = tolerance
        @test reals[2^(control-1)+2^(target-1)+1] ≈ real(β) atol = tolerance
        @test imags[2^(control-1)+2^(target-1)+1] ≈ imag(β) atol = tolerance
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end
