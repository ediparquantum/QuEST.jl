# QuEST.jl/test/runtests.jl
#
# Authors:
#  - Dirk Oliver Theis, Ketita Labs
#
# Copyright (c) 2020 Ketita Labs oü, Tartu, Estonia
#
# MIT License
#
# Copyright (c) 2020 Ketita Labs oü, Tartu, Estonia
#
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

using QuEST_jl
import QuEST_jl.QuEST64
QuEST = QuEST_jl.QuEST64
Float64 = QuEST.QuEST_Types.Float64
using Test
using RandomMatrices
using LinearAlgebra
using Printf
tolerance = 1e-10



function test_env()
    env1 = QuEST.createQuESTEnv()
    env2 = QuEST.createQuESTEnv()
    QuEST.destroyQuESTEnv(env1)
    QuEST.destroyQuESTEnv(env2)
end

function test_qureg()
    env= QuEST.createQuESTEnv()
    for i=1:10
    qureg1 = QuEST.createQureg(i, env)
    qureg2 = QuEST.createQureg(i, env)    
    @test qureg1.numQubitsRepresented == i
    @test qureg2.numQubitsRepresented == i
    QuEST.destroyQureg(qureg1, env) 
    QuEST.destroyQureg(qureg2, env)
    end
    QuEST.destroyQuESTEnv(env)

end

function test_qureg_density()
    env= QuEST.createQuESTEnv()
    for i=1:10
    qureg1 = QuEST.createDensityQureg(i, env)
    qureg2 = QuEST.createDensityQureg(i, env)    
    @test qureg1.numQubitsRepresented == i
    @test qureg2.numQubitsRepresented == i
    QuEST.destroyQureg(qureg1, env) 
    QuEST.destroyQureg(qureg2, env)
    end
    QuEST.destroyQuESTEnv(env)

end

function check_createCloneQureg()
    num_tests = 10
    env= QuEST.createQuESTEnv()
    for i =1:num_tests
        num_qubits = rand(2:12)
        qureg1 = QuEST.createQureg(num_qubits, env)
        for qubit =0:num_qubits-1
            QuEST.rotateX(qureg1, qubit, rand(Float64))
        end
        for qubit = 0:num_qubits-2
            QuEST.controlledNot(qureg1, qubit, qubit+1)
        end
        qureg2 = QuEST.createCloneQureg(qureg1, env)
        @test qureg1.numQubitsRepresented == qureg2.numQubitsRepresented
        @test qureg1.numAmpsTotal == qureg2.numAmpsTotal
        @test QuEST.getNumAmps(qureg1) == QuEST.getNumAmps(qureg2)
        @test QuEST.getNumQubits(qureg1) == QuEST.getNumQubits(qureg2)
        for ind=1:2^num_qubits-1
            @test QuEST.getProbAmp(qureg1, ind) == QuEST.getProbAmp(qureg2, ind)
            @test QuEST.getRealAmp(qureg1, ind) == QuEST.getRealAmp(qureg2, ind)
            @test QuEST.getImagAmp(qureg1, ind) == QuEST.getImagAmp(qureg2, ind)
        end 
    end
    QuEST.destroyQuESTEnv(env)
end

function check_ComplexMatrixN()
    for i=1:10
        M = QuEST.createComplexMatrixN(i)
        @test M.numQubits == i
        M_j_real = rand(qreal, 2^i * 2^i)
        M_j_imag = rand(qreal, 2^i * 2^i)
        QuEST.initComplexMatrixN(M, M_j_real, M_j_imag)
        real_ptrs = unsafe_wrap(Vector{Ptr{qreal}}, M.real, 2^i)
        imag_ptrs = unsafe_wrap(Vector{Ptr{qreal}}, M.imag, 2^i)
        for r=1:2^i
            real_row = unsafe_wrap(Vector{qreal}, real_ptrs[r], 2^i)
            imag_row = unsafe_wrap(Vector{qreal}, imag_ptrs[r], 2^i)
            for c=1:2^i
                @test M_j_real[(r-1)*2^i + c] ≈ real_row[c] atol = tolerance
                @test M_j_imag[(r-1)*2^i + c] ≈ imag_row[c] atol = tolerance
            end
        end
        QuEST.destroyComplexMatrixN(M)
    end

    for i=1:10
        M_j = rand(Base.Complex{qreal}, 2^i, 2^i)
        M = QuEST.make_QuEST_matrix(M_j)
        real_ptrs = unsafe_wrap(Vector{Ptr{qreal}}, M.real, 2^i)
        imag_ptrs = unsafe_wrap(Vector{Ptr{qreal}}, M.imag, 2^i)
        for r=1:2^i
            real_row = unsafe_wrap(Vector{qreal}, real_ptrs[r], 2^i)
            imag_row = unsafe_wrap(Vector{qreal}, imag_ptrs[r], 2^i)
            for c=1:2^i
                @test real(M_j[r, c]) ≈ real_row[c] atol = tolerance
                @test imag(M_j[r, c]) ≈ imag_row[c] atol = tolerance
            end
        end
        QuEST.destroyComplexMatrixN(M)
    end

    for i=1:10
        m_func(r, c) = 0.12121r +0.343434c*im
        M = QuEST.createComplexMatrixN(i)
        QuEST.fill_ComplexMatrix!(M, m_func)
        real_ptrs = unsafe_wrap(Vector{Ptr{qreal}}, M.real, 2^i)
        imag_ptrs = unsafe_wrap(Vector{Ptr{qreal}}, M.imag, 2^i)
        for r=1:2^i
            real_row = unsafe_wrap(Vector{qreal}, real_ptrs[r], 2^i)
            imag_row = unsafe_wrap(Vector{qreal}, imag_ptrs[r], 2^i)
            for c=1:2^i
                @test real(m_func(r, c)) ≈ real_row[c] atol = tolerance
                @test imag(m_func(r, c)) ≈ imag_row[c] atol = tolerance
            end
        end
        QuEST.destroyComplexMatrixN(M)
    end
end

function test_DiagonalOp()
    env= QuEST.createQuESTEnv()
    for i=1:10
        num_qubits = rand(1:12)
        op = QuEST.createDiagonalOp(num_qubits, env)
        @test op.numQubits == num_qubits
        op_j = rand(Complex{qreal}, 2^num_qubits)
        QuEST.initDiagonalOp(op, real(op_j), imag(op_j))
        reals = unsafe_wrap(Vector{qreal}, op.real, 2^num_qubits)
        imags = unsafe_wrap(Vector{qreal}, op.imag, 2^num_qubits)
        for e = 1:2^num_qubits
            @test reals[e] ≈ real(op_j[e]) atol = tolerance
            @test imags[e] ≈ imag(op_j[e]) atol = tolerance
        end

        
        ind1 = rand(1:2^num_qubits)
        ind2 = rand(1:2^num_qubits)
        start_ind = min(ind1, ind2)
        len_arr = abs(ind1-ind2)+1
        op_j = rand(Complex{qreal}, len_arr)
        test_re = real(op_j)
        test_im = imag(op_j)
        QuEST.setDiagonalOpElems(op, start_ind-1, test_re, test_im, len_arr)
        for ind = start_ind:start_ind+len_arr-1
            @test reals[ind] ≈ test_re[ind-start_ind+1] atol = tolerance
            @test imags[ind] ≈ test_im[ind-start_ind+1] atol = tolerance
        end
        QuEST.syncDiagonalOp(op)
        QuEST.destroyDiagonalOp(op, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_PauliHaiml()
    env= QuEST.createQuESTEnv()
    paulies = [QuEST.QuEST_Types.PAULI_I, 
               QuEST.QuEST_Types.PAULI_X, 
               QuEST.QuEST_Types.PAULI_Y, 
               QuEST.QuEST_Types.PAULI_Z]
    for i=1:10
        numQubits = rand(1:12)
        numSumTerms = rand(1:20)
        hamil = QuEST.createPauliHamil(numQubits, numSumTerms)
        @test hamil.numSumTerms == numSumTerms
        @test hamil.numQubits   == numQubits
        QuEST.destroyPauliHamil(hamil)
        
        hamil = QuEST.createPauliHamil(numQubits, numSumTerms)
        codes = rand(paulies, numQubits*numSumTerms)
        coeffs = rand(qreal, numSumTerms)
        QuEST.initPauliHamil(hamil, coeffs, codes)
        hamil_codes = unsafe_wrap(Vector{QuEST.QuEST_Types.pauliOpType}, hamil.pauliCodes, numQubits*numSumTerms)
        hamil_coeffs = unsafe_wrap(Vector{qreal}, hamil.termCoeffs, numSumTerms)
        for ind =1:numSumTerms
            @test coeffs[ind] ≈ hamil_coeffs[ind]  atol = tolerance
            for qubit = 1:numQubits
                @test codes[(ind-1)*numQubits+qubit] == hamil_codes[(ind-1)*numQubits+qubit]
            end
        end
        QuEST.destroyPauliHamil(hamil)

        
        
        numQubits = rand(1:12)
        numSumTerms = rand(1:20)
        codes = rand(paulies, numQubits*numSumTerms)
        coeffs = rand(qreal, numSumTerms)

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

        hamil = QuEST.createPauliHamilFromFile("hamil")
        hamil_codes = unsafe_wrap(Vector{QuEST.QuEST_Types.pauliOpType}, hamil.pauliCodes, numQubits*numSumTerms)
        hamil_coeffs = unsafe_wrap(Vector{qreal}, hamil.termCoeffs, numSumTerms)
        for ind=1:numSumTerms
            @test coeffs[ind] ≈ hamil_coeffs[ind] atol = tolerance
            for qubit = 1:numQubits
                @test codes[(ind-1)*numQubits+qubit] == hamil_codes[(ind-1)*numQubits+qubit]
            end
        end

        rm("hamil")



    end
    QuEST.destroyQuESTEnv(env)
end


function test_compactUnitary()
    env= QuEST.createQuESTEnv()
    for t=1:10
        α=rand(Complex{qreal})
        β=rand(Complex{qreal})
        norm = sqrt(α*conj(α) + β*conj(β))
        α /= norm
        β /= norm
        numQubits = rand(1:12)
        target = rand(1:numQubits)
        qureg = QuEST.createQureg(numQubits, env)
        QuEST.compactUnitary(qureg, target-1, α, β)
        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
        @test reals[1] ≈ real(α) atol = tolerance
        @test imags[1] ≈ imag(α) atol = tolerance
        @test reals[2^(target-1)+1] ≈ real(β) atol = tolerance
        @test imags[2^(target-1)+1] ≈ imag(β) atol = tolerance
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
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

function test_controlledMultiQubitUnitary()
    env= QuEST.createQuESTEnv()
    for t=1:10
        numQubits = rand(2:12)
        num_targs = rand(1:numQubits-1)

        targs = Cint.([x for x in 1:num_targs])
        targs_rev = Cint.([x for x in num_targs:-1:1])
        control = 0

        
        M_j = rand(Haar(2), 2^num_targs)
        if qreal == Float32
            M_j = Matrix{Complex{qreal}}(M_j)
        end
        M = QuEST.make_QuEST_matrix(M_j)
    
        qureg = QuEST.createQureg(numQubits, env)
        
        QuEST.pauliX(qureg, 0)
        
        QuEST.controlledMultiQubitUnitary(qureg, control, targs, M)
        
        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)

        X=Matrix([0 1.0;1 0])
        res = numQubits-num_targs-1
        Idm = 1.0*Matrix(I,2^res, 2^res)
        U=kron(Idm, M_j, X)
        v = zeros(2^numQubits)
        v[1]=1.0
        v=U*v
        #print(reals)
        #print(imags)
        for ind = 1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end
        QuEST.destroyComplexMatrixN(M)
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
    
end

function test_controlledNot()
    env= QuEST.createQuESTEnv()
    for t=1:10
        
        numQubits = rand(2:12)
        target = rand(1:numQubits)
        control = rand(1:numQubits)
        while target == control
            target = rand(1:numQubits)
        end

        qureg = QuEST.createQureg(numQubits, env)

        QuEST.pauliX(qureg, control-1)
        QuEST.controlledNot(qureg, control-1, target-1)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
        one_ind=2^(control-1)+2^(target-1)+1
        for ind =1:2^numQubits
            if ind == one_ind
                @test reals[ind] ≈ 1 atol = tolerance
                @test imags[ind] ≈ 0 atol = tolerance
            else
                @test reals[ind] ≈ 0 atol = tolerance
                @test imags[ind] ≈ 0 atol = tolerance
            end
        end
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_controlledPauliY()
    env= QuEST.createQuESTEnv()
    for t=1:10
        
        numQubits = rand(2:12)
        target = rand(1:numQubits)
        control = rand(1:numQubits)
        while target == control
            target = rand(1:numQubits)
        end

        qureg = QuEST.createQureg(numQubits, env)

        QuEST.pauliX(qureg, control-1)
        QuEST.controlledPauliY(qureg, control-1, target-1)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
        one_ind=2^(control-1)+2^(target-1)+1
        for ind =1:2^numQubits
            if ind == one_ind
                @test reals[ind] ≈ 0 atol = tolerance
                @test imags[ind] ≈ 1 atol = tolerance
            else
                @test reals[ind] ≈ 0 atol = tolerance
                @test imags[ind] ≈ 0 atol = tolerance
            end
        end
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_controlledPhaseFlip()
    env= QuEST.createQuESTEnv()
    for t=1:10
        
        numQubits = rand(2:12)
        target = rand(1:numQubits)
        control = rand(1:numQubits)
        while target == control
            target = rand(1:numQubits)
        end

        qureg = QuEST.createQureg(numQubits, env)

        QuEST.pauliX(qureg, control-1)
        QuEST.pauliX(qureg, target-1)
        QuEST.controlledPhaseFlip(qureg, control-1, target-1)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
        one_ind=2^(control-1)+2^(target-1)+1
        for ind =1:2^numQubits
            if ind == one_ind
                @test reals[ind] ≈ -1 atol = tolerance
                @test imags[ind] ≈ 0 atol = tolerance
            else
                @test reals[ind] ≈ 0 atol = tolerance
                @test imags[ind] ≈ 0 atol = tolerance
            end
        end
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_controlledPhaseShift()
    env= QuEST.createQuESTEnv()
    for t=1:10
        
        numQubits = rand(2:12)
        target = rand(1:numQubits)
        control = rand(1:numQubits)
        while target == control
            target = rand(1:numQubits)
        end
        θ = rand(qreal)
        qureg = QuEST.createQureg(numQubits, env)

        QuEST.pauliX(qureg, control-1)
        QuEST.pauliX(qureg, target-1)
        QuEST.controlledPhaseShift(qureg, control-1, target-1, θ)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
        one_ind=2^(control-1)+2^(target-1)+1
        for ind =1:2^numQubits
            if ind == one_ind
                @test reals[ind] ≈ real(ℯ^(im*θ)) atol = tolerance
                @test imags[ind] ≈ imag(ℯ^(im*θ)) atol = tolerance
            else
                @test reals[ind] ≈ 0 atol = tolerance
                @test imags[ind] ≈ 0 atol = tolerance
            end
        end
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end



function test_controlledRotateAroundAxis()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(2:12)
        target = rand(1:numQubits)
        control = rand(1:numQubits)
        while target == control
            target = rand(1:numQubits)
        end
        θ = rand(qreal)
        n̂ = rand(qreal, 3)
        n̂ /= norm(n̂)
        U = ℯ^(-im*θ*(n̂[1]*X+n̂[2]*Y+n̂[3]*Z)/2)
        qureg = QuEST.createQureg(numQubits, env)
        v = [1.0, 0.0]
        v = U*v
        QuEST.pauliX(qureg, control-1)
        
        QuEST.controlledRotateAroundAxis(qureg, control-1, target-1, θ, n̂)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
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
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_controlledRotateX()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(2:12)
        target = rand(1:numQubits)
        control = rand(1:numQubits)
        while target == control
            target = rand(1:numQubits)
        end
        θ = rand(qreal)
        
        U = ℯ^(-im*θ*X/2)
        qureg = QuEST.createQureg(numQubits, env)
        v = [1.0, 0.0]
        v = U*v
        QuEST.pauliX(qureg, control-1)
        
        QuEST.controlledRotateX(qureg, control-1, target-1, θ)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
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
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_controlledRotateY()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(2:12)
        target = rand(1:numQubits)
        control = rand(1:numQubits)
        while target == control
            target = rand(1:numQubits)
        end
        θ = rand(qreal)
        
        U = ℯ^(-im*θ*Y/2)
        qureg = QuEST.createQureg(numQubits, env)
        v = [1.0, 0.0]
        v = U*v
        QuEST.pauliX(qureg, control-1)
        
        QuEST.controlledRotateY(qureg, control-1, target-1, θ)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
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
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_controlledRotateZ()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(2:12)
        target = rand(1:numQubits)
        control = rand(1:numQubits)
        while target == control
            target = rand(1:numQubits)
        end
        θ = rand(qreal)
        
        U = ℯ^(-im*θ*Z/2)
        qureg = QuEST.createQureg(numQubits, env)
        v = [1.0, 0.0]
        v = U*v
        QuEST.pauliX(qureg, control-1)
        
        QuEST.controlledRotateZ(qureg, control-1, target-1, θ)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
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
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_controlledTwoQubitUnitary()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(3:12)
    
        qureg = QuEST.createQureg(numQubits, env)
        v = zeros(qreal, 2^numQubits)
        v[1]=1.0
        M_j = rand(Haar(2), 4)
        if qreal == Float32
            M_j = Matrix{Complex{qreal}}(M_j)
        end
        res = numQubits-3
        Idm = 1.0*Matrix(I,2^res, 2^res)
        U=kron(Idm, M_j, X)
        v=U*v
        
        QuEST.pauliX(qureg, 0)
        QuEST.controlledTwoQubitUnitary(qureg, 0, 1, 2, M_j)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
        
        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_controlledUnitary()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(2:12)
        qureg = QuEST.createQureg(numQubits, env)
        target = rand(1:numQubits)
        control = rand(1:numQubits)
        while target == control
            target = rand(1:numQubits)
        end
        M_j = rand(Haar(2), 2)
        if qreal == Float32
            M_j = Matrix{Complex{qreal}}(M_j)
        end
        v = [1.0, 0.0]
        v = M_j*v
        QuEST.pauliX(qureg, control-1)
        
        QuEST.controlledUnitary(qureg, control-1, target-1, M_j)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
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
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_hadamard()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    H=Matrix([1.0 1.0;1.0 -1.0])/sqrt(2)
    for t=1:10
        
        numQubits = rand(1:12)
        qureg = QuEST.createQureg(numQubits, env)
        target = rand(0:numQubits-1)

        v = zeros(qreal, 2^numQubits)
        v[1]=1
        i_right=1.0*Matrix(I, 2^target, 2^target)
        res = numQubits-target-1
        i_left =1.0*Matrix(I, 2^res, 2^res)
        U=kron(i_left, H, i_right)
        v=U*v

        QuEST.hadamard(qureg, target)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)

        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_multiControlledMultiQubitUnitary()
    print("warning check test_multiControlledMultiQubitUnitary\n")
    env= QuEST.createQuESTEnv()
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
        if qreal == Float32
            M_j = Matrix{Complex{qreal}}(M_j)
        end
        M = QuEST.make_QuEST_matrix(M_j)
    
        qureg = QuEST.createQureg(numQubits, env)
        
        for ind in ctrls
            QuEST.pauliX(qureg, ind)
        end
        #QuEST.pauliX(qureg, 0)
        QuEST.multiControlledMultiQubitUnitary(qureg, ctrls, targs, M)
        
        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)

        X=Matrix([0 1.0;1 0])
        big_X=Matrix([0 1.0;1 0])
        for iter = 1:num_ctrls-1
            big_X=kron(big_X, X)
        end
        res = numQubits-num_targs-num_ctrls
        Idm = 1.0*Matrix(I,2^res, 2^res)
        U=kron(Idm, M_j, big_X)
        v = zeros(qreal, 2^numQubits)
        v[1]=1.0
        v=U*v
        #print(reals)
        #print(imags)
        #print(v)
        for ind = 1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end
        QuEST.destroyComplexMatrixN(M)
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
    
end


function test_multiControlledPhaseFlip()
    env= QuEST.createQuESTEnv()
    for t=1:10
        
        numQubits = rand(2:12)
        num_ctrls = rand(1:numQubits)

        ctrls = Cint.([x for x in 0:num_ctrls-1])

        qureg = QuEST.createQureg(numQubits, env)

        for ind = 0:num_ctrls-1
            QuEST.pauliX(qureg, ind)
        end

        QuEST.multiControlledPhaseFlip(qureg, ctrls)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
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
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_multiControlledPhaseShift()
    env= QuEST.createQuESTEnv()
    for t=1:10
        
        numQubits = rand(2:12)
        num_ctrls = rand(1:numQubits)

        ctrls = Cint.([x for x in 0:num_ctrls-1])
        θ = rand(qreal)
        qureg = QuEST.createQureg(numQubits, env)

        for ind in ctrls
            QuEST.pauliX(qureg, ind)
        end
        
        QuEST.multiControlledPhaseShift(qureg, ctrls, θ)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
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
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_multiControlledTwoQubitUnitary()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
    
        
        numQubits = rand(3:12)
        qureg = QuEST.createQureg(numQubits, env)
        num_ctrls = rand(1:numQubits-2)

        ctrls = Cint.([x for x in 0:num_ctrls-1])
        v = zeros(qreal, 2^numQubits)
        v[1]=1.0
        M_j = rand(Haar(2), 4)
        if qreal == Float32
            M_j = Matrix{Complex{qreal}}(M_j)
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
            QuEST.pauliX(qureg, ind)
        end
        QuEST.multiControlledTwoQubitUnitary(qureg, ctrls, num_ctrls, num_ctrls+1, M_j)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
        
        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_multiControlledUnitary()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(3:12)
        qureg = QuEST.createQureg(numQubits, env)
        num_ctrls = rand(1:numQubits-1)

        ctrls = Cint.([x for x in 0:num_ctrls-1])

        M_j = rand(Haar(2), 2)
        if qreal == Float32
            M_j = Matrix{Complex{qreal}}(M_j)
        end
        
        v = zeros(qreal, 2^numQubits)
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
            QuEST.pauliX(qureg, ind)
        end
        
        QuEST.multiControlledUnitary(qureg, ctrls, num_ctrls, M_j)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
        

        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end


        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_multiQubitUnitary()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(3:12)
        qureg = QuEST.createQureg(numQubits, env)
        num_targs = rand(1:numQubits-1)

        targs = Cint.([x for x in 0:num_targs-1])

        M_j = rand(Haar(2), 2^num_targs)
        if qreal == Float32
            M_j = Matrix{Complex{qreal}}(M_j)
        end
        M = QuEST.make_QuEST_matrix(M_j)

        v = zeros(qreal, 2^numQubits)
        v[1]=1.0

        res = numQubits-num_targs
        Idm = 1.0*Matrix(I,2^res, 2^res)
        U=kron(Idm, M_j)
        v=U*v
        
        
        QuEST.multiQubitUnitary(qureg, targs, M)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
        

        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end


        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_multiRotatePauli()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    Id=Matrix([1.0 0;0 1])
    paulies = [QuEST.QuEST_Types.PAULI_I, 
               QuEST.QuEST_Types.PAULI_X, 
               QuEST.QuEST_Types.PAULI_Y, 
               QuEST.QuEST_Types.PAULI_Z]

    pauli_dict=Dict()
    pauli_dict[QuEST.QuEST_Types.PAULI_I]=Id
    pauli_dict[QuEST.QuEST_Types.PAULI_X]=X
    pauli_dict[QuEST.QuEST_Types.PAULI_Y]=Y
    pauli_dict[QuEST.QuEST_Types.PAULI_Z]=Z
    for t=1:10
        
        numQubits = rand(3:12)
        qureg = QuEST.createQureg(numQubits, env)
        
        codes = rand(paulies, numQubits)

        c_codes = Array{QuEST.QuEST_Types.pauliOpType, 1}()
        c_targs = Array{Cint, 1}()
        for i =1:numQubits
            if codes[i] != QuEST.QuEST_Types.PAULI_I
                push!(c_codes, codes[i])
                push!(c_targs, i-1)
            end
        end

        
        M = ones(1,1)
        for ind =1:numQubits
            M=kron(pauli_dict[codes[ind]], M)
        end

        v = zeros(qreal, 2^numQubits)
        v[1]=1.0

        θ = rand(qreal)
        U=ℯ^(-im*θ*M/2)
        v=U*v
        
        
        QuEST.multiRotatePauli(qureg, c_targs, c_codes, θ)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
        

        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end


        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_multiRotateZ()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    Id=Matrix([1.0 0;0 1])
    
    for t=1:10
        
        numQubits = rand(3:12)
        qureg = QuEST.createQureg(numQubits, env)
        
        codes = rand(0:1, numQubits)
        θ = rand(qreal)
        
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

        v = zeros(qreal, 2^numQubits)
        v[1]=1.0

        
        U=ℯ^(-im*θ*M/2)
        v=U*v
        
        
        QuEST.multiRotateZ(qureg, c_targs, θ)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
        

        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end


        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end


function test_multiStateControlledUnitary()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    Id=Matrix([1.0 0;0 1])
    for t=1:10
        
        numQubits = rand(3:12)
        qureg = QuEST.createQureg(numQubits, env)
        
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
        if qreal == Float32
            M_j = Matrix{Complex{qreal}}(M_j)
        end
        
        v = zeros(qreal, 2^numQubits)
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
                QuEST.pauliX(qureg, ind-1)
            end
        end
        QuEST.multiStateControlledUnitary(qureg, c_qubits, c_state, target-1, M_j)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
        

        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end


        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_pauliX()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    H=Matrix([1.0 1.0;1.0 -1.0])/sqrt(2)
    for t=1:10
        
        numQubits = rand(1:12)
        qureg = QuEST.createQureg(numQubits, env)
        target = rand(0:numQubits-1)

        v = zeros(qreal, 2^numQubits)
        v[1]=1
        i_right=1.0*Matrix(I, 2^target, 2^target)
        res = numQubits-target-1
        i_left =1.0*Matrix(I, 2^res, 2^res)
        U=kron(i_left, X, i_right)
        v=U*v

        QuEST.pauliX(qureg, target)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)

        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_pauliY()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    H=Matrix([1.0 1.0;1.0 -1.0])/sqrt(2)
    for t=1:10
        
        numQubits = rand(1:12)
        qureg = QuEST.createQureg(numQubits, env)
        target = rand(0:numQubits-1)

        v = zeros(qreal, 2^numQubits)
        v[1]=1
        i_right=1.0*Matrix(I, 2^target, 2^target)
        res = numQubits-target-1
        i_left =1.0*Matrix(I, 2^res, 2^res)
        U=kron(i_left, Y, i_right)
        v=U*v

        QuEST.pauliY(qureg, target)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)

        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_pauliZ()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    H=Matrix([1.0 1.0;1.0 -1.0])/sqrt(2)
    for t=1:10
        
        numQubits = rand(1:12)
        qureg = QuEST.createQureg(numQubits, env)
        target = rand(0:numQubits-1)

        v = zeros(qreal, 2^numQubits)
        v[1]=1
        i_right=1.0*Matrix(I, 2^target, 2^target)
        res = numQubits-target-1
        i_left =1.0*Matrix(I, 2^res, 2^res)
        U=kron(i_left, Z, i_right)
        v=U*v

        QuEST.pauliZ(qureg, target)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)

        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_phaseShift()
    env= QuEST.createQuESTEnv()
    for t=1:10
        
        numQubits = rand(2:12)
        target = rand(0:numQubits-1)
        
        θ = rand(qreal)
        qureg = QuEST.createQureg(numQubits, env)

        
        QuEST.pauliX(qureg, target)
        QuEST.phaseShift(qureg, target, θ)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
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
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_rotateAroundAxis()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(2:12)
        target = rand(0:numQubits-1)
        
        θ = rand(qreal)
        n̂ = rand(qreal, 3)
        n̂ /= norm(n̂)
        U = ℯ^(-im*θ*(n̂[1]*X+n̂[2]*Y+n̂[3]*Z)/2)
        qureg = QuEST.createQureg(numQubits, env)

        v = [1.0, 0.0]
        v = U*v
        
        QuEST.rotateAroundAxis(qureg, target, θ, n̂)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
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
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_rotateX()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(2:12)
        target = rand(0:numQubits-1)
        control = rand(1:numQubits)

        θ = rand(qreal)
        
        U = ℯ^(-im*θ*X/2)
        qureg = QuEST.createQureg(numQubits, env)

        v = [1.0, 0.0]
        v = U*v
        
        
        QuEST.rotateX(qureg, target, θ)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
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
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_rotateY()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(2:12)
        target = rand(0:numQubits-1)
        control = rand(1:numQubits)

        θ = rand(qreal)
        
        U = ℯ^(-im*θ*Y/2)
        qureg = QuEST.createQureg(numQubits, env)

        v = [1.0, 0.0]
        v = U*v
        
        
        QuEST.rotateY(qureg, target, θ)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
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
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_rotateZ()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(2:12)
        target = rand(0:numQubits-1)
        control = rand(1:numQubits)

        θ = rand(qreal)
        
        U = ℯ^(-im*θ*Z/2)
        qureg = QuEST.createQureg(numQubits, env)

        v = [1.0, 0.0]
        v = U*v
        
        
        QuEST.rotateZ(qureg, target, θ)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
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
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_sGate()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    H=Matrix([1.0 1.0;1.0 -1.0])/sqrt(2)
    S=Matrix([1.0 0;0 im])
    for t=1:10
        
        numQubits = rand(1:12)
        qureg = QuEST.createQureg(numQubits, env)
        target = rand(0:numQubits-1)

        v = zeros(qreal, 2^numQubits)
        v[1]=1
        i_right=1.0*Matrix(I, 2^target, 2^target)
        res = numQubits-target-1
        i_left =1.0*Matrix(I, 2^res, 2^res)
        U=kron(i_left, S, i_right)
        v=U*v

        QuEST.sGate(qureg, target)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)

        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_sqrtSwapGate()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    H=Matrix([1.0 1.0;1.0 -1.0])/sqrt(2)
    S=Matrix([1.0 0;0 im])
    Swap=Matrix([1 0 0 0;0 0 1 0;0 1 0 0;0 0 0 1])
    SqrtSwap=sqrt(Swap)
    for t=1:10
        
        numQubits = rand(3:12)
        qureg = QuEST.createQureg(numQubits, env)
        qubit1 = rand(0:numQubits-2)

        v = zeros(qreal, 2^numQubits)
        v[1]=1
        i_right=1.0*Matrix(I, 2^qubit1, 2^qubit1)
        res = numQubits-qubit1-2
        i_left =1.0*Matrix(I, 2^res, 2^res)
        U=kron(i_left, SqrtSwap, i_right)
        v=U*v

        QuEST.sqrtSwapGate(qureg, qubit1, qubit1+1)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)

        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_swapGate()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    H=Matrix([1.0 1.0;1.0 -1.0])/sqrt(2)
    S=Matrix([1.0 0;0 im])
    Swap=Matrix([1 0 0 0;0 0 1 0;0 1 0 0;0 0 0 1])
    
    for t=1:10
        
        numQubits = rand(3:12)
        qureg = QuEST.createQureg(numQubits, env)
        qubit1 = rand(0:numQubits-2)

        v = zeros(qreal, 2^numQubits)
        v[1]=1
        i_right=1.0*Matrix(I, 2^qubit1, 2^qubit1)
        res = numQubits-qubit1-2
        i_left =1.0*Matrix(I, 2^res, 2^res)
        U=kron(i_left, Swap, i_right)
        v=U*v

        QuEST.swapGate(qureg, qubit1, qubit1+1)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)

        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_tGate()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    H=Matrix([1.0 1.0;1.0 -1.0])/sqrt(2)
    T=Matrix([1.0 0;0 ℯ^(im*π/4)])
    for t=1:10
        
        numQubits = rand(1:12)
        qureg = QuEST.createQureg(numQubits, env)
        target = rand(0:numQubits-1)

        v = zeros(qreal, 2^numQubits)
        v[1]=1
        i_right=1.0*Matrix(I, 2^target, 2^target)
        res = numQubits-target-1
        i_left =1.0*Matrix(I, 2^res, 2^res)
        U=kron(i_left, T, i_right)
        v=U*v

        QuEST.tGate(qureg, target)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)

        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_twoQubitUnitary()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        
        M_j = rand(Haar(2), 4)
        if qreal == Float32
            M_j = Matrix{Complex{qreal}}(M_j)
        end

        numQubits = rand(3:12)
        qureg = QuEST.createQureg(numQubits, env)
        qubit1 = rand(0:numQubits-2)

        v = zeros(qreal, 2^numQubits)
        v[1]=1.0
        i_right=1.0*Matrix(I, 2^qubit1, 2^qubit1)
        res = numQubits-qubit1-2
        i_left =1.0*Matrix(I, 2^res, 2^res)
        U=kron(i_left, M_j, i_right)
        v=U*v
        
        
        QuEST.twoQubitUnitary(qureg, qubit1, qubit1+1, M_j)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
        
        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end


function test_unitary()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(2:12)
        qureg = QuEST.createQureg(numQubits, env)
        target = rand(0:numQubits-1)
        
        M_j = rand(Haar(2), 2)
        if qreal == Float32
            M_j = Matrix{Complex{qreal}}(M_j)
        end
        v = [1.0, 0.0]
        v = M_j*v
        
        
        QuEST.unitary(qureg, target, M_j)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
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
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end



function test_collapseToOutcome()
    env= QuEST.createQuESTEnv()
    for t=1:100
        numQubits = rand(2:12)
        θ=rand(qreal)*π
        qureg = QuEST.createQureg(numQubits, env)
        
        target=rand(0:numQubits-1)

        for ind =0:numQubits-1
            if ind!=target
                QuEST.rotateX(qureg, ind, rand(qreal)*π)
            end
        end

        QuEST.rotateX(qureg, target, θ)

        p = QuEST.collapseToOutcome(qureg, target, 0)

        @test p ≈ (1 + cos(θ))/2 atol = tolerance

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)


    end
    
    QuEST.destroyQuESTEnv(env)

end

function test_measure()
    env= QuEST.createQuESTEnv()
    for t=1:100
        numQubits = rand(2:12)
        θ=rand(qreal)*π
        qureg = QuEST.createQureg(numQubits, env)
        
        target=rand(0:numQubits-1)

        for ind =0:numQubits-1
            if ind!=target
                QuEST.rotateX(qureg, ind, rand(qreal)*π)
            end
        end

        QuEST.rotateX(qureg, target, θ)

        p = QuEST.measure(qureg, target)

        @test p in [0, 1]


    end
    
    QuEST.destroyQuESTEnv(env)

end

function test_measureWithStats()
    env= QuEST.createQuESTEnv()
    for t=1:100
        numQubits = 1
        θ=rand(qreal)*π
        qureg = QuEST.createQureg(numQubits, env)        

        QuEST.rotateX(qureg, 0, θ)

        cbit, prob = QuEST.measureWithStats(qureg, 0)

        @test cbit in [0, 1]
        if cbit == 0 
            @test prob ≈ (1 + cos(θ))/2 atol = tolerance
        else
            @test prob ≈ (1 - cos(θ))/2 atol = tolerance
        end

    end
    
    QuEST.destroyQuESTEnv(env)

end

function test_DiagonalOp()
    env= QuEST.createQuESTEnv()
    for i=1:100
        num_qubits = rand(1:12)
        qureg = QuEST.createQureg(num_qubits, env)

        op = QuEST.createDiagonalOp(num_qubits, env)
        op_j = rand(Complex{qreal}, 2^num_qubits)
        QuEST.initDiagonalOp(op, real(op_j), imag(op_j))

        codes=rand(0:1, num_qubits)
        #codes=zeros(num_qubits)
        basis = 1
        for ind=0:num_qubits-1
            if codes[ind+1] == 1
                QuEST.pauliX(qureg, ind)
                basis += 2^ind
            end
        end
        QuEST.applyDiagonalOp(qureg, op)
        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^num_qubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^num_qubits)
    
        @test real(op_j[basis]) ≈ reals[basis] atol = tolerance
        @test imag(op_j[basis]) ≈ imags[basis] atol = tolerance

        QuEST.destroyDiagonalOp(op, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_applyMatrix2()
    env= QuEST.createQuESTEnv()
    for t=1:10
        
        numQubits = rand(2:12)
        qureg = QuEST.createQureg(numQubits, env)
        target = rand(0:numQubits-1)
        
        M_j = rand(Haar(2), 2)
        if qreal == Float32
            M_j = Matrix{Complex{qreal}}(M_j)
        end
        v = [1.0, 0.0]
        v = M_j*v
        
        M = QuEST._quest_mtx_2(M_j)
        QuEST.applyMatrix2(qureg, target, M)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
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
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_applyMatrix4()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        
        M_j = rand(Haar(2), 4)
        if qreal == Float32
            M_j = Matrix{Complex{qreal}}(M_j)
        end

        M = QuEST._quest_mtx_4(M_j)

        numQubits = rand(3:12)
        qureg = QuEST.createQureg(numQubits, env)
        qubit1 = rand(0:numQubits-2)

        v = zeros(qreal, 2^numQubits)
        v[1]=1.0
        i_right=1.0*Matrix(I, 2^qubit1, 2^qubit1)
        res = numQubits-qubit1-2
        i_left =1.0*Matrix(I, 2^res, 2^res)
        U=kron(i_left, M_j, i_right)
        v=U*v
        
        
        QuEST.applyMatrix4(qureg, qubit1, qubit1+1, M)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
        
        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_applyMatrixN()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    for t=1:10
        
        numQubits = rand(3:12)
        qureg = QuEST.createQureg(numQubits, env)
        num_targs = rand(1:numQubits-1)

        targs = Cint.([x for x in 0:num_targs-1])

        M_j = rand(Haar(2), 2^num_targs)
        if qreal == Float32
            M_j = Matrix{Complex{qreal}}(M_j)
        end
        M = QuEST.make_QuEST_matrix(M_j)

        v = zeros(qreal, 2^numQubits)
        v[1]=1.0

        res = numQubits-num_targs
        Idm = 1.0*Matrix(I,2^res, 2^res)
        U=kron(Idm, M_j)
        v=U*v
        
        
        QuEST.applyMatrixN(qureg, targs, M)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
        

        for ind =1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end


        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_applyMultiControlledMatrixN()
    print("warning check test_multiControlledMultiQubitUnitary\n")
    env= QuEST.createQuESTEnv()
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
        if qreal == Float32
            M_j = Matrix{Complex{qreal}}(M_j)
        end
        M = QuEST.make_QuEST_matrix(M_j)
    
        qureg = QuEST.createQureg(numQubits, env)
        
        for ind in ctrls
            QuEST.pauliX(qureg, ind)
        end
        #QuEST.pauliX(qureg, 0)
        QuEST.applyMultiControlledMatrixN(qureg, ctrls, targs, M)
        
        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)

        X=Matrix([0 1.0;1 0])
        big_X=Matrix([0 1.0;1 0])
        for iter = 1:num_ctrls-1
            big_X=kron(big_X, X)
        end
        res = numQubits-num_targs-num_ctrls
        Idm = 1.0*Matrix(I,2^res, 2^res)
        U=kron(Idm, M_j, big_X)
        v = zeros(qreal, 2^numQubits)
        v[1]=1.0
        v=U*v
        #print(reals)
        #print(imags)
        #print(v)
        for ind = 1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end
        QuEST.destroyComplexMatrixN(M)
        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
    
end


function test_applyPauliHamil()
    env= QuEST.createQuESTEnv()

    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    Id=Matrix([1.0 0;0 1])

    paulies = [QuEST.QuEST_Types.PAULI_I, 
               QuEST.QuEST_Types.PAULI_X, 
               QuEST.QuEST_Types.PAULI_Y, 
               QuEST.QuEST_Types.PAULI_Z]

    pauli_dict=Dict()
    pauli_dict[QuEST.QuEST_Types.PAULI_I]=Id
    pauli_dict[QuEST.QuEST_Types.PAULI_X]=X
    pauli_dict[QuEST.QuEST_Types.PAULI_Y]=Y
    pauli_dict[QuEST.QuEST_Types.PAULI_Z]=Z
    
    for i=1:10
        numQubits = rand(1:12)
        numSumTerms = rand(1:20)
        hamil = QuEST.createPauliHamil(numQubits, numSumTerms)


        codes = rand(paulies, numQubits*numSumTerms)
        coeffs = rand(qreal, numSumTerms)
        QuEST.initPauliHamil(hamil, coeffs, codes)
        
        H=zeros(2^numQubits, 2^numQubits)
        for term=0:numSumTerms-1
            h_tmp=ones(1, 1)
            for qubit =0:numQubits-1
                h_tmp = kron(pauli_dict[codes[term*numQubits+qubit+1]], h_tmp)
            end
            H += coeffs[term+1]*h_tmp
        end

        v = zeros(qreal, 2^numQubits)
        v[1]=1.0
        v=H*v

        in_qureg = QuEST.createQureg(numQubits, env)
        out_qureg = QuEST.createQureg(numQubits, env)

        QuEST.applyPauliHamil(in_qureg, hamil, out_qureg)

        reals = unsafe_wrap(Vector{qreal}, out_qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, out_qureg.stateVec.imag, 2^numQubits)

        for ind = 1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end

        QuEST.destroyQureg(in_qureg, env)
        QuEST.destroyQureg(out_qureg, env)
        QuEST.destroyPauliHamil(hamil)


    end
    QuEST.destroyQuESTEnv(env)
end

function test_applyPauliSum()
    env= QuEST.createQuESTEnv()

    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    Id=Matrix([1.0 0;0 1])

    paulies = [QuEST.QuEST_Types.PAULI_I, 
               QuEST.QuEST_Types.PAULI_X, 
               QuEST.QuEST_Types.PAULI_Y, 
               QuEST.QuEST_Types.PAULI_Z]

    pauli_dict=Dict()
    pauli_dict[QuEST.QuEST_Types.PAULI_I]=Id
    pauli_dict[QuEST.QuEST_Types.PAULI_X]=X
    pauli_dict[QuEST.QuEST_Types.PAULI_Y]=Y
    pauli_dict[QuEST.QuEST_Types.PAULI_Z]=Z
    
    for i=1:10
        numQubits = rand(1:12)
        numSumTerms = rand(1:20)


        codes = rand(paulies, numQubits*numSumTerms)
        coeffs = rand(qreal, numSumTerms)
        
        H=zeros(2^numQubits, 2^numQubits)
        for term=0:numSumTerms-1
            h_tmp=ones(1, 1)
            for qubit =0:numQubits-1
                h_tmp = kron(pauli_dict[codes[term*numQubits+qubit+1]], h_tmp)
            end
            H += coeffs[term+1]*h_tmp
        end

        v = zeros(qreal, 2^numQubits)
        v[1]=1.0
        v=H*v

        in_qureg = QuEST.createQureg(numQubits, env)
        out_qureg = QuEST.createQureg(numQubits, env)

        QuEST.applyPauliSum(in_qureg, codes, coeffs, numSumTerms, out_qureg)

        reals = unsafe_wrap(Vector{qreal}, out_qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, out_qureg.stateVec.imag, 2^numQubits)

        for ind = 1:2^numQubits
            @test reals[ind] ≈ real(v[ind]) atol = tolerance
            @test imags[ind] ≈ imag(v[ind]) atol = tolerance
        end

        QuEST.destroyQureg(in_qureg, env)
        QuEST.destroyQureg(out_qureg, env)


    end
    QuEST.destroyQuESTEnv(env)
end

function test_applyTrotterCircuit()
    env= QuEST.createQuESTEnv()

    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    Id=Matrix([1.0 0;0 1])

    paulies = [QuEST.QuEST_Types.PAULI_I, 
               QuEST.QuEST_Types.PAULI_X, 
               #QuEST.QuEST_Types.PAULI_Y, 
               QuEST.QuEST_Types.PAULI_Z]

    pauli_dict=Dict()
    pauli_dict[QuEST.QuEST_Types.PAULI_I]=Id
    pauli_dict[QuEST.QuEST_Types.PAULI_X]=X
    pauli_dict[QuEST.QuEST_Types.PAULI_Y]=Y
    pauli_dict[QuEST.QuEST_Types.PAULI_Z]=Z
    
    for i=1:10
        numQubits = rand(1:12)
        numSumTerms = rand(1:4)
        hamil = QuEST.createPauliHamil(numQubits, numSumTerms)


        codes = rand(paulies, numQubits*numSumTerms)
        coeffs = rand(qreal, numSumTerms)
        reps = rand(1:10)
        t=rand(qreal)

        QuEST.initPauliHamil(hamil, coeffs, codes)
        qureg = QuEST.createQureg(numQubits, env)
        test_qureg = QuEST.createQureg(numQubits, env)
        
        #QuEST.startRecordingQASM(test_qureg)
        
        for trot=1:reps
            for term=0:numSumTerms-1
                z_list=Vector{Cint}()
                for qubit =0:numQubits-1
                    if codes[term*numQubits+qubit+1] == QuEST.QuEST_Types.PAULI_X
                        QuEST.hadamard(test_qureg, qubit)
                    end
                end
                for qubit =0:numQubits-1
                    if codes[term*numQubits+qubit+1] != QuEST.QuEST_Types.PAULI_I
                        push!(z_list, qubit)
                    end
                end

                for ind=1:length(z_list)-1
                    QuEST.controlledNot(test_qureg, z_list[ind], z_list[ind+1])
                end

                if length(z_list) !=0
                    QuEST.rotateZ(test_qureg, z_list[length(z_list)], 2*t*coeffs[term+1]/reps)
                end

                for ind=length(z_list)-1:-1:1
                    QuEST.controlledNot(test_qureg, z_list[ind], z_list[ind+1])
                end

                for qubit =0:numQubits-1
                    if codes[term*numQubits+qubit+1] == QuEST.QuEST_Types.PAULI_X
                        QuEST.hadamard(test_qureg, qubit)
                    end
                end
            end

        end
            
        #QuEST.stopRecordingQASM(test_qureg)
        #QuEST.printRecordedQASM(test_qureg)
        #print("##################################################")
        #QuEST.startRecordingQASM(qureg)
        QuEST.applyTrotterCircuit(qureg, hamil, t, 1, reps)
        #QuEST.stopRecordingQASM(qureg)
        #QuEST.printRecordedQASM(qureg)

        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)

        reals2 = unsafe_wrap(Vector{qreal}, test_qureg.stateVec.real, 2^numQubits)
        imags2 = unsafe_wrap(Vector{qreal}, test_qureg.stateVec.imag, 2^numQubits)

        # for ind = 1:2^numQubits
        #    @test isapprox(reals[ind], real(v[ind]);atol = 1e-2)
        #    @test isapprox(imags[ind], imag(v[ind]); atol = 1e-2)
        # end

        for ind = 1:2^numQubits
            @test reals[ind] ≈ reals2[ind] atol = tolerance
            @test imags[ind] ≈ imags2[ind] atol = tolerance
        end

        QuEST.destroyQureg(qureg, env)
        QuEST.destroyPauliHamil(hamil)


    end
    QuEST.destroyQuESTEnv(env)
end


function test_calcDensityInnerProduct()
    env= QuEST.createQuESTEnv()
    
    for i=1:20
        numQubits = rand(1:5)
        
        targs = Cint.([x for x in 0:numQubits-1])

        M1 = rand(Haar(2), 2^numQubits)
        if qreal == Float32
            M1 = Matrix{Complex{qreal}}(M1)
        end
        U1 = QuEST.make_QuEST_matrix(M1)

        M2 = rand(Haar(2), 2^numQubits)
        if qreal == Float32
            M2 = Matrix{Complex{qreal}}(M2)
        end
        U2 = QuEST.make_QuEST_matrix(M2)

        v = zeros(qreal, 2^numQubits)
        v[1]=1.0

        state1 = M1*v
        state2 = M2*v

        my_inner_prod = transpose(state1)*conj(state2)
        my_inner_prod = my_inner_prod*conj(my_inner_prod)

        qureg1 = QuEST.createDensityQureg(numQubits, env)
        qureg2 = QuEST.createDensityQureg(numQubits, env)
        
        QuEST.multiQubitUnitary(qureg1, targs, U1)
        QuEST.multiQubitUnitary(qureg2, targs, U2)

        @test my_inner_prod ≈ QuEST.calcDensityInnerProduct(qureg1, qureg2) atol = tolerance

        QuEST.destroyQureg(qureg1, env)
        QuEST.destroyQureg(qureg2, env)

    end
    QuEST.destroyQuESTEnv(env)
end


function test_calcExpecDiagonalOp()
    env= QuEST.createQuESTEnv()
    
    for i=1:20
        numQubits = rand(1:12)
        
        qureg = QuEST.createQureg(numQubits, env)

        state = rand(Complex{qreal}, 2^numQubits)
        state /= norm(state)

        QuEST.initStateFromAmps(qureg, real.(state), imag.(state))
        

        op = QuEST.createDiagonalOp(numQubits, env)
        op_j = rand(Complex{qreal}, 2^numQubits)
        QuEST.initDiagonalOp(op, real(op_j), imag(op_j))
        
        rhs = QuEST.calcExpecDiagonalOp(qureg, op)

        lhs = sum([state[x]*conj(state[x])*op_j[x] for x =1:2^numQubits])

        @test rhs.real ≈ real(lhs) atol = tolerance
        @test rhs.imag ≈ imag(lhs) atol = tolerance

        QuEST.destroyDiagonalOp(op, env)

        
        QuEST.destroyQureg(qureg, env)

    end
    QuEST.destroyQuESTEnv(env)
end

function test_calcExpecPauliHamil()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    Id=Matrix([1.0 0;0 1])

    paulies = [QuEST.QuEST_Types.PAULI_I, 
               QuEST.QuEST_Types.PAULI_X, 
               QuEST.QuEST_Types.PAULI_Y, 
               QuEST.QuEST_Types.PAULI_Z]

    pauli_dict=Dict()
    pauli_dict[QuEST.QuEST_Types.PAULI_I]=Id
    pauli_dict[QuEST.QuEST_Types.PAULI_X]=X
    pauli_dict[QuEST.QuEST_Types.PAULI_Y]=Y
    pauli_dict[QuEST.QuEST_Types.PAULI_Z]=Z
    
    for i=1:20
        numQubits = rand(1:12)
        numSumTerms = rand(1:20)

        hamil = QuEST.createPauliHamil(numQubits, numSumTerms)
        codes = rand(paulies, numQubits*numSumTerms)
        coeffs = rand(qreal, numSumTerms)
        QuEST.initPauliHamil(hamil, coeffs, codes)

        qureg = QuEST.createQureg(numQubits, env)
        work = QuEST.createQureg(numQubits, env)

        H=zeros(2^numQubits, 2^numQubits)
        for term=0:numSumTerms-1
            h_tmp=ones(1, 1)
            for qubit =0:numQubits-1
                h_tmp = kron(pauli_dict[codes[term*numQubits+qubit+1]], h_tmp)
            end
            H += coeffs[term+1]*h_tmp
        end

        state = rand(Complex{qreal}, 2^numQubits)
        state /= norm(state)

        QuEST.initStateFromAmps(qureg, real.(state), imag.(state))
        

        lhs = transpose(conj(state))*H*state

        rhs = QuEST.calcExpecPauliHamil(qureg, hamil, work)

        @test rhs ≈ real(lhs) atol = tolerance
        


        QuEST.destroyQureg(qureg, env)
        QuEST.destroyQureg(work, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_calcExpecPauliProd()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    Id=Matrix([1.0 0;0 1])

    paulies = [QuEST.QuEST_Types.PAULI_I, 
               QuEST.QuEST_Types.PAULI_X, 
               QuEST.QuEST_Types.PAULI_Y, 
               QuEST.QuEST_Types.PAULI_Z]

    pauli_dict=Dict()
    pauli_dict[QuEST.QuEST_Types.PAULI_I]=Id
    pauli_dict[QuEST.QuEST_Types.PAULI_X]=X
    pauli_dict[QuEST.QuEST_Types.PAULI_Y]=Y
    pauli_dict[QuEST.QuEST_Types.PAULI_Z]=Z
    
    for i=1:20
        numQubits = rand(1:12)

        codes = rand(paulies, numQubits)

        qureg = QuEST.createQureg(numQubits, env)
        work = QuEST.createQureg(numQubits, env)

        H=zeros(2^numQubits, 2^numQubits)
        
        H=ones(1, 1)
        c_targs = Vector{Cint}([x for x in 0:numQubits-1])
        #c_codes = Vector{QuEST.QuEST_Types.pauliOpType}()

        H=zeros(2^numQubits, 2^numQubits)
        H=ones(1, 1)
        h_tmp=ones(1, 1)
        for qubit =0:numQubits-1
            H = kron(pauli_dict[codes[qubit+1]], H)
        end


        state = rand(Complex{qreal}, 2^numQubits)
        state /= norm(state)

        QuEST.initStateFromAmps(qureg, real.(state), imag.(state))
        

        lhs = transpose(conj(state))*H*state

        rhs = QuEST.calcExpecPauliProd(qureg, c_targs, codes, work)

        @test rhs ≈ real(lhs) atol = tolerance
        


        QuEST.destroyQureg(qureg, env)
        QuEST.destroyQureg(work, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_calcExpecPauliSum()
    env= QuEST.createQuESTEnv()
    X=Matrix([0 1.0;1 0])
    Y=Matrix([0 -im;im 0])
    Z=Matrix([1.0 0;0 -1.0])
    Id=Matrix([1.0 0;0 1])

    paulies = [QuEST.QuEST_Types.PAULI_I, 
               QuEST.QuEST_Types.PAULI_X, 
               QuEST.QuEST_Types.PAULI_Y, 
               QuEST.QuEST_Types.PAULI_Z]

    pauli_dict=Dict()
    pauli_dict[QuEST.QuEST_Types.PAULI_I]=Id
    pauli_dict[QuEST.QuEST_Types.PAULI_X]=X
    pauli_dict[QuEST.QuEST_Types.PAULI_Y]=Y
    pauli_dict[QuEST.QuEST_Types.PAULI_Z]=Z
    
    for i=1:20
        numQubits = rand(1:12)
        numSumTerms = rand(1:20)

        codes = rand(paulies, numQubits*numSumTerms)
        coeffs = rand(qreal, numSumTerms)

        qureg = QuEST.createQureg(numQubits, env)
        work = QuEST.createQureg(numQubits, env)

        H=zeros(2^numQubits, 2^numQubits)
        for term=0:numSumTerms-1
            h_tmp=ones(1, 1)
            for qubit =0:numQubits-1
                h_tmp = kron(pauli_dict[codes[term*numQubits+qubit+1]], h_tmp)
            end
            H += coeffs[term+1]*h_tmp
        end

        state = rand(Complex{qreal}, 2^numQubits)
        state /= norm(state)

        QuEST.initStateFromAmps(qureg, real.(state), imag.(state))
        

        lhs = transpose(conj(state))*H*state

        rhs = QuEST.calcExpecPauliSum(qureg, codes, coeffs, work)

        @test isapprox(rhs, real(lhs); atol=1e-5)
        


        QuEST.destroyQureg(qureg, env)
        QuEST.destroyQureg(work, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_calcFidelity()
    env= QuEST.createQuESTEnv()
    
    for i=1:20
        numQubits = rand(1:12)


        qureg1 = QuEST.createQureg(numQubits, env)
        qureg2 = QuEST.createQureg(numQubits, env)

        
        state1 = rand(Complex{qreal}, 2^numQubits)
        state1 /= norm(state1)

        state2 = rand(Complex{qreal}, 2^numQubits)
        state2 /= norm(state2)

        QuEST.initStateFromAmps(qureg1, real.(state1), imag.(state1))
        QuEST.initStateFromAmps(qureg2, real.(state2), imag.(state2))

        lhs = abs(transpose(conj(state1))*state2)^2

        rhs = QuEST.calcFidelity(qureg1, qureg2)

        @test isapprox(rhs, real(lhs); atol=1e-5)
        


        QuEST.destroyQureg(qureg1, env)
        QuEST.destroyQureg(qureg2, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_calcHilbertSchmidtDistance()
    env= QuEST.createQuESTEnv()
    
    for i=1:20
        numQubits = rand(1:5)
        
        targs = Cint.([x for x in 0:numQubits-1])

        M1 = rand(Haar(2), 2^numQubits)
        if qreal == Float32
            M1 = Matrix{Complex{qreal}}(M1)
        end
        U1 = QuEST.make_QuEST_matrix(M1)

        M2 = rand(Haar(2), 2^numQubits)
        if qreal == Float32
            M2 = Matrix{Complex{qreal}}(M2)
        end
        U2 = QuEST.make_QuEST_matrix(M2)

        v = zeros(qreal, 2^numQubits)
        v[1]=1.0

        state1 = M1*v
        state2 = M2*v

        my_inner_prod = sqrt(sum([abs( state1[i]*conj(state1[j]) - state2[i]*conj(state2[j]))^2 for i=1:2^numQubits for j=1:2^numQubits]))

        qureg1 = QuEST.createDensityQureg(numQubits, env)
        qureg2 = QuEST.createDensityQureg(numQubits, env)
        
        QuEST.multiQubitUnitary(qureg1, targs, U1)
        QuEST.multiQubitUnitary(qureg2, targs, U2)

        @test real(my_inner_prod) ≈ QuEST.calcHilbertSchmidtDistance(qureg1, qureg2) atol = tolerance

        QuEST.destroyQureg(qureg1, env)
        QuEST.destroyQureg(qureg2, env)

    end
    QuEST.destroyQuESTEnv(env)
end

function test_calcInnerProduct()
    env= QuEST.createQuESTEnv()
    
    for i=1:20
        numQubits = rand(1:5)
        
        targs = Cint.([x for x in 0:numQubits-1])

        M1 = rand(Haar(2), 2^numQubits)
        if qreal == Float32
            M1 = Matrix{Complex{qreal}}(M1)
        end
        U1 = QuEST.make_QuEST_matrix(M1)

        M2 = rand(Haar(2), 2^numQubits)
        if qreal == Float32
            M2 = Matrix{Complex{qreal}}(M2)
        end
        U2 = QuEST.make_QuEST_matrix(M2)

        v = zeros(qreal, 2^numQubits)
        v[1]=1.0

        state1 = M1*v
        state2 = M2*v

        my_inner_prod = transpose(conj(state1))*state2

        qureg1 = QuEST.createQureg(numQubits, env)
        qureg2 = QuEST.createQureg(numQubits, env)
        
        QuEST.multiQubitUnitary(qureg1, targs, U1)
        QuEST.multiQubitUnitary(qureg2, targs, U2)

        @test my_inner_prod ≈ QuEST.calcInnerProduct(qureg1, qureg2) atol = tolerance

        QuEST.destroyQureg(qureg1, env)
        QuEST.destroyQureg(qureg2, env)

    end
    QuEST.destroyQuESTEnv(env)
end

function test_calcProbOfOutcome()
    env= QuEST.createQuESTEnv()
    zero_zero = Matrix([1.0 0;0 0])
    one_one = Matrix([0 0;0 1.0])
    for i=1:20
        numQubits = rand(1:12)
        
        target = rand(0:numQubits-1)

        qureg = QuEST.createQureg(numQubits, env)
        
        state = rand(Complex{qreal}, 2^numQubits)
        state /= norm(state)

        QuEST.initStateFromAmps(qureg, real.(state), imag.(state))
        
        Id_right = 1.0*Matrix(I,2^target, 2^target)
        Id_left = 1.0*Matrix(I,2^(numQubits - target -1), 2^(numQubits - target -1))
        
        M_zero = kron(Id_left, zero_zero, Id_right)
        M_one = kron(Id_left, one_one, Id_right)

        prob_zero = transpose(conj(state))*M_zero*state
        prob_one = transpose(conj(state))*M_one*state

        @test real(prob_zero) ≈ QuEST.calcProbOfOutcome(qureg, target, 0) atol = tolerance
        @test real(prob_one) ≈ QuEST.calcProbOfOutcome(qureg, target, 1) atol = tolerance

        QuEST.destroyQureg(qureg, env)

    end
    QuEST.destroyQuESTEnv(env)
end

function test_calcPurity()
    env= QuEST.createQuESTEnv()
    
    for i=1:20
        numQubits = rand(1:5)
        
        targs = Cint.([x for x in 0:numQubits-1])

        M1 = rand(Haar(2), 2^numQubits)
        if qreal == Float32
            M1 = Matrix{Complex{qreal}}(M1)
        end
        U1 = QuEST.make_QuEST_matrix(M1)

        qureg1 = QuEST.createDensityQureg(numQubits, env)
        
        QuEST.multiQubitUnitary(qureg1, targs, U1)

        @test 1 ≈ QuEST.calcPurity(qureg1) atol = tolerance

        QuEST.destroyQureg(qureg1, env)

    end
    QuEST.destroyQuESTEnv(env)
end

function test_calcTotalProb()
    env= QuEST.createQuESTEnv()
    
    for i=1:20
        numQubits = rand(1:5)
        
        targs = Cint.([x for x in 0:numQubits-1])

        M1 = rand(Haar(2), 2^numQubits)
        if qreal == Float32
            M1 = Matrix{Complex{qreal}}(M1)
        end
        U1 = QuEST.make_QuEST_matrix(M1)

        qureg1 = QuEST.createQureg(numQubits, env)
        
        QuEST.multiQubitUnitary(qureg1, targs, U1)

        @test 1 ≈ QuEST.calcTotalProb(qureg1) atol = tolerance

        QuEST.destroyQureg(qureg1, env)

    end
    QuEST.destroyQuESTEnv(env)
end

function test_amp_funcs()
    env= QuEST.createQuESTEnv()
    for i=1:20
        numQubits = rand(1:12)
        
        qureg = QuEST.createQureg(numQubits, env)
        
        state = rand(Complex{qreal}, 2^numQubits)
        state /= norm(state)

        QuEST.initStateFromAmps(qureg, real.(state), imag.(state))
        
        for ind =0:2^numQubits-1
            @test QuEST.getAmp(qureg, ind) ≈ state[ind+1] atol = tolerance
            @test QuEST.getImagAmp(qureg, ind) ≈ imag(state[ind+1]) atol = tolerance
            @test QuEST.getRealAmp(qureg, ind) ≈ real(state[ind+1]) atol = tolerance
            @test QuEST.getProbAmp(qureg, ind) ≈ abs(state[ind+1])^2 atol = tolerance
            @test QuEST.getNumQubits(qureg) == numQubits
            @test QuEST.getNumAmps(qureg) == 2^numQubits
        end
        QuEST.destroyQureg(qureg, env)

    end
    QuEST.destroyQuESTEnv(env)
end

function test_getDensityAmp()
    env= QuEST.createQuESTEnv()
    
    for i=1:20
        numQubits = rand(1:5)
        
        targs = Cint.([x for x in 0:numQubits-1])

        M1 = rand(Haar(2), 2^numQubits)
        if qreal == Float32
            M1 = Matrix{Complex{qreal}}(M1)
        end
        U1 = QuEST.make_QuEST_matrix(M1)

        v = zeros(qreal, 2^numQubits)
        v[1]=1.0

        state1 = M1*v

        qureg1 = QuEST.createDensityQureg(numQubits, env)
        
        QuEST.multiQubitUnitary(qureg1, targs, U1)

        for i=0:2^numQubits-1
            for j=0:2^numQubits-1
                @test QuEST.getDensityAmp(qureg1, i, j) ≈ state1[i+1]*conj(state1[j+1]) atol = tolerance
            end
        end
        QuEST.destroyQureg(qureg1, env)
    end
    QuEST.destroyQuESTEnv(env)
end

function test_cloneQureg()
    env= QuEST.createQuESTEnv()
    for i=1:20
        numQubits = rand(1:12)
        
        qureg1 = QuEST.createQureg(numQubits, env)
        qureg2 = QuEST.createQureg(numQubits, env)
        
        state = rand(Complex{qreal}, 2^numQubits)
        state /= norm(state)

        QuEST.initStateFromAmps(qureg1, real.(state), imag.(state))

        QuEST.cloneQureg(qureg2, qureg1)
        
        
        for ind =0:2^numQubits-1
            @test QuEST.getAmp(qureg1, ind) ≈ QuEST.getAmp(qureg2, ind) atol = tolerance
        end
        QuEST.destroyQureg(qureg1, env)
        QuEST.destroyQureg(qureg2, env)

    end
    QuEST.destroyQuESTEnv(env)
end

function test_initBlankState()
    env= QuEST.createQuESTEnv()
    for i=1:20
        numQubits = rand(1:12)
        
        qureg = QuEST.createQureg(numQubits, env)
        
        QuEST.initBlankState(qureg)
        
        
        for ind =0:2^numQubits-1
            @test QuEST.getAmp(qureg, ind) ≈ 0 atol = tolerance
        end
        QuEST.destroyQureg(qureg, env)

    end
    QuEST.destroyQuESTEnv(env)
end

function test_initClassicalState()
    env= QuEST.createQuESTEnv()
    for i=1:20
        numQubits = rand(1:12)
        
        qureg = QuEST.createQureg(numQubits, env)
        state = rand(0:2^numQubits-1)

        QuEST.initClassicalState(qureg, state)
        
        
        for ind =0:2^numQubits-1
            if ind == state
                @test QuEST.getAmp(qureg, ind) ≈ 1 atol = tolerance
            else
                @test QuEST.getAmp(qureg, ind) ≈ 0 atol = tolerance
            end
        end
        QuEST.destroyQureg(qureg, env)

    end
    QuEST.destroyQuESTEnv(env)
end

function test_initPlusState()
    env= QuEST.createQuESTEnv()
    for i=1:20
        numQubits = rand(1:12)
        
        qureg = QuEST.createQureg(numQubits, env)

        QuEST.initPlusState(qureg)
        
        
        for ind =0:2^numQubits-1
            @test QuEST.getAmp(qureg, ind) ≈ 1/sqrt(2^numQubits) atol = tolerance
        end
        QuEST.destroyQureg(qureg, env)

    end
    QuEST.destroyQuESTEnv(env)
end

function test_initPureState()
    env= QuEST.createQuESTEnv()
    for i=1:20
        numQubits = rand(1:12)
        
        qureg1 = QuEST.createQureg(numQubits, env)
        qureg2 = QuEST.createQureg(numQubits, env)
        
        state = rand(Complex{qreal}, 2^numQubits)
        state /= norm(state)

        QuEST.initStateFromAmps(qureg1, real.(state), imag.(state))

        QuEST.initPureState(qureg2, qureg1)
        
        
        for ind =0:2^numQubits-1
            @test QuEST.getAmp(qureg1, ind) ≈ QuEST.getAmp(qureg2, ind) atol = tolerance
        end
        QuEST.destroyQureg(qureg1, env)
        QuEST.destroyQureg(qureg2, env)

    end
    QuEST.destroyQuESTEnv(env)
end

function test_initStateFromAmps()
    env= QuEST.createQuESTEnv()
    for i=1:20
        numQubits = rand(1:12)
        
        qureg1 = QuEST.createQureg(numQubits, env)
        
        state = rand(Complex{qreal}, 2^numQubits)
        state /= norm(state)

        QuEST.initStateFromAmps(qureg1, real.(state), imag.(state))
        
        for ind =0:2^numQubits-1
            @test QuEST.getAmp(qureg1, ind) ≈ state[ind+1] atol = tolerance
        end
        QuEST.destroyQureg(qureg1, env)

    end
    QuEST.destroyQuESTEnv(env)
end

function test_initZeroState()
    env= QuEST.createQuESTEnv()
    for i=1:20
        numQubits = rand(1:12)
        
        qureg1 = QuEST.createQureg(numQubits, env)
        
        

        QuEST.initZeroState(qureg1)
        
        @test QuEST.getAmp(qureg1, 0) ≈ 1 atol = tolerance

        for ind =1:2^numQubits-1
            @test QuEST.getAmp(qureg1, ind) ≈ 0 atol = tolerance
        end
        QuEST.destroyQureg(qureg1, env)

    end
    QuEST.destroyQuESTEnv(env)
end

function test_setAmps()
    env= QuEST.createQuESTEnv()
    for i=1:20
        numQubits = rand(1:12)
        
        qureg1 = QuEST.createQureg(numQubits, env)
        
        state = rand(Complex{qreal}, 2^numQubits)
        state /= norm(state)

        QuEST.initStateFromAmps(qureg1, real.(state), imag.(state))

        start_ind = rand(0:2^numQubits-2)
        finish_ind = rand(start_ind:2^numQubits-1)
        numAmps = finish_ind - start_ind +1

        new_amps = rand(Complex{qreal}, numAmps)

        QuEST.setAmps(qureg1, start_ind, real.(new_amps), imag.(new_amps), numAmps)
        
        for ind =0:2^numQubits-1
            if ind < start_ind
                @test QuEST.getAmp(qureg1, ind) ≈ state[ind+1] atol = tolerance
            elseif ind <= finish_ind
                @test QuEST.getAmp(qureg1, ind) ≈ new_amps[ind - start_ind+1] atol = tolerance
            else
                @test QuEST.getAmp(qureg1, ind) ≈ state[ind+1] atol = tolerance
            end
        end
        QuEST.destroyQureg(qureg1, env)

    end
    QuEST.destroyQuESTEnv(env)
end

function test_setWeightedQureg()
    env= QuEST.createQuESTEnv()
    for i=1:10
        numQubits = rand(1:12)
        
        qureg1 = QuEST.createQureg(numQubits, env)
        qureg2 = QuEST.createQureg(numQubits, env)
        qureg3 = QuEST.createQureg(numQubits, env)
        
        state1 = rand(Complex{qreal}, 2^numQubits)
        state1 /= norm(state1)

        state2 = rand(Complex{qreal}, 2^numQubits)
        state2 /= norm(state2)

        state3 = rand(Complex{qreal}, 2^numQubits)
        state3 /= norm(state3)

        QuEST.initStateFromAmps(qureg1, real.(state1), imag.(state1))
        QuEST.initStateFromAmps(qureg2, real.(state2), imag.(state2))
        QuEST.initStateFromAmps(qureg3, real.(state3), imag.(state3))

        fac1 = rand(Complex{qreal})
        fac2 = rand(Complex{qreal})
        fac3 = rand(Complex{qreal})

        QuEST.setWeightedQureg(fac1, qureg1, fac2, qureg2, fac3, qureg3)
        
        for ind =0:2^numQubits-1
            rhs  = fac1 * state1[ind+1]
            rhs += fac2 * state2[ind+1]
            rhs += fac3 * state3[ind+1]    
            @test QuEST.getAmp(qureg3, ind) ≈ rhs atol = tolerance
        end
        
        QuEST.destroyQureg(qureg1, env)
        QuEST.destroyQureg(qureg2, env)
        QuEST.destroyQureg(qureg3, env)

    end
    QuEST.destroyQuESTEnv(env)
end

function test_QASM()
    env= QuEST.createQuESTEnv()
    for i=1:10
        numQubits = rand(3:12)
        

        qureg = QuEST.createQureg(numQubits, env)
        QuEST.startRecordingQASM(qureg)
        qubits = rand(0:numQubits-1, 10)
        θ = qreal(rand(0.0:0.00001:1))
        α=rand(Complex{qreal})
        β=rand(Complex{qreal})
        norm = sqrt(α*conj(α) + β*conj(β))
        α /= norm
        β /= norm

        QuEST.hadamard(qureg, qubits[1])
        QuEST.tGate(qureg, qubits[2])
        QuEST.sGate(qureg, qubits[3])
        QuEST.pauliX(qureg, qubits[4])
        QuEST.pauliY(qureg, qubits[5])
        QuEST.pauliZ(qureg, qubits[6])
        QuEST.rotateX(qureg, qubits[7], θ)
        QuEST.rotateY(qureg, qubits[8], θ)
        QuEST.rotateZ(qureg, qubits[9], θ)

        original_stdout = stdout
        read_pipe, write_pipe = redirect_stdout()
        lines = Vector{String}()
        QuEST.printRecordedQASM(qureg)
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
        
        QuEST.stopRecordingQASM(qureg)
        QuEST.pauliX(qureg, qubits[1])

        original_stdout = stdout
        read_pipe, write_pipe = redirect_stdout()
        lines = Vector{String}()
        QuEST.printRecordedQASM(qureg)
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

        QuEST.startRecordingQASM(qureg)
        QuEST.pauliZ(qureg, qubits[1])
        QuEST.controlledNot(qureg, 0, 1)

        original_stdout = stdout
        read_pipe, write_pipe = redirect_stdout()
        lines = Vector{String}()
        QuEST.printRecordedQASM(qureg)
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
        
        QuEST.stopRecordingQASM(qureg)
        QuEST.clearRecordedQASM(qureg)
        QuEST.startRecordingQASM(qureg)

        QuEST.hadamard(qureg, qubits[1])
        QuEST.tGate(qureg, qubits[2])
        QuEST.sGate(qureg, qubits[3])

        original_stdout = stdout
        read_pipe, write_pipe = redirect_stdout()
        lines = Vector{String}()
        QuEST.printRecordedQASM(qureg)
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


        QuEST.writeRecordedQASMToFile(qureg, "test_QASM_file.txt")

        io = open("test_QASM_file.txt", "r")

        @test readline(io) == "h q[$(qubits[1])];"
        @test readline(io) == "t q[$(qubits[2])];"
        @test readline(io) == "s q[$(qubits[3])];"
        @test readline(io) == ""

        close(io)

        rm("test_QASM_file.txt")

        QuEST.destroyQureg(qureg, env)
    end
    QuEST.destroyQuESTEnv(env)
end
#

function test_GPU()
    env= QuEST.createQuESTEnv()
    for i=1:1
        numQubits = rand(3:12)
        qureg = QuEST.createQureg(numQubits, env)
        state = rand(Complex{qreal}, 2^numQubits)
        state /= norm(state)
        QuEST.initStateFromAmps(qureg, real.(state), imag.(state))
        QuEST.copyStateToGPU(qureg)
        QuEST.copyStateFromGPU(qureg)
        #env_str = QuEST.getEnvironmentString(env, qureg)
        #print(env_str)
        print("Add more stuff here to test GPU")

    end

    QuEST.destroyQuESTEnv(env)
end

function test_initDebugState()
    env= QuEST.createQuESTEnv()
    for i=1:10
        numQubits = rand(3:12)
        qureg = QuEST.createQureg(numQubits, env)
        QuEST.initDebugState(qureg)
        reals = unsafe_wrap(Vector{qreal}, qureg.stateVec.real, 2^numQubits)
        imags = unsafe_wrap(Vector{qreal}, qureg.stateVec.imag, 2^numQubits)
        for i=0:2^numQubits-1
            @test reals[i+1] ≈ 2i/10 atol = tolerance
            @test imags[i+1] ≈ (2i+1)/10 atol = tolerance
        end

    end

    QuEST.destroyQuESTEnv(env)
end

function test_reportPauliHamil()
    paulies = [QuEST.QuEST_Types.PAULI_I, 
               QuEST.QuEST_Types.PAULI_X, 
               QuEST.QuEST_Types.PAULI_Y, 
               QuEST.QuEST_Types.PAULI_Z]
    env= QuEST.createQuESTEnv()
    for i=1:10

        numQubits = rand(3:12)
        numSumTerms = rand(1:20)
        
        
        hamil = QuEST.createPauliHamil(numQubits, numSumTerms)
        codes = rand(paulies, numQubits*numSumTerms)
        coeffs = rand(qreal, numSumTerms)
        QuEST.initPauliHamil(hamil, coeffs, codes)


        original_stdout = stdout
        read_pipe, write_pipe = redirect_stdout()
        lines = Vector{String}()
        QuEST.reportPauliHamil(hamil)
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
            @test isapprox(parse(qreal, coeff_str), coeffs[row+1], atol=1e-5)
            for q = 1:numQubits
                @test paulies[parse(Cint, codes_list[q])+1] == codes[row*numQubits + q]
            end
        end

        @test lines[end] == ""

        

        QuEST.destroyPauliHamil(hamil)
    end

    QuEST.destroyQuESTEnv(env)
end

function test_reportQuregParams()

    env= QuEST.createQuESTEnv()
    numQubits = rand(3:12)
    qureg = QuEST.createQureg(numQubits, env)


    original_stdout = stdout
    read_pipe, write_pipe = redirect_stdout()
    lines = Vector{String}()

    QuEST.reportQuregParams(qureg)
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

    QuEST.destroyQureg(qureg, env)
    QuEST.destroyQuESTEnv(env)
end

function test_reporQuESTEnv()
    
    env= QuEST.createQuESTEnv()
    numQubits = rand(3:12)
    qureg = QuEST.createQureg(numQubits, env)

    original_stdout = stdout
    read_pipe, write_pipe = redirect_stdout()
    lines = Vector{String}()

    QuEST.reportQuESTEnv(env)
    Base.Libc.flush_cstdio()
    redirect_stdout(original_stdout)

    close(write_pipe)
    
    for l =1:7
        push!(lines, readline(read_pipe))
    end

    @test lines[1] == "EXECUTION ENVIRONMENT:"
    @test lines[end] == ""

    QuEST.destroyQureg(qureg, env)
    QuEST.destroyQuESTEnv(env)
end

function test_reportState()
    env= QuEST.createQuESTEnv()
    for i=1:1
        numQubits = rand(3:5)
        qureg = QuEST.createQureg(numQubits, env)
        state = rand(Complex{qreal}, 2^numQubits)
        state /= norm(state)
        QuEST.initStateFromAmps(qureg, real.(state), imag.(state))

        original_stdout = stdout
        read_pipe, write_pipe = redirect_stdout()
        lines = Vector{String}()
        QuEST.reportStateToScreen(qureg, env, 1)
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
            this_re = parse(qreal, split(lines[l+2], " ")[1][1:end-1])
            this_im = parse(qreal, split(lines[l+2], " ")[2])
            @test this_im ≈ imag(state[l]) atol = tolerance
            @test this_re ≈ real(state[l]) atol = tolerance
        end

        QuEST.reportState(qureg)

        io = open("state_rank_0.csv", "r")

        @test readline(io) == "real, imag"

        for l=1:2^numQubits
            this_line = readline(io)
            this_re = parse(qreal, split(this_line, " ")[1][1:end-1])
            this_im = parse(qreal, split(this_line, " ")[2])
            @test this_im ≈ imag(state[l]) atol = tolerance
            @test this_re ≈ real(state[l]) atol = tolerance
        end

        @test readline(io) == ""
        
        close(io)

        rm("state_rank_0.csv")

    end

    QuEST.destroyQuESTEnv(env)
    
end

function test_seed()
    env= QuEST.createQuESTEnv()

    num_seeds = rand(1:30)
    seeds = rand(Culong, num_seeds)
    QuEST.seedQuESTDefault()
    QuEST.seedQuEST(seeds)

    QuEST.destroyQuESTEnv(env)

end

function test_sync()
    env= QuEST.createQuESTEnv()
    numQubits = rand(3:12)
    qureg = QuEST.createQureg(numQubits, env)
    state = rand(Complex{qreal}, 2^numQubits)
    state /= norm(state)
    QuEST.initStateFromAmps(qureg, real.(state), imag.(state))

    QuEST.syncQuESTEnv(env)
    @test QuEST.syncQuESTSuccess(1) == 1

    QuEST.destroyQureg(qureg, env)
    QuEST.destroyQuESTEnv(env)

end
