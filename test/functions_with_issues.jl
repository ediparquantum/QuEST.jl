function check_ComplexMatrixN()
    num_qubits = 3
    M = createComplexMatrixN(num_qubits)
    @test M.numQubits == num_qubits
    M_j_real = rand(Float64, 2^num_qubits * 2^num_qubits)
    M_j_imag = rand(Float64, 2^num_qubits * 2^num_qubits)
    initComplexMatrixN(M, M_j_real, M_j_imag)
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