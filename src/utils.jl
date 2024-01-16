"""
    Gets all function, structs and @enums by names and returns a string
"""
function get_function_struct_names(filename)
    # Read the file content into a string
    content = read(filename, String)
  
    # Define a regular expression pattern to capture function and struct names
    pattern = r"function\s+(\w+)*\(|struct\s+(\w+)|@enum\s+(\w+)"
  
  
    res = []
    # Extract all matches from the file content
    matches = eachmatch(pattern, content)
    for m in matches
      results = m.captures
      results = filter(x->!isnothing(x),results)[1]
      push!(res,results)
    end
    
    res
    return res
end


"""
    Takes an integer index and returns the corresponding index in C,
    namely index - 1
"""
function c_shift_index(index::Integer)
  index -= 1
end


"""
    unsafe_load_state_vec(state_vec_pointer, num_elements)

Load a state vector from memory.

This function takes a pointer to a state vector and the number of elements in the vector. 
It then loads the real and imaginary parts of the state vector from memory and returns a 
complex array representing the state vector.

# Arguments
- `state_vec_pointer`: A pointer to a state vector. This should be a struct with `real` and `imag` fields, 
  each of which is a pointer to a `Float64`.
- `num_elements`: The number of elements in the state vector.

# Returns
- A `Complex{Float64}` array representing the state vector.

# Warning
This function uses the `unsafe_load` function to directly load data from memory, which can be unsafe 
if the pointers do not point to valid memory. Use with caution.

"""
function unsafe_load_state_vec(state_vec_pointer,num_elements)
  reals = [unsafe_load(state_vec_pointer.real,q) for q in Base.OneTo(num_elements)]
  imags = [unsafe_load(state_vec_pointer.imag,q) for q in Base.OneTo(num_elements)]
  Complex.(reals,imags) 
end


function complex(c::QComplex)
  Complex(c.real, c.imag)
end

function qComplex(c::Complex)
  QComplex(c.re, c.im)
end

function vector(v::QVector)
  [v.x, v.y, v.z]
end

function qVector(v::Vector)
  QVector(v[1], v[2], v[3])
end


function make_QuEST_matrix_2x2(U ::Matrix{Complex{Float64}}) ::ComplexMatrix2
  size(U) != (2,2) && error("U must be 2x2")
  u = ComplexMatrix2(
      ( (real(U[1,1]), real(U[1,2])), (real(U[2,1]), real(U[2,2])) ),
      ( (imag(U[1,1]), imag(U[1,2])), (imag(U[2,1]), imag(U[2,2])) )  )
  return u
end

function make_QuEST_matrix_4x4(U ::Matrix{Complex{Float64}}) ::ComplexMatrix4
  size(U) != (4,4) && error("U must be 4x4")
  u = ComplexMatrix4( ( ( real(U[1,1]), real(U[1,2]), real(U[1,3]), real(U[1,4]) ),
                        ( real(U[2,1]), real(U[2,2]), real(U[2,3]), real(U[2,4]) ),
                        ( real(U[3,1]), real(U[3,2]), real(U[3,3]), real(U[3,4]) ),
                        ( real(U[4,1]), real(U[4,2]), real(U[4,3]), real(U[4,4]) ) ),
                      ( ( imag(U[1,1]), imag(U[1,2]), imag(U[1,3]), imag(U[1,4]) ),
                        ( imag(U[2,1]), imag(U[2,2]), imag(U[2,3]), imag(U[2,4]) ),
                        ( imag(U[3,1]), imag(U[3,2]), imag(U[3,3]), imag(U[3,4]) ),
                        ( imag(U[4,1]), imag(U[4,2]), imag(U[4,3]), imag(U[4,4]) ) )  )
end

function make_QuEST_matrix_NxN(mat ::Matrix{Complex{Float64}}) ::ComplexMatrixN
  (R,C) = size(mat)
  @assert R==C
  @assert R ≥ 2
  num_qubits = Int(round(log2(R)))
  @assert 2^num_qubits == R
  @assert 1 ≤ num_qubits ≤ 50
  mat_Q = createComplexMatrixN(num_qubits)
  reals = []
  imags = []
  for c = 1:R                       # Julia is column major
      for r = 1:R
          push!(reals,real(mat[r,c]))
          push!(imags,imag(mat[r,c]))
      end
  end

  initComplexMatrixN(mat_Q, pointer(reals), pointer(imags))
  mat_Q
end

function make_QuEST_matrix(M ::Matrix{Base.Complex{Float64}}) ::ComplexMatrixN
  (R,C) = size(M)
  @assert R==C
  @assert R ≥ 2
  numQubits = Int(round( log2(R) ))
  @assert 2^numQubits == R
  @assert 1 ≤ numQubits ≤ 50

  MQ = createComplexMatrixN(numQubits)
  for c = 1:R                       # Julia is column major
      for r = 1:R
          re_MQ_r = unsafe_load(MQ.real,r)
          im_MQ_r = unsafe_load(MQ.imag,r)
          unsafe_store!(re_MQ_r, real(M[r,c]), c)
          unsafe_store!(im_MQ_r, imag(M[r,c]), c)
      end
  end
  return MQ
end

pauliX_matrix()=Matrix([0 1.0;1 0])
pauliY_matrix()=Matrix([0 -im;im 0])
pauliZ_matrix()=Matrix([1.0 0;0 -1.0])
hadamard_matrix()=Matrix([1.0 1.0;1.0 -1.0])/sqrt(2)

function create_state_vec_zero_state(num_qubits)
  v = zeros(Float64, 2^num_qubits)
  v[1]=1
  v
end

function apply_hadamard_to_state_vector_at_target_index(state_vec, target)
  H = hadamard_matrix()
  num_qubits = log2(length(state_vec)) |> Int
  i_right=1.0*Matrix(I, 2^(target-1), 2^(target-1))
  un_changed_qubits = num_qubits-target
  i_left =1.0*Matrix(I, 2^un_changed_qubits, 2^un_changed_qubits)
  U=kron(i_left, H, i_right)
  U*state_vec
end

function get_pauli_code(str)
  str == "I" ? pauliOpType(0) :
  str == "X" ? pauliOpType(1) :
  str == "Y" ? pauliOpType(2) :
  str == "Z" ? pauliOpType(3) :
  error("Invalid Pauli code: only I, X, Y, Z are allowed")
end

function test_is_probability(prob)
  !(0.0 <= prob <= 1.0) && error("prob must be a probability, between 0 and 1")
end

function test_is_in_qubit_range(qureg::Qureg, targetQubit)
  qureg.numQubitsRepresented < targetQubit && error("Qubit $targetQubit is not in the qubit range of the qureg")
end

function test_is_valid_probability(prob,threshold)
  test_is_probability(prob)
  prob > threshold && error("probability, ($(prob)) must be less than $(threshold)")
end

function test_is_density_matrix(qureg::Qureg)
  qureg.isDensityMatrix == 0  && error("Qureg must be a density matrix")
end

function test_kraus_sum_is_identity(operators)
  num_ops = length(operators)
  dims = size(operators[1])
  identity = Matrix(I,dims...)
  M_ident = map(x->operators[x]'*operators[x], 1:num_ops) |> sum
  !isapprox(M_ident,identity,atol=1e-15) && error("Kraus sum is not identity")
end

function test_max_kraus_operators(operators,max_number)
  num_ops = length(operators)
  num_ops > max_number && error("Number of operators, $(num_ops), must be less than $(max_number)")
end