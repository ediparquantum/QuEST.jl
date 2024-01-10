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
