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
