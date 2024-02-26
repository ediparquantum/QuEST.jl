##################################################################
# Filename  : partial_trace.jl
# Author    : Jonathan Miller
# Date      : 2024-02-26
# Aim       : aim_script
#           : Using the package TensorOperations
#           : to calculate the partial trace of a density matrix
##################################################################


"""
    per_qubit_partial_trace(density_matrix::AbstractArray,qubits_to_keep::Vector{Int})

Calculate the partial trace of a density matrix for a given set of qubits to keep, can choose qubits by index, startin with 1.
"""
function per_qubit_partial_trace(density_matrix::AbstractArray,qubits_to_keep::Vector{Int})
    system_total_size = size(density_matrix,1) # Get n for n x n matrix
    num_qubits = Int(log2(system_total_size)) # Number of qubits is log2 of n
    substem_sizes_per_qubit = fill(2,num_qubits) # Each qubit is a 2x2 matrix so 2 for each qubit
    qubit_indices = Base.OneTo(num_qubits) # Create a range of qubit indices
    traceidx = [qubit_indices; qubit_indices] # Create a trace index for the tensor trace
    traceidx[qubits_to_keep] .+= num_qubits # Add n to the qubit indices to keep
    tensor = reshape(density_matrix, [substem_sizes_per_qubit; substem_sizes_per_qubit]...) # Reshape the density matrix into a tensor
    keepdim = prod([size(tensor, x) for x in qubits_to_keep]) # Get the dimension of the traced out density matrix
    reshape(tensortrace(tensor, Tuple(traceidx)), keepdim, keepdim) # Perform the tensor trace and reshape the result
end

function get_partial_trace(qureg::Qureg,qubits_to_keep::Vector{Int})
    density_matrix = get_qureg_matrix(qureg)
    per_qubit_partial_trace(density_matrix,qubits_to_keep)
end

function get_per_qubit_trace(qureg::Qureg)
    num_qubits = get_num_qubits(qureg)
    density_matrix = get_qureg_matrix(qureg)
    [per_qubit_partial_trace(density_matrix,Int64[i]) for i in Base.OneTo(num_qubits)]
end


