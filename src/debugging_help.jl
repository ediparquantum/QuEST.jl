struct ErrorQubitIndexOutofBoundsFromQureg end
function throw_error(::ErrorQubitIndexOutofBoundsFromQureg,qubit_index,qureg)
    num_qubits = qureg.numQubitsRepresented
    error("Qubit index ($(qubit_index)) is out of bounds of quantum state ($(num_qubits))")
end

function check_qubits_in_range(qureg,qubits)
    max_qubit = qureg.numQubitsRepresented
    if !all(qubits .< max_qubit && qubits .> 0)
        throw_error(ErrorQubitIndexOutofBoundsFromQureg(),qubit_index,qureg)
    end
end