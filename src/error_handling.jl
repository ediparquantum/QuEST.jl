abstract type ErrorMessage end
abstract type Julia2CSyntexError <: ErrorMessage end
struct QubitsNotInQuregError <: Julia2CSyntexError end


function throw_error(::QubitsNotInQuregError)
    error("Qubits not in qureg, Segmentation fault will occur. Throwing error to prevent this.")
end


function test_qubit_present(qureg::Qureg,qubits::Union{Int,Array{Int}})
    
    if !(all(qubits .∈ Ref(Base.OneTo(qureg.numQubitsRepresented))))
        throw_error(QubitsNotInQuregError())
    end
end






