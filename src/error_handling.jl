abstract type AbstractErrorMessage end
abstract type AbstractJulia2CSyntexError <: AbstractErrorMessage end
abstract type AbstractDecoherenceError <: AbstractErrorMessage end




struct QubitsNotInQuregError <: AbstractJulia2CSyntexError end
function throw_error(::QubitsNotInQuregError)
    error("Qubits not in qureg, Segmentation fault will occur. Throwing error to prevent this.")
end
function test_qubit_present(qureg::Qureg,qubits)
    
    if !(all(qubits .∈ Ref(Base.OneTo(qureg.numQubitsRepresented))))
        throw_error(QubitsNotInQuregError())
    end
end


struct QubitsNotInQuregMultiControllError <: AbstractJulia2CSyntexError end
function throw_error(::QubitsNotInQuregMultiControllError)
    error("Qubits not in qureg - 2, Segmentation fault will occur. Throwing error to prevent this.")
end
function test_qubit_present_multi_control(qureg,controlQubits)
    if !(0 ≤ length(controlQubits) ≤ qureg.numQubitsRepresented -2)
        throw_error(QubitsNotInQuregMultiControllError())
    end
end


struct RowColNotInQuregError <: AbstractJulia2CSyntexError end
function throw_error(::RowColNotInQuregError)
    error("Qubits not in qureg, Segmentation fault will occur. Throwing error to prevent this.")
end
function test_row_col_in_size(qureg::Qureg,row_col::Int)
    num_qubits = qureg.numQubitsRepresented
    num_row_or_col = 2*2^num_qubits
    if !(row_col ≤ num_row_or_col)
        throw_error(RowColNotInQuregError())
    end
end




struct TwoQubitsMustBeDifferent <: AbstractJulia2CSyntexError end
function throw_error(::TwoQubitsMustBeDifferent)
    error("Situation: qubits must be different, but they are the same")
end
function test_qubits_different(qubit_1::Int,qubit_2::Int)
    if qubit_1 == qubit_2
        throw_error(TwoQubitsMustBeDifferent())
    end
end




struct PauliCodeLengthErrorPauliProd <: AbstractJulia2CSyntexError end
function throw_error(::PauliCodeLengthErrorPauliProd)
    error("targetQubits and pauliCodes must have the same length")
end
function test_length_pauli_prod(targetQubits,pauliCodes)
    if !(length(targetQubits) == length(pauliCodes))
        throw_error(PauliCodeLengthErrorPauliProd())
    end
end




struct PauliCodeLengthErrorPauliSum <: AbstractJulia2CSyntexError end
function throw_error(::PauliCodeLengthErrorPauliSum)
    error("There must be the same number of Pauli codes (I,X,Y or Z) as length total sum terms multiplied by the number of qubits represented by qureg (e.g. 2 qubits and 2 terms = 4 codes)")
end
function test_length_pauli_sum(qureg,numSumTerms,allPauliCodes)
    if !(numSumTerms*qureg.numQubitsRepresented == length(allPauliCodes))
        throw_error(PauliCodeLengthErrorPauliProd())
    end
end




struct KrausNotSquareError <: AbstractDecoherenceError end
function throw_error(::KrausNotSquareError)
    error("Kraus operators must be square") 
end
function test_kraus_operator_dimension_square(rows,cols)
    if !(rows == cols)
        throw_error(KrausNotSquareError())
    end
end



struct VectorHasRepetitions <: AbstractJulia2CSyntexError end
function throw_error(::VectorHasRepetitions)
    error("Vector has repeitions") 
end
function test_if_vector_has_repetitions(v)
    if length(v) != length(Set(v))
        throw_error(VectorHasRepetitions())
    end
end


struct VectorsNotMutuallyExclusive <: AbstractJulia2CSyntexError end
function throw_error(::VectorsNotMutuallyExclusive)
    error("Vectors are not mutually exclusive") 
end
function test_if_vec1_is_in_vec2_vise_versa(vec1,vec2)
    v1v2 = [i .∈ vec1 for i in vec2][1] |> any
    v2v1 = [i .∈ vec2 for i in vec1][1] |> any
    if v1v2 || v2v1
        throw_error(VectorsNotMutuallyExclusive())
    end
end


struct LengthVectorNotSufficientMustBeAtLeastOne  <: AbstractJulia2CSyntexError end
function throw_error(::LengthVectorNotSufficientMustBeAtLeastOne)
    error("Length must be at least one") 
end
function test_if_vec_has_length_at_least_one(vec)
  
    if length(vec) < 1
        throw_error(LengthVectorNotSufficientMustBeAtLeastOne())
    end
end