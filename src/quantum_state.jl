##################################################################
# Filename  : quantum_state.jl
# Author    : Jonathan Miller
# Date      : 2024-02-26
# Aim       : aim_script
#           : Functions to get state vector and density matrix from Qureg
#           : based on Qureg
##################################################################


struct DensityMatrix end
struct StateVector end

function state_type(qureg::Qureg)
    if Bool(qureg.isDensityMatrix)
        return DensityMatrix()
    else
        return StateVector()
    end
end


function get_qureg_matrix(::StateVector,qureg::Qureg)
    reals = unsafe_wrap(Vector{Float64},qureg.stateVec.real,qureg.numAmpsTotal) |> x -> chop.(x)
    imags = unsafe_wrap(Vector{Float64},qureg.stateVec.imag,qureg.numAmpsTotal) |> x -> chop.(x)
    Complex.(reals,imags)
end

function get_qureg_matrix(::DensityMatrix,qureg::Qureg)
    reals = unsafe_wrap(Vector{Float64},qureg.stateVec.real,qureg.numAmpsTotal) |> x -> chop.(x)
    imags = unsafe_wrap(Vector{Float64},qureg.stateVec.imag,qureg.numAmpsTotal) |> x -> chop.(x)
    reshape(Complex.(reals,imags),Int(sqrt(qureg.numAmpsTotal)),Int(sqrt(qureg.numAmpsTotal)))
end

function get_qureg_matrix(qureg::Qureg)
    state_t = state_type(qureg)
    get_qureg_matrix(state_t,qureg)
end