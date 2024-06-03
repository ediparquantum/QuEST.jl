##################################################################
# Filename  : multiRotatePauli_fix_errors_remove_after_branch_merged_tested.jl
# Author    : Jonathan Miller
# Date      : 2024-05-14
# Aim       : aim_script
#           :
#           :
#           : Dummy script, just used in dev to test minor changes
#           : Last minor change was on 2024-05-14 working on 
#           : multiRotatePauli function
#           :
#           :
##################################################################



using Pkg
Pkg.activate(".")
include("src/QuEST.jl")
using .QuEST

num_qubits = 5
env = createQuESTEnv()
  
# create a 2 qubit register in the zero state
qureg = createQureg(num_qubits, env)
initZeroState(qureg)

multiRotatePauli(qureg, [1,2,3,4], [PAULI_I,PAULI_X,PAULI_Y,PAULI_Z], .1)  #<-- Use enum directly6
multiRotatePauli(qureg, [1,2,3,4], [1,1,2,3], .1)  #<-- Or use integer values

# unload QuEST
destroyQureg(qureg, env)
destroyQuESTEnv(env)

