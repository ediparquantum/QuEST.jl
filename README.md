# QuEST [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://fieldofnodes.github.io/QuEST.jl/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://fieldofnodes.github.io/QuEST.jl/dev/) [![Build Status](https://github.com/fieldofnodes/QuEST.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/fieldofnodes/QuEST.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Build Status](https://travis-ci.com/fieldofnodes/QuEST.jl.svg?branch=main)](https://travis-ci.com/fieldofnodes/QuEST.jl) [![Coverage](https://codecov.io/gh/fieldofnodes/QuEST.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/fieldofnodes/QuEST.jl)

## Details

1. Quantum Exact Simulation Toolkit is a high performance simulator of quantum circuits, state-vectors and density matrices.
2. `QuEST.jl` is a wrapper for [`QuEST`](https://github.com/QuEST-Kit/QuEST), which is written in `C`. `QuEST.jl` was compiled using [`BinaryBuilder.jl`](https://github.com/JuliaPackaging/BinaryBuilder.jl/tree/master), which is a system for the compilation of binary dependencies --whose aim it to *just work* anywhere the official Julia distribution does.
3. Function calls, enumerators, structs, etc were wrapped automatically using [`Clang.jl`](https://github.com/JuliaInterop/Clang.jl), which provided the needed `C` bindings for the Julia interface.
4. **Note:** Julia is indexed starting at $1$ and goes until $N$ items in an $N-$ element list. `C` indexes an $N-$ element list from $0$ to $N-1$.
5. To give the Julia user uniform interface, every `QuEST` call using an integer index, will first shift the index to be in the $0$ indexed form inside the Julia function. The user should be aware of this, but it is the aim of the `QuEST.jl` parckage wrapper writer(s) to handle this for the user.

## Getting started

Whilst we are working on getting this package registered with the official Julia registry, to use `QuEST.jl` 

```julia
using Pkg
Pkg.add("QuEST")
```

or

```julia
] add QuEST
```

```julia
using QuEST
num_qubits = 2
env = createQuESTEnv()
  
# create a 2 qubit register in the zero state
qubits = createQureg(num_qubits, env)
initZeroState(qubits)

# apply circuit
hadamard(qubits, 1)
controlledNot(qubits, 1, 2)
measure(qubits, 2)

# unload QuEST
destroyQureg(qubits, env)
destroyQuESTEnv(env)
```

## Information direct from the QuEST developers
For more in depth tutorials see (`C` based)

1. [QuEST's github](https://github.com/QuEST-Kit/QuEST)
2. [QuEST C based tutorial](https://github.com/QuEST-Kit/QuEST/blob/master/examples/README.md) - there are some code examples to get an idea of the flow
3. [QuEST website](https://quest.qtechtheory.org)
4. [QuEST's webiste documentation](https://quest.qtechtheory.org/docs/)
5. [Decoherence (noise) models](https://quest.qtechtheory.org/docs/decoherence/)
6. [QuEST's code documentation](https://quest-kit.github.io/QuEST/modules.html)
7. [QuESTlink - mathematica](https://github.com/QTechTheory/QuESTlink)
8. [QuEST whitepaper](https://www.nature.com/articles/s41598-019-47174-9)





## Acknowledgements

Some tests and functions were inspired and directly used from the public repository [Dighvijay: QuEST.jl](https://github.com/Dighvijay/QuEST.jl) made possible from the efforts of [Dighvijay](https://github.com/Dighvijay) on GitHub. We thank them for their valuable contributions.

