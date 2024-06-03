# QuEST.jl/src/QuEST_jl.jl
#
# Authors:
#  - Dirk Oliver Theis, Ketita Labs & Uni Tartu
#  - Bahman Ghandchi, Uni Tartu
#
# MIT License
#
# (c) Ketita Labs, Uni Tartu, and the authors.
#
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
module QuEST_jl

export QuEST32, QuEST64


module _QuEST_jl_Build_Setup

    include(joinpath("..","deps","build_setup.jl"))

end #^ module _QuEST_jl_Build_Setup

using Libdl: dlopen, dlclose, RTLD_LAZY, RTLD_DEEPBIND, RTLD_GLOBAL

const _QUEST_LIB = Ref{Ptr{Cvoid}}(0)

"""
Function `_dlopen_QuEST()` — loads QuEST dll, returns precision code.
"""
function _dlopen_QuEST(num::Int) ::Cint
    @assert num ∈ [32,64]
    let numstr = string(num)

        if _QuEST_jl_Build_Setup.EXPERT_BUILD
            global _QUEST_LIB[] =
                dlopen("libQuEST_"*numstr,
                       RTLD_LAZY|RTLD_DEEPBIND|RTLD_GLOBAL)
        else
            global _QUEST_LIB[] =
                dlopen(joinpath(@__DIR__,"..","deps","QuEST","build"*numstr,"QuEST","libQuEST"),
                       RTLD_LAZY|RTLD_DEEPBIND|RTLD_GLOBAL)
        end

    end #^ let

    return ccall(:getQuEST_PREC , Cint, (),)
end #^ fn _dlopen_QuEST

function _dlclose_QuEST() ::Nothing
    dlclose(_QUEST_LIB[])
    _QUEST_LIB[] = 0
    nothing;
end

"""
Functions `_quest_mtx_𝑛()` for 𝑛 ∈ {2,4} — convert Julia matrix into QuEST matrix
"""
function _quest_mtx_2(U ::Matrix{Complex{FLT}}) ::QuEST_Types.ComplexMatrix2 where FLT <: Real
    @assert size(U) == (2,2)
    u = QuEST_Types.ComplexMatrix2(
        ( (real(U[1,1]), real(U[1,2])), (real(U[2,1]), real(U[2,2])) ),
        ( (imag(U[1,1]), imag(U[1,2])), (imag(U[2,1]), imag(U[2,2])) )  )
    return u
end

"""
Functions `_quest_mtx_𝑛()` for 𝑛 ∈ {2,4} — convert Julia matrix into QuEST matrix
"""
function _quest_mtx_4(U ::Matrix{Complex{FLT}}) ::QuEST_Types.ComplexMatrix4 where FLT <: Real
    @assert size(U) = (4,4)
    u = QuEST_Types.ComplexMatrix4(
        ( ( real(U[1,1]), real(U[1,2]), real(U[1,3]), real(U[1,4]) ),
          ( real(U[2,1]), real(U[2,2]), real(U[2,3]), real(U[2,4]) ),
          ( real(U[3,1]), real(U[3,2]), real(U[3,3]), real(U[3,4]) ),
          ( real(U[4,1]), real(U[4,2]), real(U[4,3]), real(U[4,4]) ) ),
        ( ( imag(U[1,1]), imag(U[1,2]), imag(U[1,3]), imag(U[1,4]) ),
          ( imag(U[2,1]), imag(U[2,2]), imag(U[2,3]), imag(U[2,4]) ),
          ( imag(U[3,1]), imag(U[3,2]), imag(U[3,3]), imag(U[3,4]) ),
          ( imag(U[4,1]), imag(U[4,2]), imag(U[4,3]), imag(U[4,4]) ) )  )
end


################################################################################
# 32 bit QuEST
################################################################################

"""
Module `QuEST`𝑥𝑦 — Julia wrapper for QuEST with 𝑥𝑦-bit floating point precision.
"""
module QuEST32
export QuEST_Types
export QubitIdx

import .._dlopen_QuEST

function QuEST_init() ::Nothing
    prec = _dlopen_QuEST(32)
    @assert prec == 1 "Wrong precision. Please rebuild the package"
end

module QuEST_Types  # 32 bit qreal
    using CEnum

    const MPI_FLOAT  = Cfloat   # This is needed, because ...
    const MPI_DOUBLE = Cdouble  # ...

    include(joinpath("..","deps","questclang_common_32.jl"))

end #^ module QuEST_Types

const QubitIdx = Cint
const Qreal    = Float32
@assert sizeof(Qreal)==sizeof(QuEST_Types.qreal)

#
# Include C-wrappers, which will use the correct data types
#

using .QuEST_Types

include("data_structure_functions.jl")
include("QASM_logging.jl")
include("debugging.jl")
include("operators.jl")
include("decoherence.jl")
include("state_init.jl")
include("unitaries.jl")
include("calculations.jl")
include("gates.jl")

end #^ module QuEST32

################################################################################
# 64 bit QuEST
################################################################################

"""
Module `QuEST`𝑥𝑦 — Julia wrapper for QuEST with 𝑥𝑦-bit floating point precision.
"""
module QuEST64
export QuEST_Types
export QubitIdx

import .._dlopen_QuEST

function QuEST_init() ::Nothing
    prec = _dlopen_QuEST(64)
    @assert prec == 2 "Wrong precision. Please rebuild the package"
end

module QuEST_Types  # 64 bit qreal
    using CEnum

    const MPI_FLOAT  = Cfloat   # This is needed, because ...
    const MPI_DOUBLE = Cdouble  # ...

    include(joinpath("..","deps","questclang_common_64.jl"))

end #^ module QuEST_Types

const QubitIdx = Cint
const Qreal    = Float64
@assert sizeof(Qreal)==sizeof(QuEST_Types.qreal)

#
# Include C-wrappers, which will use the correct data types
#

using .QuEST_Types

include("data_structure_functions.jl")
include("QASM_logging.jl")
include("debugging.jl")
include("operators.jl")
include("decoherence.jl")
include("state_init.jl")
include("unitaries.jl")
include("calculations.jl")
include("gates.jl")

end #^ module QuEST64

end # module QuEST_jl
# EOF
