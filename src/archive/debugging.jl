# QuEST_jl/src/base/debugging.jl
#

function copyStateFromGPU(qureg ::QuEST_Types.Qureg) ::Nothing
    ccall(:copyStateFromGPU, Cvoid, (QuEST_Types.Qureg,), qureg)
    return nothing
end

function copyStateToGPU(qureg ::QuEST_Types.Qureg) ::Nothing
    ccall(:copyStateToGPU, Cvoid, (QuEST_Types.Qureg,), qureg)
    return nothing
end

function getEnvironmentString(env ::QuEST_Types.QuESTEnv, qureg ::QuEST_Types.Qureg) ::String
    cstr = Vector{Cchar}(undef,232)
    ccall(:getEnvironmentString, Cvoid, (QuEST_Types.QuESTEnv, QuEST_Types.Qureg, Ptr{Cchar}), env, qureg, cstr)

    return unsafe_string(pointer(cstr))
end

function initDebugState(qureg ::QuEST_Types.Qureg) ::Nothing
    ccall(:initDebugState, Cvoid, (QuEST_Types.Qureg,), qureg)
    return nothing
end

#
#void 	invalidQuESTInputError (const char *errMsg, const char *errFunc)
#An internal function called when invalid arguments are passed to a QuEST API call, which the user can optionally override by redefining.  More...
#

function reportPauliHamil(hamil     ::QuEST_Types.PauliHamil)   ::Nothing
    ccall(:reportPauliHamil, Cvoid, (QuEST_Types.PauliHamil,), hamil)
    return nothing
end

function reportQuESTEnv(env ::QuEST_Types.QuESTEnv) ::Nothing
    ccall(:reportQuESTEnv, Cvoid, (QuEST_Types.QuESTEnv,), env)
    return nothing
end

function reportQuregParams(qureg ::QuEST_Types.Qureg) ::Nothing
    ccall(:reportQuregParams, Cvoid, (QuEST_Types.Qureg,), qureg)
    return nothing
end

function reportState(qureg ::QuEST_Types.Qureg) ::Nothing
    ccall(:reportState, Cvoid, (QuEST_Types.Qureg,), qureg)
    return nothing
end

function reportStateToScreen(qureg      ::QuEST_Types.Qureg,
                             env        ::QuEST_Types.QuESTEnv,
                             reportRank ::Integer)              ::Nothing
    ccall(:reportStateToScreen,
          Cvoid,
          (QuEST_Types.Qureg, QuEST_Types.QuESTEnv, Cint),
          qureg,              env,                  reportRank)
    nothing
end

function seedQuEST(seedarray ::Base.Vector{Culong}) ::Nothing

    @assert  ! isempty(seedarray)

    ccall(:seedQuEST,
          Cvoid, (Ptr{Culong}, Cint),
          seedarray,
          Cint(length(seedarray)))
    return nothing
end

function seedQuESTDefault() ::Nothing
    ccall(:seedQuESTDefault, Cvoid, () )
    return nothing
end

function syncQuESTEnv(env ::QuEST_Types.QuESTEnv) ::Nothing
    ccall(:syncQuESTEnv, Cvoid,
          (QuEST_Types.QuESTEnv,),
          env)
    return nothing
end

function syncQuESTSuccess(successCode       ::Integer)       ::Cint
    return ccall(:syncQuESTSuccess, Cint,
                 (Cint,),
                 successCode)
end

#EOF
