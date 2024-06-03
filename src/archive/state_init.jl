# QuEST_jl/src/base/state_init.jl
#

function cloneQureg(targetQureg ::QuEST_Types.Qureg, copyQureg ::QuEST_Types.Qureg) ::Nothing

    return ccall(:cloneQureg, Cvoid, (QuEST_Types.Qureg, QuEST_Types.Qureg), targetQureg, copyQureg)

end

function initBlankState(qureg ::QuEST_Types.Qureg) ::Nothing

    ccall(:initBlankState, Cvoid,
          (QuEST_Types.Qureg,),
          qureg)
end

function initClassicalState(qureg    ::QuEST_Types.Qureg,
                            stateIdx ::Integer)           ::Nothing

    ccall(:initClassicalState, Cvoid,
          (QuEST_Types.Qureg, Clonglong),
          qureg,              stateIdx)
end

function initPlusState(qureg ::QuEST_Types.Qureg) ::Nothing

    ccall(:initPlusState, Cvoid,
          (QuEST_Types.Qureg,),
          qureg)
end

function initPureState(qureg ::QuEST_Types.Qureg,
                       pure  ::QuEST_Types.Qureg) ::Nothing

    @assert qureg.numQubitsRepresented == pure.numQubitsRepresented
    @assert pure.isDensityMatrix == 0

    ccall(:initPureState, Cvoid,
          (QuEST_Types.Qureg,  QuEST_Types.Qureg),
          qureg,               pure)
end

function initStateFromAmps(qureg     ::QuEST_Types.Qureg,
                           amps_real ::Vector{Qreal},
                           amps_imag ::Vector{Qreal}      ) ::Nothing

    @assert length(amps_real) == qureg.numAmpsTotal
    @assert length(amps_imag) == qureg.numAmpsTotal

    ccall(:initStateFromAmps,
          Cvoid,
          (QuEST_Types.Qureg, Ptr{Qreal}, Ptr{Qreal}),
          qureg,              amps_real,  amps_imag)
end

function initZeroState(qureg ::QuEST_Types.Qureg) ::Nothing

    ccall(:initZeroState, Cvoid,
          (QuEST_Types.Qureg,),
          qureg)
end

function setAmps(qureg        ::QuEST_Types.Qureg,
                 startIdx     ::Integer,
                 amps_real    ::Vector{Qreal},
                 amps_imag    ::Vector{Qreal},
                 numAmps      ::Clonglong)   :: Nothing

    @assert numAmps             == length(amps_real)
    @assert numAmps             == length(amps_imag)
    @assert startIdx + numAmps  <= qureg.numAmpsTotal

    ccall(:setAmps, Cvoid,
          (QuEST_Types.Qureg, Clonglong, Ptr{Qreal}, Ptr{Qreal}, Clonglong),
          qureg,              startIdx,  amps_real,  amps_imag,  numAmps)
end



const __ℂ     = QuEST_Types.Complex
const __Qureg = QuEST_Types.Qureg

function setWeightedQureg(factor1    ::Complex{Qreal},
                          qureg1     ::QuEST_Types.Qureg,
                          factor2    ::Complex{Qreal},
                          qureg2     ::QuEST_Types.Qureg,
                          factorOut  ::Complex{Qreal},
                          quregOut   ::QuEST_Types.Qureg)  ::Nothing

    @assert Bool(qureg1.isDensityMatrix) == Bool(qureg2.isDensityMatrix) == Bool(quregOut.isDensityMatrix)
    @assert qureg1.numQubitsRepresented  == qureg2.numQubitsRepresented  == quregOut.numQubitsRepresented

    local factor1_quest   = __ℂ(real(factor1  ),imag(factor1  ))
    local factor2_quest   = __ℂ(real(factor2  ),imag(factor2  ))
    local factorOut_quest = __ℂ(real(factorOut),imag(factorOut))

    ccall(:setWeightedQureg, Cvoid,
          (__ℂ,          __Qureg,  __ℂ,           __Qureg, __ℂ,             __Qureg),
          factor1_quest, qureg1,   factor2_quest, qureg2,  factorOut_quest, quregOut)
end
