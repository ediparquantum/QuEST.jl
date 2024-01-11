# QuEST_jl/src/base/gates.jl
#

function  collapseToOutcome(qureg        ::QuEST_Types.Qureg,
                            measureQubit ::Integer,
                            outcome      ::Integer)   ::Qreal

    return ccall(:collapseToOutcome, Qreal,
                 (QuEST_Types.Qureg, Cint,         Cint),
                 qureg,              measureQubit, outcome)
end

function measure(qureg        ::QuEST_Types.Qureg,
                 measureQubit ::Integer)            ::Int
    return ccall(:measure, Cint,
                 (QuEST_Types.Qureg, Cint),
                 qureg,              measureQubit)
end

function measureWithStats(qureg             ::QuEST_Types.Qureg,
                          measureQubit      ::Integer)            :: Tuple{Int,Qreal}

    outcomeProb = Ref{Qreal}(-1)

    outcome = ccall(:measureWithStats, Cint,
                    (QuEST_Types.Qureg, Cint,            Ref{Qreal}),
                    qureg,              measureQubit,    outcomeProb)

    return outcome, outcomeProb[]

end

#EOF
