# QuEST_jl/src/base/QASM_logging.jl
#

function clearRecordedQASM(qureg ::QuEST_Types.Qureg) ::Nothing
    ccall(:clearRecordedQASM, Cvoid, (QuEST_Types.Qureg,), qureg)
    return nothing
end

function printRecordedQASM(qureg ::QuEST_Types.Qureg) ::Nothing
    ccall(:printRecordedQASM, Cvoid, (QuEST_Types.Qureg,), qureg)
    return nothing
end

function startRecordingQASM(qureg ::QuEST_Types.Qureg) ::Nothing
    ccall(:startRecordingQASM, Cvoid, (QuEST_Types.Qureg,), qureg)
    return nothing
end

function stopRecordingQASM(qureg ::QuEST_Types.Qureg) ::Nothing
    ccall(:stopRecordingQASM, Cvoid, (QuEST_Types.Qureg,), qureg)
    return nothing
end

function writeRecordedQASMToFile(qureg    ::QuEST_Types.Qureg,
                                 filename ::String) ::Nothing
    @assert begin
        ios = open(filename, "w")
        close(ios) === nothing
    end

    ccall(:writeRecordedQASMToFile, Cvoid, (QuEST_Types.Qureg, Cstring), qureg, filename)
    return nothing
end
