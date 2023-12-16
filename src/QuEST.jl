module QuEST

    using Reexport
    @reexport using QuEST_jll
    @reexport using Clang
    @reexport using JuliaFormatter
    include("utils.jl")
   
    wrappers_file = joinpath(@__DIR__,"C/quest_julia_c_wrapper.jl")
    exports = Symbol.(get_function_struct_names(wrappers_file))
    
export 
    exports

    include("C/quest_julia_c_wrapper.jl")
end
