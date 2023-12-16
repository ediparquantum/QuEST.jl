module QuEST

    using Reexport
    @reexport using QuEST_jll
    @reexport using Clang
    @reexport using JuliaFormatter
    include("utils.jl")
   
    wrappers_file = "src/C/quest_julia_c_wrapper.jl"
    exports = Symbol.(get_function_struct_names(wrappers_file))
export 
    exports

    include("gen/wrapper.jl")
    include("C/quest_julia_c_wrapper.jl")
end
