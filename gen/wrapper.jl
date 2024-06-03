##################################################################
# Filename  : wrapper.jl
# Author    : Jonathan Miller
# Date      : 2024-02-26
# Aim       : aim_script
#           : Script is really only used to generate functions for the C wrapper
#           : once that has been done, this script can be deleted (or saved and not used)
##################################################################


using QuEST_jll 
using Clang
using Clang.Generators
using JuliaFormatter



cd(@__DIR__)
julia_C_wrapper = "quest_julia_c_wrapper"
artifact_dir = QuEST_jll.artifact_dir
include_dir = normpath(artifact_dir, "include") 
headers = include_dir |> readdir
headers_dir = joinpath.(include_dir,headers)
options = load_options(joinpath(@__DIR__, "QuEST.toml"))
options["general"]["output_file_path"] = joinpath("..", "..","src", "C", "$(julia_C_wrapper).jl")
args = get_default_args()
push!(args, "-I$include_dir")
ctx = create_context(headers_dir, args, options)
build!(ctx)

