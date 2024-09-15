module SymPath

using LinearAlgebra
using ForwardDiff: derivative
using GAP: GapObj, gap_to_julia as g2j, @gap, evalstr, Globals as GG, Packages
using NLsolve: nlsolve, NLsolve
using Optim: optimize, Options, Optim
using Printf: @sprintf
using Base: @kwdef, kron
using RuntimeGeneratedFunctions
import TOML

RuntimeGeneratedFunctions.init(@__MODULE__)

include("types.jl")
include("path.jl")
include("initialize.jl")
include("matrices.jl")
include("action.jl")
include("projectors.jl")
include("minimization.jl")

export find_orbits, initialize_and_find_orbits, initialize
export print_path_to_file, path_animation, refine_path, read_path_from_file
export Newton, BFGS, ConjugateGradient, Methods, CompoundMethod, OneMethod, TwoMethods, MinimizationResult
export morse_index

end
