module MinPath

    using LinearAlgebra
    using DelimitedFiles: readdlm, writedlm
    using ForwardDiff: derivative
    using GAP: GapObj, gap_to_julia as g2j, @gap, evalstr, Globals as GG, Packages
    using JSON: parsefile
    using NLsolve: nlsolve, NLsolve
    using OffsetArrays: OffsetVector, OffsetMatrix, OffsetArray, Origin
    using Optim: optimize, Options, Optim
    using Printf: @sprintf
    using Base: @kwdef

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

end
