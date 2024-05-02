module SymorbJL

using GAP
using OffsetArrays
using Optim
using LinearAlgebra

include("globals.jl")
include("path.jl")
include("initdata.jl")
include("matrices.jl")
include("action.jl")
include("projectors.jl")


export symorb_minimize


function symorb_minimize(config::AbstractDict)
    global Ω, K

    println("--- SymorbJL ---")
    println(">> Generating LGS...")
    LGS_from_config(config)

    println(">> Minimizing...")

    Γ = (project ∘ random_starting_path)()
    v = flatten(Γ)

    show_infos = true

    res = optimize(action, ∇action, v, BFGS(), Optim.Options(g_tol = 1e-8, show_trace=false, extended_trace=true,
            callback=x -> begin
                if (show_infos && (x.iteration % 10 == 0))
                    println(x.iteration, "\t", x.value, "\t", x.g_norm)
                end
                return false
            end); inplace=false)
    return (build_path ∘ project ∘ emboss)(res.minimizer)
end



end # module symorb
