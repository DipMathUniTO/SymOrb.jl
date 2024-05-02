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
    global Ω, K, iterations, action_grad, action_value, trajectory , group_order, fourier_coeff


    LGS_from_config(config)

    Γ = (project ∘ random_starting_path)()
    v = flatten(Γ)

    show_infos = false

    res = optimize(action, ∇action, v, BFGS(), Optim.Options(g_tol = 1e-8, show_trace=false, extended_trace=true,
            callback=x -> begin
                if (show_infos && (x.iteration % 10 == 0))
                    println(x.iteration, "\t", x.value, "\t", x.g_norm)
                end
                return false
            end); inplace=false)
     
    fourier_coeff = coeff_from_minimizer(res.minimizer)
    trajectory = path_from_minimizer(res.minimizer)
    iterations = res.iterations
    action_grad = res.g_residual
    action_value =  res.minimum
    return res.iteration_converged ||  res.x_converged || res.f_converged || res.g_converged
end



end # module symorb
