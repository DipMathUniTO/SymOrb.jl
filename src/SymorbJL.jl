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

check_convergence(res)::Bool = res.iteration_converged ||  res.x_converged || res.f_converged || res.g_converged

function symorb_minimize(config::AbstractDict)

    LGS_from_config(config)

    starting_path = random_starting_path() 
    
    Γ = project(starting_path)
    v = flatten(Γ)

    method = :BFGS
    method_callable = getfield(Optim, method)

    show_infos = false

    options = Optim.Options(
        g_tol = 1e-8,
        show_trace = false,
        extended_trace = true,
        callback = x -> begin
            if (show_infos && (x.iteration % 10 == 0))
                println(x.iteration, "\t", x.value, "\t", x.g_norm)
            end
            return false
        end
    )

    res = Optim.optimize(action, ∇action, v, method_callable(), options; inplace=false)
     

    return MinimizationResult( 
        fourier_coeff = coeff_from_minimizer(res.minimizer),
        trajectory = path_from_minimizer(res.minimizer),
        iterations =  res.iterations,
        action_value =  res.minimum,
        gradient_norm = res.g_residual,
        converged = check_convergence(res),
        initial = starting_path,
        minimization_method = method,  
        minimization_options = MinimizationOptions(options)    
    )

end




end # module symorb
