module SymorbJL

using GAP
using OffsetArrays
using Optim
using LinearAlgebra
using Plots

include("globals.jl")
include("path.jl")
include("initdata.jl")
include("matrices.jl")
include("action.jl")
include("projectors.jl")


export symorb_minimize

check_convergence(res)::Bool = res.iteration_converged ||  res.x_converged || res.f_converged || res.g_converged

function symorb_minimize(config::AbstractDict; show_infos = false)

    LGS_from_config(config)

    starting_path = random_starting_path() 
    
    Γ = project(starting_path)
    v = flatten(Γ)

    method = :BFGS
    method_callable = getfield(Optim, method)


    options = Optim.Options(
        g_tol = 1e-8,
        iterations = Int(1e5),
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
     
    minimizer = (project ∘ emboss)(res.minimizer)
    return MinimizationResult( 
        fourier_coeff = minimizer,
        trajectory = build_path(minimizer),
        iterations =  res.iterations,
        action_value =  res.minimum,
        gradient_norm = res.g_residual,
        converged = check_convergence(res),
        initial = starting_path,
        minimization_method = method,  
        minimization_options = MinimizationOptions(options)    
    );

end


function plot_path(path)
    pl = plot(aspect_ratio=:equal)
    for i ∈ 1:N
        plot!(pl, [[path[j][i][h] for j ∈ 0:steps+1] for h ∈ 1:dim]..., label="Body $i", linewidth=3,  aspect_ratio=:equal)
    end
    display(plot(pl))
    return false
end



end # module symorb
