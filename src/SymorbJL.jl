module SymorbJL

using GAP, OffsetArrays, Optim, LinearAlgebra, Plots, NLsolve

include("globals.jl")
include("path.jl")
include("initdata.jl")
include("matrices.jl")
include("action.jl")
include("projectors.jl")

export symorb_minimize, plot_path 

check_convergence(res::Optim.OptimizationResults)::Bool = res.iteration_converged ||  res.x_converged || res.f_converged || res.g_converged
check_convergence(res::NLsolve.SolverResults)::Bool = res.x_converged || res.f_converged 

function symorb_minimize(config::AbstractDict, method=:BFGS; 
                    iter_max=Int(1e4),
                    starting_path_type = nothing,
                    starting_path = nothing, 
                    show_infos = false, 
                    plot_steps = false)

    LSG_from_config(config)

    if isnothing(starting_path_type) && isnothing(starting_path)
        Γ = random_starting_path()
    elseif starting_path_type == :circle
        Γ = circular_starting_path()
    else 
        Γ = starting_path
    end

    Γ = project(Γ)
    
    minimizer = nothing
    res = nothing
    res_options = MinimizationOptions()

    if method == :Newton    
        res = nlsolve(∇action, Haction, flatten(Γ), method=:trust_region, factor=10, autoscale=true, iterations=iter_max, show_trace = true)
        minimizer = (project ∘ emboss)(res.zero)
    else
        method_callable = getfield(Optim, method)   
        options = Optim.Options(g_tol = 1e-8, iterations = iter_max, show_trace = true)

        res = Optim.optimize(action, ∇action, Haction, flatten(Γ), method_callable(), options; inplace=false)
        minimizer = (project ∘ emboss)(res.minimizer)
        res_options = MinimizationOptions(options)
    end

    return MinimizationResult( 
        fourier_coeff   = minimizer,
        trajectory      = build_path(minimizer),
        iterations      = res.iterations,
        action_value    = _action(minimizer),
        gradient_norm   = norm(_∇action(minimizer)),
        initial         = Γ,
        converged       = check_convergence(res),
        method          = method,
        options         = res_options  
    );
end


function plot_path(path)
    pl = plot(aspect_ratio=:equal)
    for i ∈ 1:N
        pl = plot!(pl, [[path[j][i][h] for j ∈ 0:steps+1] for h ∈ 1:dim]..., label="Body $i", linewidth=3,  aspect_ratio=:equal);
    end
    display(pl)
end


end # module symorb
