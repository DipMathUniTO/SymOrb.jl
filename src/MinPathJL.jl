module MinPathJL

using GAP, OffsetArrays, Optim, LinearAlgebra, GLMakie, NLsolve, LineSearches

include("globals.jl")
include("path.jl")
include("initdata.jl")
include("matrices.jl")
include("action.jl")
include("projectors.jl")

export minimize, plot_path, print_path_to_file, path_animation


check_convergence(res::Optim.OptimizationResults)::Bool = res.iteration_converged ||  res.x_converged || res.f_converged || res.g_converged
check_convergence(res::NLsolve.SolverResults)::Bool = res.x_converged || res.f_converged 

function minimize(config::AbstractDict, method=:BFGS; 
                    iter_max=Int(1e4),
                    starting_path_type = nothing,
                    starting_path = nothing, 
                    show_infos = false, 
                    plot_steps = false, ϵ_in=0)

    LSG_from_config(config)
    global ϵ = ϵ_in
    
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

    if method == :Newton || method == :NewtonOnly
        while isnothing(res) || !check_convergence(res)
            Γ = project(random_starting_path())
            if method != :NewtonOnly  
                options = Optim.Options(g_tol = 1e-8, allow_f_increases = true, iterations = 10, show_trace = true)
                res = Optim.optimize(action, ∇action, Haction, flatten(Γ), BFGS(), options; inplace=false)
                v = res.minimizer
            else 
                v = flatten(Γ)
            end
        
            res = NLsolve.nlsolve(∇action, Haction, v, method=:trust_region, factor = 5, ftol = 1e-8, iterations=iter_max, show_trace = true)
        end

        minimizer = res.zero
    else
        method_callable = getfield(Optim, method)   
        options = Optim.Options(g_tol = 1e-8, iterations = iter_max, show_trace = true)

        res = Optim.optimize(action, ∇action, Haction, flatten(Γ), method_callable(), options; inplace=false)
        minimizer = res.minimizer
        res_options = MinimizationOptions(options)
    end
    Γ_min =  (project ∘ emboss)(minimizer)

    result = MinimizationResult( 
        initial         = Γ,
        fourier_coeff   = Γ_min,
        path            = reconstruct_path(build_path(Γ_min)), 
        iterations      = res.iterations,
        action_value    = action(minimizer),
        gradient_norm   = norm(∇action(minimizer)),
        converged       = check_convergence(res),
        method          = method,
        options         = res_options  
    );
    println(result)
    return result
end

end
