module MinPathJL

using GAP, OffsetArrays, Optim, LinearAlgebra, GLMakie, NLsolve, LineSearches, ForwardDiff

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
                    starting_path_type = :random,
                    starting_path = nothing,
                    show_trace = true, 
                    denominator = nothing)

    LSG_from_config(config)

    if isnothing(denominator)
        global f = x -> x
        global df = x -> 1
        global d2f = x -> 0
    else
        global f = denominator
        global df = x -> ForwardDiff.derivative(denominator, x)
        global d2f = x -> ForwardDiff.derivative(df, x)
    end

    if isnothing(starting_path)
        Γ = get_starting_path(starting_path_type)
    else 
        Γ = starting_path
    end

    Γ = project(Γ)
    
    minimizer = nothing
    res = nothing
    res_options = MinimizationOptions()

    println("Starting minimization...")
    if method == :Newton || method == :NewtonOnly
        while isnothing(res) || !check_convergence(res)
            if method != :NewtonOnly  
                options = Optim.Options(g_tol = 1e-8, allow_f_increases = true, iterations = 5, show_trace = show_trace)
                res = Optim.optimize(action, ∇action, Haction, flatten(Γ), BFGS(), options; inplace=false)
                v = res.minimizer
            else 
                v = flatten(Γ)
            end
        
            res = NLsolve.nlsolve(∇action, Haction, v, method=:trust_region, factor = 10, ftol = 1e-8, iterations=iter_max, show_trace = show_trace)
            println("\tAction:\t\t",action(res.zero)) 
            println("\tGradient:\t",norm(∇action(res.zero))) 
            print("\tConverged:\t")
            printstyled( check_convergence(res), "\n", color = check_convergence(res) ? :green : :red)
            println()
            if ! check_convergence(res)
                
                if action(res.zero) > 1.0 
                    method = :Newton
                    Γ = emboss(res.zero)
                    println("Action value is high, continuing minimization")
                else 
                    Γ = project(get_starting_path(starting_path_type))
                    println("Action value is low, starting a new minimization")
                end
            end 
        end

        minimizer = res.zero
    else
        method_callable = getfield(Optim, method)   
        options = Optim.Options(g_tol = 1e-8, iterations = iter_max, show_trace = show_trace)

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
