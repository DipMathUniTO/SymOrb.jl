module MinPathJL

using GAP, OffsetArrays, Optim, LinearAlgebra, GLMakie, NLsolve, LineSearches, ForwardDiff

include("globals.jl")
include("path.jl")
include("initdata.jl")
include("matrices.jl")
include("action.jl")
include("projectors.jl")

export find_orbit, plot_path, print_path_to_file, path_animation
export Newton, BFGS, ConjugateGradient, Methods

check_convergence(res::Optim.OptimizationResults)::Bool = res.x_converged || res.f_converged || res.g_converged
check_convergence(res::NLsolve.SolverResults)::Bool = res.x_converged || res.f_converged


function find_critical_point(Γ::Coefficients, method::NLSolveMethod)
    res = NLsolve.nlsolve(∇action, Haction, flatten(Γ), method=method.f_name, iterations=method.max_iter, show_trace=method.show_trace, factor=10, ftol=1e-8)
    MinimizationResult(Γ, res.zero, method, res.iterations, check_convergence(res))
end

function find_critical_point(Γ::Coefficients, method::OptimMethod)
    method_callable = getfield(Optim, method.f_name)
    options = Optim.Options(g_tol=1e-8, iterations=method.max_iter, show_trace=method.show_trace)
    res = Optim.optimize(action, ∇action, Haction, flatten(Γ), method_callable(), options; inplace=false)
    MinimizationResult(Γ, res.minimizer, method, res.iterations, check_convergence(res))
end

function find_orbit(config::AbstractDict, methods=[(:BFGS, 200)]; starting_path_type=:random, starting_path=nothing, denominator=nothing, options...)

    LSG_from_config(config)
    max_repetitions = 100
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

    println("Starting minimization...")

    result = nothing
    
    
    @show methods
    if ! isnothing(methods.init)
        result = find_critical_point(Γ, methods.init)
    end
    
    repetitions = 0
    while (isnothing(result) || ! result.converged) && repetitions < max_repetitions
        println(result)
        
        method = if iseven(repetitions) && ! isnothing(methods.second) methods.second else methods.first end

        result = find_critical_point(Γ, method)

        if !result.converged
            println(result)
            if result.action_value > 1.0
                println("Action value is high\nPerturbing path and continuing minimization")
                method = :NewtonOnly
                Γ = perturbed_path(result.fourier_coeff, 5e-3)
            else
                println("Action value is low, starting a new minimization")
                Γ = project(get_starting_path(starting_path_type))
            end
        end
        repetitions += 1
    end

    return result
end

end
