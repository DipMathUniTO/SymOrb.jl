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

check_convergence(res::Optim.OptimizationResults)::Bool = (res.x_converged || res.f_converged || res.g_converged)
check_convergence(res::NLsolve.SolverResults)::Bool =  (res.x_converged || res.f_converged)


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


function find_critical_point_loop(Γ::Coefficients, methods::Methods; repetitions=10, perturb::Bool = false, perturbation::Float64 = 1e-3, action_threshold=2.0)::MinimizationResult
    result = nothing
    last_result = nothing
    for i ∈ 1:repetitions

        printstyled("#$i ", color=:yellow, bold=true)

        result = find_critical_point(Γ, methods.first)

        if ! isnothing(methods.second) && ! result.converged
            result = find_critical_point(result.fourier_coeff, methods.second)
        end

        println(result)

        if result.converged
            return result
        end

        if  result.action_value < action_threshold
            return result
        end

        println("Action value $(result.action_value) > threshold $(action_threshold)")
        
        if i < repetitions 
            if perturb
                printstyled("==> Perturbing path and continuing minimization\n\n", color=:yellow, bold=true)
                Γ = (project ∘ perturbed_path)(result.fourier_coeff, perturbation)
            else 
                printstyled("==> Continuing minimization\n\n", color=:yellow, bold=true)
                Γ = project(result.fourier_coeff)
            end
            last_result = result
        end
    end
    return result
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

    println(methods)

    while true 
        if isnothing(starting_path)
            Γ = get_starting_path(starting_path_type)
        else
            Γ = starting_path
        end
        Γ = project(Γ)

        printstyled("Starting a new minimization...\n\n", color=:blue, bold=true)

        if ! isnothing(methods.init)
            printstyled("Init: ", color=:yellow, bold=true)
            result = find_critical_point(Γ, methods.init)
            println(result)
        end
        
        result = find_critical_point_loop(Γ, methods)
        if result.converged && ! isnan(result.action_value) return result end
        printstyled("==> Path did not converge\n\n", color=:red, bold=true)
        
    end    

    return result
end

end
