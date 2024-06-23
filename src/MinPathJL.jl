module MinPathJL

using GAP
using OffsetArrays, Reexport, JSON
@reexport using LinearAlgebra
using Optim, NLsolve, LineSearches, ForwardDiff
using GLMakie, LaTeXStrings

include("globals.jl")
include("path.jl")
include("initialize.jl")
include("graphics.jl")
include("matrices.jl")
include("action.jl")
include("projectors.jl")

export find_orbits, initialize_and_find_orbits, initialize
export print_path_to_file, path_animation, refine_path, read_path_from_file
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

function find_critical_point_loop(Γ::Coefficients, methods::Methods; max_repetitions::Int64=10, perturb::Bool = false, perturbation::Float64 = 1e-3, action_threshold::Float64=2.0)::MinimizationResult
    result = nothing
    for i ∈ 1:max_repetitions

        printstyled("#$i ", color=:yellow, bold=true)

        result = find_critical_point(Γ, methods.first)

        if ! isnothing(methods.second) && ! result.converged
            result = find_critical_point(result.fourier_coeff, methods.second)
        end

        println(result)
 
        if result.converged || result.action_value < action_threshold
            return result 
        end

        println("Action value $(result.action_value) > threshold $(action_threshold)")
        
        if i < max_repetitions 
            if perturb
                printstyled("==> Perturbing path and continuing minimization\n\n", color=:yellow, bold=true)
                Γ = (project ∘ perturbed_path)(result.fourier_coeff, perturbation)
            else 
                printstyled("==> Continuing minimization\n\n", color=:yellow, bold=true)
                Γ = project(result.fourier_coeff)
            end
        end
    end
    return result
end 


function find_orbits(methods::Methods=Methods(BFGS()); starting_path_type::Symbol=:random, starting_path::Union{Path, Nothing}=nothing, options...)

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
            result = find_critical_point(Γ, methods.init; options...)
            println(result)
        end
        
        result = find_critical_point_loop(Γ, methods; options...)
        if (result.converged  || result.gradient_norm < 1e-2) && ! isnan(result.action_value) 
            print_path_to_file(result.fourier_coeff, string(result.action_value)*".txt")
        end
        printstyled("==> Optimization did not converge\n\n", color=:red, bold=true)
        
    end    

    return result
end

function initialize_and_find_orbits(file::String, methods::Methods=Methods(BFGS()), options...)::MinimizationResult
    initialize(file)
    find_orbits(methods; options...)
end



function minimize_lattice(f_opt, f_starting, lb, ub, n_div)
    x0 = f_starting(lb, ub)
    mins = f_opt(x0, lb, ub)

    if n_div == 0
        return mins
    end
    
    cutoff, index  = lat.lattice(mins)

    new_ub, new_lb = lat.split(lb, ub, index, cutoff)

    lower_mins = minimize_lattice(f, f_starting, lb, new_ub, n_div-1)
    upper_mins = minimize_lattice(f, f_starting, new_lb, ub, n_div-1)

    return vcat(lower_mins, upper_mins)
end

end