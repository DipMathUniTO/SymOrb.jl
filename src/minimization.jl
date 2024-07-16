
function if_log(variable, expr...; color=:normal, bold=false)
    if variable printstyled(expr..., color=color, bold=bold) end
end

has_converged(res::Optim.OptimizationResults)::Bool = (res.x_converged || res.f_converged || res.g_converged)
has_converged(res::NLsolve.SolverResults)::Bool = (res.x_converged || res.f_converged)

""" 
    perform_optimization(Γ0::Coefficients, method::NLSolveMethod)::Tuple{Coefficients, Int64, Bool}
    perform_optimization(Γ0::Coefficients, method::OptimMethod)::Tuple{Coefficients, Int64, Bool}

Wrappers around the minimization libraries
""" 
function perform_optimization(problem::Problem, Γ0::Coefficients, method::NLSolveMethod)::Tuple{Coefficients, Int64, Bool}
    dimensions = dims(Γ0)
    ∇func = x -> emboss(x, dimensions) |> (x -> ∇action(problem, x)) |> flatten
    Hfunc = x -> emboss(x, dimensions) |> (x -> Haction(problem, x)) |> flatten
    res = nlsolve(∇func, Hfunc, flatten(Γ0), method=method.f_name, iterations=method.max_iter, show_trace=method.show_trace, factor=10, ftol=1e-8)
    return emboss(res.zero, dimensions), res.iterations, has_converged(res)
end

function perform_optimization(problem::Problem, Γ0::Coefficients, method::OptimMethod)::Tuple{Coefficients, Int64, Bool}
    method_callable = getfield(Optim, method.f_name)
    dimensions = dims(Γ0)
    func  = x -> emboss(x, dimensions) |> (x -> Haction(problem, x))
    ∇func = x -> emboss(x, dimensions) |> (x -> ∇action(problem, x)) |> flatten
    Hfunc = x -> emboss(x, dimensions) |> (x -> Haction(problem, x)) |> flatten
    options = Options(g_tol=1e-8, iterations=method.max_iter, show_trace=method.show_trace)
    res = optimize(func, ∇func, Hfunc, flatten(Γ0), method_callable(), options; inplace=false)
    return emboss(res.minimizer, dimensions), res.iterations, has_converged(res)
end

function show_action(show_steps, Γ)
    if_log(show_steps, "     Action  = $(action(Γ))\n", bold = true)
    if_log(show_steps, "     ∇action = $(norm(∇action(Γ)))\n\n", bold = true)
end

"""
    new_orbit(starting_path::Coefficients, method::CompoundMethod; max_repetitions::Int64=10, perturb::Bool=false, perturbation::Float64=1e-3, action_threshold::Float64=2.0, show_steps=true)::MinimizationResult

Find a new orbit using the given `CompoundMethod` `method` and starting path. 

It first executes the `init` method. Then it alternates the `first` method and the `second` method
until the minimization converges, the action is below `action_threshold` or the maximum number of repetitions is reached.

To implement new optimization methods, define new methods for this function. 
New methods must accept 2 arguments, the starting path and the method to use, and a keyword argument `show_steps`.
They must return a `MinimizationResult`.

# Keyword Arguments
- `max_repetitions::Int64=10`: how many times to minimize using first and second method
- `perturb::Bool=false`: whether to perturb the path if minimization fails 
- `perturbation::Float64=1e-3`: the strength of the perturbation to apply to the path
- `action_threshold::Float64=2.0`: the threshold for the action value above which to continue minimization
- `show_steps::Bool=true`: whether to show the steps of the minimization
"""
function new_orbit(problem, starting_path::Coefficients, method::CompoundMethod; max_repetitions::Int64=10, perturb::Bool=false, perturbation::Float64=1e-3, action_threshold::Float64=2.0, show_steps=true)::MinimizationResult
    Γ = project(starting_path[:])
    converged = false
    # Perform the init steps of minimization to get closer to a minimum
    # This is required because second order methdos (Newton) are very unstable when far away from the minimum
    if !isnothing(method.init)
        if_log(show_steps, "  Init: "; color=:yellow, bold=true)
        if_log(show_steps, method.init, "\n")
        Γ, _, _ = perform_optimization(problem, Γ, method.init)
        show_action(show_steps, Γ)
    end

    # Start alternating between the first and second method
    for i ∈ 1:max_repetitions
        if_log(show_steps, "  #$i \n", color = :yellow, bold = true)

        # Try with the firsy optimization method (usually, first order method)
        if_log(show_steps, "     First method: "; color=:yellow, bold=true)
        if_log(show_steps, method.first, "\n")
        Γ, _, converged = perform_optimization(problem, Γ, method.first)
        show_action(show_steps, Γ)

        # If the method converged or the action is below the threshold, we don't need to continue
        if converged || action(Γ) < action_threshold
            break
        end

        # If a second method is provided (usually, second order method), try with it
        if !isnothing(method.second)
            if_log(show_steps, "     Second method: "; color=:yellow, bold=true)
            if_log(show_steps, method.second, "\n")
            Γ, _, converged = perform_optimization(problem, Γ, method.second)
            show_action(show_steps, Γ)
        end

        # If the method converged or the action is below the threshold, we don't need to continue
        if converged || action(Γ) < action_threshold
            break
        end

        # If there are still repetitions left, prepare to continue
        if i < max_repetitions
            if perturb
                Γ = perturbed_path(Γ, perturbation)
            end
            
            Γ = project(Γ)
            if_log(show_steps, "  ==> Continuing minimization\n\n", color = :yellow, bold = true)
        end
    end

    MinimizationResult(starting_path, project(Γ), method, converged)
end

""" 
    find_orbits(method::AbstractMethod=OneMethod(BFGS()), number_of_orbits::Float64=Inf; starting_path_type::Symbol=:random, starting_path::Union{Path,Nothing}=nothing, show_steps=true, print_path=true, options...)

Find `number_of_orbits` orbits using the given `method` and starting path. 

# Arguments
- `method::AbstractMethod=OneMethod(BFGS())`: the minimization method to use
- `number_of_orbits::Float64=Inf`: the number of orbits to find

# Keyword Arguments
- `starting_path_type::Symbol=:random`: the type of starting path to use
- `starting_path::Union{Path,Nothing}=nothing`: the starting path to use
- `show_steps::Bool=true`: whether to show the steps of the minimization
- `print_path::Bool=true`: whether to print the path to a file
- `options...`: additional options to pass to the optimization method
"""
function find_orbits(problem::Problem, method::AbstractMethod=OneMethod(BFGS()), number_of_orbits::Int64=Inf; starting_path_type::Symbol=:random, starting_path::Union{Path,Nothing}=nothing, show_steps=true, print_path=true, options...)
    # Set the correct starting path type if the starting path is user-provided
    if ! isnothing(starting_path)
        starting_path_type = :given
    end

    i = 0
    results = []
    
    while i < number_of_orbits
        if_log(show_steps, "#$i: Searching new orbit...\n", color = :blue, bold = true)
        # If the starting path is not user-provided, generate it
        if starting_path_type != :given
            starting_path = get_starting_path(starting_path_type, problem.dims)
        end


        # Find a new orbit
        result = new_orbit(problem, starting_path, method; show_steps=show_steps, options...)
        if_log(show_steps, result)

        # Check if the minimization converged and, if so, add the result to the list
        if result.converged && !isnan(result.action_value) && print_path
            # Print the path to a file if required
            if print_path
                print_path_to_file(result.fourier_coeff, @sprintf "%.5f.txt" result.action_value)
            end
            if_log(show_steps, "==> Optimization converged\n\n", color = :green, bold = true)
            push!(results, result)
        else
            if_log(show_steps, "==> Optimization did not converge\n\n", color = :red, bold = true)
        end

        i += 1
    end

    return results
end

