
function log_if(variable, expr...; color=:normal, bold=false)
    if variable printstyled(expr..., color=color, bold=bold) end
end

has_converged(res::Optim.OptimizationResults)::Bool = (res.x_converged || res.f_converged || res.g_converged)
has_converged(res::NLsolve.SolverResults)::Bool = (res.x_converged || res.f_converged)

""" 
    perform_optimization(P::Problem, Γ0::Vector, method::NLSolveMethod)::Tuple{Vector, Int, Bool}
    perform_optimization(P::Problem, Γ0::Vector, method::OptimMethod)::Tuple{Vector, Int, Bool}

Wrappers around the minimization libraries
""" 
function perform_optimization(P::Problem, Γ0::Vector, method::NLSolveMethod)::Tuple{Vector, Int, Bool}

    res = nlsolve(x -> ∇action(P, x), x -> Haction(P, x), Γ0, method=method.f_name, iterations=method.max_iter, show_trace=method.show_trace, factor=10, ftol=1e-8)
    return res.zero, res.iterations, has_converged(res)
end

function perform_optimization(P::Problem, Γ0::Vector, method::OptimMethod)::Tuple{Vector, Int, Bool}
    method_callable = getfield(Optim, method.f_name)
    
    options = Options(g_tol=1e-8, iterations=method.max_iter, show_trace=method.show_trace)
    res = optimize(x-> action(P, x), x -> ∇action(P, x), x -> Haction(P, x), Γ0, method_callable(), options; inplace=false)
    return res.minimizer, res.iterations, has_converged(res)
end

function show_action(show_steps, Γ, p)
    log_if(show_steps, "     Action  = $(action(p, Γ))\n", bold = true)
    log_if(show_steps, "     ∇action = $(norm(∇action(p, Γ)))\n\n", bold = true)
end

"""
    new_orbit(p::Problem, starting_path::Vector, method::CompoundMethod; max_repetitions::Int=10, perturb::Bool=false, perturbation::Float64=1e-3, action_threshold::Float64=2.0, show_steps=true)::MinimizationResult

Find a new orbit using the given `CompoundMethod` and starting path. 

It first executes the `init` method. Then it alternates the `first` method and the `second` method
until the minimization converges, the action is below `action_threshold` or the maximum number of repetitions is reached.

To implement new optimization algorithms, define new methods for this function. 
New methods must accept 2 arguments, the starting path and the method to use, and a keyword argument `show_steps`.
They must return a `MinimizationResult`.

# Keyword Arguments
- `max_repetitions::Int=10`: how many times to minimize using first and second method
- `perturb::Bool=false`: whether to perturb the path if minimization fails 
- `perturbation::Float64=1e-3`: the strength of the perturbation to apply to the path
- `action_threshold::Float64=2.0`: the threshold for the action value above which to continue minimization
- `show_steps::Bool=true`: whether to show the steps of the minimization
"""
function new_orbit(p, starting_path::Vector, method::CompoundMethod; max_repetitions::Int=10, perturb::Bool=false, perturbation::Float64=1e-3, action_threshold::Float64=2.0, show_steps=true)::MinimizationResult
    Γ = p.Π*starting_path[:]

    converged = false
    # Perform the init steps of minimization to get closer to a minimum
    # This is required because second order methdos (Newton) are very unstable when far away from the minimum
    if !isnothing(method.init)
        log_if(show_steps, "  Init: ", color=:yellow, bold=true)
        log_if(show_steps, method.init, "\n")
        Γ, _, _ = perform_optimization(p, Γ, method.init)
        show_action(show_steps, Γ, p)
    end

    # Start alternating between the first and second method
    for i ∈ 1:max_repetitions
        log_if(show_steps, "  #$i \n", color = :yellow, bold = true)

        # Try with the firsy optimization method (usually, first order method)
        log_if(show_steps,  "     First method: ", color=:yellow, bold=true)
        log_if(show_steps, method.first, "\n")
        Γ, _, converged = perform_optimization(p, Γ, method.first)
        show_action(show_steps, Γ, p)

        # If the method converged or the action is below the threshold, we don't need to continue
        if converged || action(p, Γ) < action_threshold
            break
        end

        # If a second method is provided (usually, second order method), try with it
        if !isnothing(method.second)
            log_if(show_steps, "     Second method: "; color=:yellow, bold=true)
            log_if(show_steps, method.second, "\n")
            Γ, _, converged = perform_optimization(p, Γ, method.second)
            show_action(show_steps, Γ, p)
        end

        # If the method converged or the action is below the threshold, we don't need to continue
        if converged || action(p, Γ) < action_threshold
            break
        end

        # If there are still repetitions left, prepare to continue
        if i < max_repetitions
            if perturb
                Γ = perturbe_path(Γ, (P.F, P.N, P.dim), perturbation)
            end
            
            Γ = p.Π* Γ
            log_if(show_steps, "  ==> Continuing minimization\n\n", color = :yellow, bold = true)
        end
    end

    MinimizationResult(starting_path, p.Π*Γ, norm(∇action(p, Γ)), action(p, Γ), converged && action(p, Γ) > action_threshold, method)
end

""" 
    find_orbits(P::Problem, method::AbstractMethod=OneMethod(BFGS());; starting_path_type::Symbol=:random, starting_path::Union{Path,Nothing}=nothing, show_steps=true, print_path=true,  print_path_folder::String="/.", options...)

Find `number_of_orbits` orbits using the given `method` and starting path. 

# Arguments
- `method::AbstractMethod=OneMethod(BFGS())`: the minimization method to use. Defaults to BFGS for 200 steps.
- `number_of_orbits=Inf`: the number of orbits to find. Defaults to infinity.

# Keyword Arguments
- `starting_path_type::Symbol=:random`: the type of starting path to use
- `starting_path::Union{Path,Nothing}=nothing`: the starting path to use. If provided, `starting_path_type` is ignored
- `show_steps::Bool=true`: whether to show the steps of the minimization.
- `print_path::Bool=true`: whether to print the path to a file.
- `print_path_folder::String="./"`: the folder where to print the path.
- `options...`: additional options to pass to the optimization method.
"""
function find_orbits(P::Problem, method::AbstractMethod=OneMethod(BFGS()); number_of_orbits::Union{Int, Float64}=1, starting_path_type::Symbol=:random, starting_path::Union{Vector,Nothing}=nothing, show_steps=true, print_path=true, print_path_folder::String="./", options...)
    # Set the correct starting path type if the starting path is user-provided
    if ! isnothing(starting_path)
        starting_path_type = :given
    end

    if haskey(P.meta, "symmetry_name")
        print_path_folder = joinpath(print_path_folder, P.meta["symmetry_name"])*"/"
        if !isdir(print_path_folder)
            mkdir(print_path_folder)
        end
    end

    i = 0
    results = []
    @show number_of_orbits

    
    while i < number_of_orbits
        log_if(show_steps, "#$i: Searching new orbit...\n", color = :blue, bold = true)
        # If the starting path is not user-provided, generate it
        if starting_path_type != :given
            starting_path = get_starting_path(P, starting_path_type)
        end

        # Find a new orbit
        result = new_orbit(P, starting_path, method; show_steps=show_steps, options...)
        log_if(show_steps, result)

        # Check if the minimization converged and, if so, add the result to the list
        if result.converged && !isnan(result.action_value) && print_path
            # Print the path to a file if required
            if print_path
                print_path_to_file(P, result.fourier_coeff, @sprintf "%s%.4f.toml" print_path_folder result.action_value)
            end
            log_if(show_steps, "==> Optimization converged\n\n", color = :green, bold = true)
            push!(results, result)
        else
            log_if(show_steps, "==> Optimization did not converge\n\n", color = :red, bold = true)
        end

        i += 1
    end
    return results
end

