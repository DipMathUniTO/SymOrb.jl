"""
    build_path(P::Problem, Γ::Vector{T})::Array{T, 3}
    build_path(P::Problem, Γ::Vector{T}, steps::Int64)::Array{T, 3}

Build the path from the Fourier coefficients `Γ`. If `steps` is provided, the path is built with `steps` steps.
NOTE: If `steps` is provided, the function will compute the Jacobian of the path in terms of the Fourier coefficients.
Otherwise, it will use the precomputed Jacobian.
"""
function build_path(P::Problem, Γ::Vector{T})::Array{T, 3} where T
    steps = size(P.A_to_x, 1) ÷ (P.N*P.dim)
    x = P.A_to_x * Γ
    reshape(x, P.dim, P.N, steps)
end

function build_path(P::Problem, Γ::Vector{T}, steps::Int64)::Array{T, 3} where T
    F = size(Γ, 1) ÷ ((P.N-1)*P.dim) - 2
    A_to_x = kron(compute_dx_dA(F, steps-2), I(P.dim * P.N)) * P.R
    x = A_to_x * Γ
    reshape(x, P.dim, P.N, steps)
end

"""
    fourier_coefficients(P::Problem, x::Array{T, 3})::Vector{T}

Compute the Fourier coefficients of the path ``x``.
"""
function fourier_coefficients(P::Problem, x::Array{T, 3})::Vector{T} where T
    return P.x_to_A * x[:]
end

function fourier_coefficients(P::Problem, x::Array{T, 3}, steps::Int64)::Vector{T} where T
    x_to_A = P.Ri * kron(compute_dA_dx(P.F, steps-2), I(P.dim * P.N)) 
    return x_to_A * x[:]
end


"""
    get_starting_path(path_type::Symbol, dimensions)::Vector

Generate the starting path for the minimization problem according to the given ``path_type``.
"""
function get_starting_path(P::Problem, path_type::Symbol)::Vector
    if path_type == :circular
        return circular_starting_path(P::Problem)
    elseif path_type == :perturbed_circular
        return perturbed_circular_starting_path(P::Problem)
    else
        return random_starting_path(P::Problem)
    end
end

"""
    random_starting_path(dimensions::NTuple{3, Int64})::Vector

Generate a random starting path for the minimization problem.
"""
function random_starting_path(P::Problem)::Vector
    rand(Float64, P.dim * (P.N-1) *(P.F+2))
end


"""
    circular_starting_path(dimensions)::Vector

Generate a circular starting path for the minimization problem.
"""
function circular_starting_path(P::Problem)

    x = zeros(P.dim, P.N, P.steps+2)

    for h ∈ 0:P.steps+1, i ∈ 1:P.N
        x[1, i, h+1] = cos(h * π / (P.steps + 1) + (i - 1) * 2 * π / P.N)
        x[2, i, h+1] = sin(h * π / (P.steps + 1) + (i - 1) * 2 * π / P.N)
    end
    return fourier_coefficients(P, x)
end

"""
   perturbe_path(Γ::Vector{T}, dims::NTuple{3, Int64}, λ=0.001)::Vector{T}

Perturb the given path `Γ`` by a factor ``λ``.
"""
function perturbe_path(P::Problem, Γ::Vector{T}, λ=0.001)::Vector{T} where T 
    Γ + λ * random_starting_path(P)
end

"""
    perturbed_circular_starting_path(dimensions, λ::Float64=0.001)::Vector

Generate a perturbed circular starting path for the minimization problem.
"""
perturbed_circular_starting_path(P::Problem, λ::Float64=0.001)::Vector = perturbe_path(P, circular_starting_path(P::Problem), λ)


"""
    extend_to_period(P::Problem, x::Array{T, 3})::Array{T, 3} 

Extend the path ``x`` form the fundamental domain I = [0, π] to a full period.
"""
function extend_to_period(P::Problem, x::Array{T, 3})::Array{T, 3} where T<:Real

    dim, N, steps = size(x)
    n = steps - 1
    b = P.G.action_type == Cyclic ? 1 : 2
    
    complete_path = zeros(dim, N, b*n*cyclic_order(P.G) + 1) 
    complete_path[:, :, 1:steps] .= x

    # If the group is not cyclic, we first need to add the reflection of the path
    if P.G.action_type != Cyclic
        for k ∈ 1:n
            complete_path[:, :, steps + k] = ϕ(P.G.H1, x[:, :, steps - k])
        end
        n *= b
    end

    # If the cyclic order is greater than 1, we need to add the images of the path under the group action
    for j ∈ 2:cyclic_order(P.G), k ∈ 1:n
        complete_path[:, :, n * (j-1) + k] .=  ϕ(P.G.g[j-1], complete_path[:, :, k])
    end

    # Finally, we add the starting point of the path to close the loop
    complete_path[:, :, end] = complete_path[:, :, 1] 

    return complete_path
end

"""
    extend_to_period(P::Problem, f::T)::Array{T, 3}

Extend the function ``f`` form the fundamental domain I = [0, π] to a full period.
"""
function extend_to_period(P::Problem, f::Vector{T})::Vector{T} where {T}
    steps = length(f)
    n = steps - 1
    b = P.G.action_type == Cyclic ? 1 : 2
    
    complete_f = zeros(b*n*cyclic_order(P.G) + 1) 
    complete_f[1:steps] .= f

    # If the group is not cyclic, we first need to add the reflection of the function
    if P.G.action_type != Cyclic
        for k ∈ 1:n
            complete_f[steps + k] = f[steps - k]
        end
        n *= b
    end

    # If the cyclic order is greater than 1, we need to copy the function
    for j ∈ 2:cyclic_order(P.G), k ∈ 1:n
        complete_f[n * (j-1) + k] =  complete_f[k]
    end

    # Finally, we add the starting point
    complete_f[end] = complete_f[1] 

    return complete_f
end


"""
   print_path_to_file(P::Problem, Γ::Vector, filename::String)::Nothing

Print the information about `P` and the Fourier coefficients of the path ``Γ`` to the file named `filename`.
"""
function print_path_to_file(P, Γ, filename)
    data = copy(P.meta)
    data["path"] = Γ
    open(filename, "w") do io
        TOML.print(io, data)
    end
end

"""
    read_path_from_file(filename::String)::Tuple{Problem, Vector}

Read the Fourier coefficients of a path  and the problem configuration from the file named `filename`.
"""
function read_path_from_file(filename::String)::Tuple{Problem, Vector}
    data = TOML.parsefile(filename)
    P = initialize(data)
    return P, data["path"]
end

function path_animation()
    error("GLMakie must be loaded")
end

function SymOrb.path_animation(filename::String, opts...)
    P, Γ = SymOrb.read_path_from_file(filename)
    SymOrb.path_animation(P, Γ; opts...)
    return P, Γ
end