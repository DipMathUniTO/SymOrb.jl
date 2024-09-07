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

    dim, N, n = size(x)

    initial_path = zeros(dim, N, 2n-1)
    
    # The first half of the path is the same as the original path
    initial_path[:, :, 1:n] = copy(x[:, :, 1:n])
    
    
    # The second half of the path is the reflection of the first half if the action is not Cyclic
    if P.G.action_type == Cyclic
        max_index = n
    else
        max_index = 2n - 1 
        for k ∈ 1:n-1
            initial_path[:, :, n+k] = reshape(ϕ(P.G.H1) * initial_path[:, :, n-k][:], dim, N)
        end
    end

    complete_path =  zeros(dim, N, max_index * cyclic_order(P.G)) 
    
    # If the cyclic order is greater than 1, we need to add the images of the path under the group action
    for j ∈ 1:cyclic_order(P.G), k ∈ 1:max_index
        complete_path[:, :, max_index * (j-1) + k] =  reshape(ϕg_n(P.G.g)[j] *  initial_path[:, :, k][:], dim, N)
    end


    # Finally, we add the starting point of the path to close the loop
    complete_path[:, :, end] = reshape(ϕg_n(P.G.g)[1] *  initial_path[:, :, 1][:], dim, N) 

    return complete_path
end

"""
    extend_to_period(P::Problem, f::T)::Array{T, 3}

Extend the function ``f`` form the fundamental domain I = [0, π] to a full period.
"""
function extend_to_period(P::Problem, f::Vector{T})::Vector{T} where {T}
    n = lastindex(f)

    initial_f = zeros(2n - 1)
    
    initial_f[1:n] = copy(f[1:n])
    
    if P.G.action_type == Cyclic
        max_index = n
    else
        max_index = 2n-1
        for k ∈ 1:n-1
            initial_f[n+k] = initial_f[n-k]
        end
    end
    complete_f = zeros( max_index * cyclic_order(P.G))

    for j ∈ 1:cyclic_order(P.G), k ∈ 1:max_index
        complete_f[ max_index * (j-1) + k ] = initial_f[k]
    end

    complete_f[end] = initial_f[1]
    @show typeof(complete_f)
    return complete_f
end


"""
   print_path_to_file(P::Problem, Γ::Vector, filename::String)::Nothing

Print the Fourier coefficients of the path ``Γ`` to the file named `filename`.
"""
function print_path_to_file(P::Problem, Γ::Vector, filename::String)::Nothing
    Γr = P.R * Γ
    
    F = length(Γr) ÷ (P.N*P.dim)

    Γr = reshape(Γr, P.dim, P.N, F)

    data = Matrix{Float64}(undef, F*P.N, P.dim)

    for k ∈ axes(data, 1)7
        h = (k - 1) % F + 1
        i = (k - 1) ÷ F + 1
        data[k, :] = Γr[:, i, h]
    end

    writedlm(filename, data, ' ')
    nothing
end

"""
    read_path_from_file(P::Problem, filename::String)::Vector

Read the Fourier coefficients of a path from the file named `filename`.
"""
function read_path_from_file(P::Problem, filename::String)::Vector
    
    data = readdlm(filename, ' ')
    F = size(data, 1) ÷ P.N
    Γ = zeros(P.dim, P.N, 26)

    for k in axes(data, 1)
        h = (k - 1) % F + 1
        i = (k - 1) ÷ F + 1
        Γ[:, i, h] = data[k, :]
    end
    return P.Ri * Γ[:]
end

function path_animation()
    error("GLMakie must be loaded")
end
