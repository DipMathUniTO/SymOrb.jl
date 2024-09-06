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

function fourier_coefficients(P::Problem, x::Array{T, 3})::Vector{T} where T
    steps = size(x, 3)
    return P.x_to_A * reshape(x, P.dim*P.N*(steps), 1)
end

"""
    get_starting_path(path_type::Symbol)::OffsetArray

Generate the starting path for the minimization problem according to the given ``path_type``.
"""
function get_starting_path(path_type::Symbol, dimensions)::Vector
    if path_type == :circular
        return circular_starting_path(dimensions)
    elseif path_type == :perturbed_circular
        return perturbed_circular_starting_path(dimensions)
    else
        return random_starting_path(dimensions)
    end
end

"""
    random_starting_path()::OffsetArray

Generate a random starting path for the minimization problem.
"""
function random_starting_path(dimensions)::Vector
    F, N, dim = dimensions
    rand(Float64, dim * (N-1) *(F+2))
end


"""
    circular_starting_path()::OffsetArray

Generate a circular starting path for the minimization problem.
"""
function circular_starting_path(dimensions)::Coefficients
    dim, F, N = dimensions
    steps = 2 * F
    x = zeros(dim, N, steps+2)

    for h ∈ 0:steps+1, i ∈ 1:N-1
        x[1, i, h] = cos(h * π / (steps + 1) + (i - 1) * 2 * π / N)
        x[2, i, h] = sin(h * π / (steps + 1) + (i - 1) * 2 * π / N)
    end

    return FromZero(fourier_coefficients(x, F))
end

"""
    perturbed_path(x::Path, λ = 0.001)::Coefficients

Perturb the given path ``x`` by a factor ``λ``.
"""
perturbe_path(x, λ=0.001)::Coefficients = x + λ * random_starting_path(size(x))

"""
    perturbed_circular_starting_path(λ = 0.001)::Coefficients

Generate a perturbed circular starting path for the minimization problem.
"""
perturbed_circular_starting_path(dimensions, λ::Float64=0.001)::Coefficients = perturbe_path(circular_starting_path(dimensions), λ)


"""
    extend_to_period(x::Path)::Path

Extend the path ``x`` form the fundamental domain I = [0, π] to a full period.
"""
function extend_to_period(P::Problem, x::Array{T, 3})::Array{T, 3} where T<:Real

    dim, N, n = size(x)

    initial_path = zeros(dim, N, 2n-1)
    complete_path =  zeros(dim, N, (2n-1) * cyclic_order(P.G)) 

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

    # If the cyclic order is greater than 1, we need to add the images of the path under the group action
    for j ∈ 1:cyclic_order(P.G), k ∈ 1:max_index
        complete_path[:, :, max_index * (j-1) + k] =  reshape(ϕg_n(P.G.g)[j] *  initial_path[:, :, k][:], dim, N)
    end


    # Finally, we add the starting point of the path to close the loop
    complete_path[:, :, end] = reshape(ϕg_n(P.G.g)[1] *  initial_path[:, :, 1][:], dim, N) 

    return complete_path
end

"""
    extend_to_period(f::T) where {T <: AbstractVector{<:Real}}

Extend the function ``f`` form the fundamental domain I = [0, π] to a full period.
"""
function extend_to_period(P::Problem, f::T) where {T<:AbstractVector{<:Real}}
    n = lastindex(f)

    initial_f = zeros(2n - 1)
    complete_f = zeros((2n-1) * cyclic_order(P.G))

    initial_f[1:n] = copy(f[1:n])

    if P.G.action_type == Cyclic
        max_index = n
    else
        max_index = 2n-1
        for k ∈ 1:n-1
            initial_f[n+k] = initial_f[n-k]
        end
    end

    for j ∈ 1:cyclic_order(P.G), k ∈ 1:max_index
        complete_f[ max_index * (j-1) + k ] = initial_f[k]
    end

    complete_f[end] = initial_f[1]

    return complete_f
end


"""
    print_path_to_file(Γ::Coefficients, filename::String)

Print the Fourier coefficients of the path ``Γ`` to a file.
"""
function print_path_to_file(P::Problem, Γ::Vector, filename::String)
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

end

"""
    read_path_from_file(filename::String)::Coefficients

Read the Fourier coefficients of a path from a file.
"""
function read_path_from_file(P::Problem, filename::String)
    
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
