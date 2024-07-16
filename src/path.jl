"""
    fourier_series(A::Coefficients, n=steps+1)::Path

Compute the Fourier series for a given configuration ``A`` and using ``n`` points along the path.
"""
fourier_series(A::Coefficients, n)::Path = FromZero([sum(A[k] * sin(π * k * h / n) for k ∈ 1:lastindex(A)-1) for h ∈ 0:n])

"""
    inverse_fourier(x::Path, n=steps+1)::Coefficients

Compute the inverse Fourier series for a given path ``x`` obtaining ``F`` coefficients.
"""
inverse_fourier(x::Path, F)::Coefficients = 2 / lastindex(x) * FromZero([sum(x[h] * sin(k * h * π / lastindex(x)) for h ∈ 0:lastindex(x)) for k ∈ 0:F+1])

"""
    segment(a::Vector{T}, b::Vector{T}, n=steps+1)::Vector{T}

Compute the segment between two points ``a`` and ``b`` using ``n`` points.
"""
segment(a, b, n) = FromZero([a + k * (b - a) / n for k ∈ 0:n])

"""
    build_path(A::Coefficients, n=steps+1)::Path

Build a path from the Fourier coefficients ``A`` using ``n`` points.
"""
function build_path(A::Coefficients, n::Int)::Path
    return segment(A[0], A[end], n) + fourier_series(A, n)
end

"""
    fourier_coefficients(y::Path, n=steps+1)::Coefficients

Compute ``F`` Fourier coefficients for a given path ``y``.
"""
function fourier_coefficients(y::Path, F::Int)::Coefficients
    n = lastindex(y)
    x = y - segment(y[0], y[end], n)
    A = inverse_fourier(x, F)
    A[0] = y[0]
    A[end] = y[end]
    return A
end


"""
    emboss(v::Vector{T})::Coefficients

Convert a 1D vector ``v`` into a nested vector.
"""
function emboss(v::Vector{T}, dimensions)::Coefficients where {T}
    F, N, dim = dimensions
    Γ = FromZero([[zeros(T, dim) for i ∈ 1:N] for k ∈ 0:F+1])

    for j ∈ 1:dim, i ∈ 1:N, k ∈ 0:F+1
        Γ[k][i][j] = v[k*N*dim+(i-1)*dim+j]
    end

    return Γ
end

"""
    emboss(M::Matrix{T})::OffsetArray

Convert a 2D matrix ``M`` into a nested matrix.
"""
function emboss(M::Matrix{T}, dimensions)::OffsetArray where {T}
    F, N, dim = dimensions   
    H = FromZero([[zeros(T, dim, dim) for _ ∈ 1:N, _ ∈ 1:N] for _ ∈ 0:F+1, _ ∈ 0:F+1])

    for k1 ∈ 0:F+1, k2 ∈ 0:F+1, i1 ∈ 1:N, i2 ∈ 1:N, j1 ∈ 1:dim, j2 ∈ 1:dim
        H[k1, k2][i1, i2][j1, j2] = M[k1*N*dim+(i1-1)*dim+j1, k2*N*dim+(i2-1)*dim+j2]
    end

    return H
end

"""
    flatten(Γ::AbstractVector{Vector{Vector{T}}})::Vector{T}

Convert a nested vector ``Γ`` into a 1D vector.
"""
function flatten(Γ::AbstractVector{Vector{Vector{T}}})::Vector{T} where {T}
    F, N, dim = dims(Γ)
    v = zeros(T, (F + 2) * N * dim)

    for k ∈ 0:F+1, i ∈ 1:N, j ∈ 1:dim
        v[k*N*dim+(i-1)*dim+j] = Γ[k][i][j]
    end

    return v
end


"""
    flatten(H::AbstractMatrix{Matrix{Matrix{T}}})::Matrix{T}

Convert a nested matrix ``H`` into a 2D matrix.
"""
function flatten(H::AbstractMatrix{Matrix{Matrix{T}}})::Matrix{T} where {T}
    F, N, dim = dims(H)
    M = zeros(T, (F + 2) * N * dim, (F + 2) * N * dim)

    for k1 ∈ 0:F+1, k2 ∈ 0:F+1, i1 ∈ 1:N, i2 ∈ 1:N, j1 ∈ 1:dim, j2 ∈ 1:dim
        M[k1*N*dim+(i1-1)*dim+j1, k2*N*dim+(i2-1)*dim+j2] = H[k1, k2][i1, i2][j1, j2]
    end

    return M
end

"""
    get_starting_path(path_type::Symbol)::OffsetArray

Generate the starting path for the minimization problem according to the given ``path_type``.
"""
function get_starting_path(path_type::Symbol, dimensions)::OffsetArray
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
function random_starting_path(dimensions)::OffsetArray
    F, N, dim = dimensions
    FromZero([[rand(Float64, dim) .- 0.5 for i ∈ 1:N] for j ∈ 0:F+1])
end


"""
    circular_starting_path()::OffsetArray

Generate a circular starting path for the minimization problem.
"""
function circular_starting_path(dimensions)::Coefficients
    F, N, dim = dimensions
    steps = 2 * F
    x::Path = FromZero([[zeros(dim) for i ∈ 1:N] for j ∈ 0:steps+1])

    for h ∈ 0:steps+1, i ∈ 1:N
        x[h][i][1] = cos(h * π / (steps + 1) + (i - 1) * 2 * π / N)
        x[h][i][2] = sin(h * π / (steps + 1) + (i - 1) * 2 * π / N)
    end

    return FromZero(fourier_coefficients(x, F))
end

"""
    perturbed_path(x::Path, λ = 0.001)::Coefficients

Perturb the given path ``x`` by a factor ``λ``.
"""
perturbe_path(x::Path, λ=0.001)::Coefficients = x + λ * random_starting_path(dims(x))

"""
    perturbed_circular_starting_path(λ = 0.001)::Coefficients

Generate a perturbed circular starting path for the minimization problem.
"""
perturbed_circular_starting_path(dimensions, λ::Float64=0.001)::Coefficients = perturbe_path(circular_starting_path(dimensions), λ)


"""
    extend_to_period(x::Path)::Path

Extend the path ``x`` form the fundamental domain I = [0, π] to a full period.
"""
function extend_to_period(x::Path)::Path

    n = lastindex(x) - 1

    initial_path = FromZero([[zeros(dim) for i ∈ 1:N] for j ∈ 0:(2*n+1)])
    complete_path = FromZero([])

    # The first half of the path is the same as the original path
    initial_path[0:n+1] = copy(x[0:n+1])

    # The second half of the path is the reflection of the first half if the action is not Cyclic
    if G.action_type == Cyclic
        max_index = n + 1
    else
        max_index = 2 * n + 1
        for k ∈ 0:n
            initial_path[n+k+1] = ϕ(G.H1, initial_path[n-k+1])
        end
    end

    # If the cyclic order `m = length(P.g)` is greater than 1, we need to add the images of the path under the group action
    for j ∈ 1:length(G.g), k ∈ 0:max_index
        push!(complete_path, ϕg_n(initial_path[k], j))
    end

    # Finally, we add the starting point of the path to close the loop
    push!(complete_path, ϕg_n(initial_path[0], 1))

    return complete_path
end

"""
    extend_to_period(f::T) where {T <: AbstractVector{<:Real}}

Extend the function ``f`` form the fundamental domain I = [0, π] to a full period.
"""
function extend_to_period(f::T) where {T<:AbstractVector{<:Real}}
    n = lastindex(f) - 1

    initial_f = FromZero(zeros(2 * n + 2))
    complete_f = FromZero([])
    initial_f[0:n+1] = copy(f[0:n+1])

    if G.action_type == Cyclic
        max_index = n + 1
    else
        max_index = 2 * n + 1
        for k ∈ 0:n
            initial_f[n+k+1] = initial_f[n-k+1]
        end
    end

    for _ ∈ 1:length(G.g), k ∈ 0:max_index
        push!(complete_f, initial_f[k])
    end

    push!(complete_f, initial_f[0])

    return complete_f
end


"""
    print_path_to_file(Γ::Coefficients, filename::String)

Print the Fourier coefficients of the path ``Γ`` to a file.
"""
function print_path_to_file(Γ::Coefficients, filename::String)
    F, N, dim = dims(Γ)

    data = Matrix{Float64}(undef, (length(Γ[0]) * (F + 2)), dim)

    for k ∈ axes(data, 1)
        h = (k - 1) % (F + 2)
        i = div((k - 1), F + 2) + 1
        data[k, :] = Γ[h][i]
    end

    writedlm(filename, data, ' ')

end

"""
    read_path_from_file(filename::String)::Coefficients

Read the Fourier coefficients of a path from a file.
"""
function read_path_from_file(filename::String, F::Int, N::Int, dim::Int)::Coefficients
    Γ = FromZero([[zeros(Float64, dim) for i ∈ 1:N] for j ∈ 0:F+1])
    data = readdlm(filename, ' ')

    for k in axes(data, 1)
        h = (k - 1) % (F + 2)
        i = div((k - 1), F + 2) + 1
        Γ[h][i] = data[k, :]
    end
    return Γ
end

function path_animation()
    error("GLMakie must be loaded")
end
