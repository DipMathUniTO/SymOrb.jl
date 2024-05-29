
fourier_series(A::Coefficients, n=steps+1)::Path  = [sum(A[k] * sin(π * k * h  / n ) for k ∈ 1:F) for h ∈ 0:n]

inverse_fourier(x::Path, n=steps+1)::Coefficients =  2 / n * [sum(x[h] * sin(k * h * π / n) for h ∈ 1:n) for k ∈ 0:F+1]

segment(a, b, n=steps+1) = [a + k * (b - a) / n for k ∈ 0:n]


function build_path(A::Coefficients, n::Int = steps+1)::Path
    y = segment(A[0], A[F+1], n) + fourier_series(A, n)
    return FromZero(y)
end


function fourier_coefficients(y::Path, n::Int = steps+1)::Coefficients
    x = y - FromZero(segment(y[0], y[n]))
    A = FromZero(inverse_fourier(x))
    A[0] = y[0]
    A[F+1] = y[n]
    return FromZero(A)
end


function emboss(v::Vector{T})::Coefficients where {T}
    Γ = FromZero([[zeros(T, dim) for i ∈ 1:N] for k ∈ 0:F+1])

    for j ∈ 1:dim, i ∈ 1:N, k ∈ 0:F+1
        Γ[k][i][j] = v[k*N*dim+(i-1)*dim+j]
    end

    return Γ
end


function flatten(Γ::AbstractVector{Vector{Vector{T}}})::Vector{T} where {T}
    v = zeros(T, (F+2) * N * dim)

    for k ∈ 0:F+1, i ∈ 1:N, j ∈ 1:dim
        v[k*N*dim+(i-1)*dim+j] = Γ[k][i][j]
    end

    return v
end


function emboss(M::Matrix{T})::OffsetArray where {T}
    H =  FromZero([[zeros(T, dim, dim) for _ ∈ 1:N, _ ∈ 1:N] for _ ∈ 0:F+1, _ ∈ 0:F+1 ])

    for k1 ∈ 0:F+1, k2 ∈ 0:F+1, i1 ∈ 1:N, i2 ∈ 1:N, j1 ∈ 1:dim, j2 ∈ 1:dim
        H[k1, k2][i1, i2][j1, j2] = M[k1*N*dim+(i1-1)*dim+j1, k2*N*dim+(i2-1)*dim+j2]
    end

    return H
end


function flatten(H::AbstractMatrix{Matrix{Matrix{T}}})::Matrix{T} where {T}
    M = zeros(T, (F+2) * N * dim, (F+2) * N * dim)

    for k1 ∈ 0:F+1, k2 ∈ 0:F+1, i1 ∈ 1:N, i2 ∈ 1:N, j1 ∈ 1:dim, j2 ∈ 1:dim
        M[k1*N*dim+(i1-1)*dim+j1, k2*N*dim+(i2-1)*dim+j2] = H[k1, k2][i1, i2][j1, j2]
    end

    return M
end

function get_starting_path(path_type::Symbol)::OffsetArray
    if path_type == :random
        return random_starting_path()
    elseif path_type == :circular
        return circular_starting_path()
    elseif path_type == :perturbed_circular
        return perturbed_circular_starting_path()
    else
        return random_starting_path()
    end
end

random_starting_path()::OffsetArray =  FromZero([[rand(Float64, dim) .- 0.5 for i ∈ 1:N] for j ∈ 0:F+1])

circular_starting_path()::Coefficients = begin
    x::Path = FromZero([[zeros(dim) for i ∈ 1:N] for j ∈ 1:steps+2])

    for h ∈ 0:steps+1, i ∈ 1:N
        x[h][i][1] = cos(h * dt + (i-1) * 2 * π / N)
        x[h][i][2] = sin(h * dt + (i-1) * 2 * π / N)
    end

    return (FromZero ∘ fourier_coefficients)(x)
end

perturbed_path(x::Path, λ = 0.001 )::Coefficients = x + λ * random_starting_path()

perturbed_circular_starting_path(λ::Float64 = 0.001)::Coefficients = circular_starting_path() + λ * random_starting_path()

function reconstruct_path(x::Path)
    n = lastindex(x)-1
    initial_path = FromZero([[zeros(dim) for i ∈ 1:N] for j ∈ 0:(2*n+1)])
    complete_path = FromZero([])
    initial_path[0:n+1] = copy(x[0:n+1])

    if action_type == Cyclic 
        max_index = n + 1
    else
        max_index = 2*n + 1 
        for k ∈ 0:n
            initial_path[n + k + 1] = ϕ(H_1, initial_path[n - k + 1])     
        end
    end
    
    for j ∈ 1:cyclic_order, k ∈ 0:max_index
        push!(complete_path, ϕg_n(initial_path[k], j))
    end 

    push!(complete_path, ϕg_n(initial_path[0], 1))

    return complete_path
end



function print_path_to_file(Γ, filename)
    open(filename, "w") do file
        for i in 1:N, h in 0:F+1
            write(file, join(string.(Γ[h][i]), " "), "\n")
        end
    end
end

function read_path_from_file(filename, config)
    LSG_from_config(config)
    Γ = FromZero([[zeros(Float64, dim) for i ∈ 1:N] for j ∈ 0:F+1])
    data = []
    open(filename, "r") do file
        data = readlines(file)
    end 

    for i in 1:length(data)
        t = (i-1) % (F+2)
        n = div((i-1), F+2) + 1
        Γ[t][n] = parse.(Float64, split(data[i]))
    end
    return Γ

end

 refine_path(Γ::Coefficients, n::Int = steps+1)::Path = reconstruct_path(build_path(Γ, n))
