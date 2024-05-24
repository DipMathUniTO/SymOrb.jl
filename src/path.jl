
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


function random_starting_path()::OffsetArray
    A = [[rand(Float64, dim) .- 0.5 for i ∈ 1:N] for j ∈ 0:F+1]
    return OffsetArrays.Origin(0)(A)
end


function circular_starting_path()::OffsetArray
    x::Path = FromZero([[zeros(dim) for i ∈ 1:N] for j ∈ 1:steps+2])

    for h ∈ 0:steps+1
        for i ∈ 1:N
            x[h][i][1] = cos(h * dt + (i-1) * 2 * π / N)
            x[h][i][2] = sin(h * dt + (i-1) * 2 * π / N)
        end
    end

    A = fourier_coefficients(x)
    plot_path(x)
    return FromZero(A)
end

function perturbed_circular_starting_path()::OffsetArray
    x::Path = FromZero([[zeros(dim) for i ∈ 1:N] for j ∈ 1:steps+2])

    for h ∈ 0:steps+1
        for i ∈ 1:N
            x[h][i][1] = cos(h * dt + (i-1) * 2 * π / N) + 0.1 * randn()
            x[h][i][2] = sin(h * dt + (i-1) * 2 * π / N) + 0.1 * randn()
        end
    end

    A = fourier_coefficients(x)
    plot_path(x)
    return FromZero(A)
end


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


function plot_path(path)
    theme = theme_dark()
    set_theme!(theme)
    f = Figure(size=(800, 800))
    ax = nothing
    if dim == 2
        ax  = Axis(f[1,1],autolimitaspect = 1)
    elseif dim == 3
        ax  = Axis3(f[1,1],aspect = :data)
    else 
        return
    end

    hidedecorations!(ax)
    l = length(path)-1
    @show l
    for i ∈ 1:N
        lines!(ax, [[path[j][i][h] for j ∈ 0:l] for h ∈ 1:dim]..., color=i, colormap=:lightrainbow, colorrange = (1,N), label="Body $i");
    end

    display(f)
    return ax
end


function path_animation(path, fps::Int = 30)
    
    ax = plot_path(path)

    if dim == 2
        bodies = Observable([Point2(zeros(2)) for _ ∈ 1:N])
    elseif dim == 3
        bodies = Observable([Point3(zeros(3)) for _ ∈ 1:N])
    else 
        return
    end
    scatter!(ax, bodies, markersize=40, color=1:N, colormap=:lightrainbow, colorrange=(1,N))
    t = 0
    while true
        if(t == length(path)) t = 0 end
        bodies[] = path[t]
        t += 1
        sleep(1.0/fps)
    end
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

    for i in axis(data)
        t = (i-1) % (F+2)
        n = div((i-1), F+2) + 1
        @show t, n
        Γ[t][n] = parse.(Float64, split(data[i]))
    end
    return Γ

end

