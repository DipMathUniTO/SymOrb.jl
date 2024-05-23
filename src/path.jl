
fourier_series(A::Coefficients)::Path  = [sum(A[k] * sin(k * h * dt) for k ∈ 1:F) for h ∈ 0:steps+1]

inverse_fourier(x::Path)::Coefficients =  2/(steps+1)*[sum(x[h] * sin(k * h * dt) for h ∈ 1:steps) for k ∈ 0:F+1]

segment(a, b) = [a + k * (b - a) / (steps + 1) for k ∈ 0:steps+1]


function build_path(A::Coefficients)::Path
    y = segment(A[0], A[F+1]) + fourier_series(A)
    return OffsetArrays.Origin(0)(y)
end


function fourier_coefficients(y::Path)::Coefficients
    x = y - OffsetArrays.Origin(0)(segment(y[0], y[steps+1]))
    A = inverse_fourier(x)
    A[1] = y[0]
    A[F+2] = y[steps+1]
    return OffsetArrays.Origin(0)(A)
end


function emboss(v::Vector{T})::Coefficients where {T}
    Γ = OffsetArray([[zeros(T, dim) for i ∈ 1:N] for k ∈ 0:F+1], 0:F+1)

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
    H =  OffsetArray([[zeros(T, dim, dim) for _ ∈ 1:N, _ ∈ 1:N] for _ ∈ 0:F+1, _ ∈ 0:F+1 ], 0:F+1, 0:F+1)

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


function random_starting_path()::OffsetArray
    A = [[rand(Float64, dim) .- 0.5 for i ∈ 1:N] for j ∈ 0:F+1]
    return OffsetArrays.Origin(0)(A)
end


function circular_starting_path()::OffsetArray
    x::Path = OffsetArrays.Origin(0)([[zeros(dim) for i ∈ 1:N] for j ∈ 1:steps+2])

    for h ∈ 0:steps+1
        for i ∈ 1:N
            x[h][i][1] = cos(h * dt + (i-1) * 2 * π / N)
            x[h][i][2] = sin(h * dt + (i-1) * 2 * π / N)
        end
    end

    A = fourier_coefficients(x)
    plot_path(x)
    return OffsetArrays.Origin(0)(A)
end



function reconstruct_path(x::Path)
    cyclic_x = OffsetArray([[zeros(dim) for i ∈ 1:N] for j ∈ 0:(2*steps+1)], 0:(2*steps+1))
    
    cyclic_x[0:steps+1] = copy(x[0:steps+1])

    max_index = if action_type == Cyclic 
        steps+1
    else 
        2*steps +1
    end


    if action_type != Cyclic
        for k ∈ 0:steps
            cyclic_x[steps + k + 1] = ϕ(H_1, cyclic_x[steps - k + 1])     
        end
    end
    
    #plot_path(cyclic_x)
    #return cyclic_x
    

    reconstructed_x = OffsetArrays.Origin(0)([])

    for j ∈ 1:cyclic_order, k ∈ 0:max_index
        push!(reconstructed_x, ϕg_n(cyclic_x[k], j))
    end 
    push!(reconstructed_x, ϕg_n(cyclic_x[0], 1))

    plot_path(reconstructed_x)
    return reconstructed_x
end


function plot_path(path)
    pl = plot(aspect_ratio=:equal)
    l = length(path)-1
    for i ∈ 1:N
        pl = plot!(pl, [[path[j][i][h] for j ∈ 0:l] for h ∈ 1:dim]..., label="Body $i", linealpha=0.5, linewidth=3,  aspect_ratio=:equal);
    end
    # anim = [Animation() for i ∈ 1:N]
    # for i ∈ 1:N
    #     anim[i] = @animate for j ∈ 1:l
    #         plot([[path[j][i][h] for h ∈ 1:dim]...], label="Body $i", linealpha=0.5, linewidth=3, aspect_ratio=:equal)
    #     end
    # end

    display(pl)

end

function print_path_to_file(Γ, filename)
    open(filename, "w") do file
        for i in 1:N
              for h in 0:F
                for d  in 1:dim
                    write(file, string(Γ[h][i][d]) * " ")
                end
                write(file, "\n")
              end
        end
        end
end
