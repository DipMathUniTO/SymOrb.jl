
fourier_series(A::Coefficients)::Path = [sum(A[k] * sin(π * k * h / (steps + 1)) for k in 1:F) for h in 0:steps+1]

inverse_fourier(A::Path)::Coefficients = [sum(A[h] * sin(π * k * h / (steps + 1)) for h in 1:steps) for k in 0:F+1]
segment(a, b) = [a + k * (b - a) / (steps + 1) for k ∈ 0:steps+1]



function build_path(A::Coefficients)::Path
    y = segment(A[0], A[F+1]) + fourier_series(A)
    return OffsetArrays.Origin(0)(y)
end


function emboss(v::Vector{T})::Coefficients where {T}
    Γ = OffsetArray([[zeros(T, dim) for i ∈ 1:N] for k ∈ 0:F+1], 0:F+1)

    for j in 1:dim, i in 1:N, k in 0:F+1
        Γ[k][i][j] = v[k*N*dim+(i-1)*dim+j]
    end

    return Γ
end

function flatten(Γ::OffsetVector{Vector{Vector{T}}}) where {T}
    return [((Γ...)...)...]
end



function random_starting_path()::OffsetArray
    A = [[rand(Float64, dim) .- 0.5 for i in 1:N] for j in 0:F+1]
    return OffsetArrays.Origin(0)(A)
end


function circular_starting_path()::OffsetArray
    A::Path = OffsetArrays.Origin(0)([[zeros(dim) for i in 1:N] for j in 1:steps+2])

    for h in 0:steps+1
        for i in 1:N
            A[h][i][1] = cos(π * h / (steps + 1) + (i - 1) * 2 * π / N)
            A[h][i][2] = sin(π * h / (steps + 1) + (i - 1) * 2 * π / N)
        end
    end
    return OffsetArrays.Origin(0)(inverse_fourier(A))
end



function path_from_minimizer(v::Vector{T})::OffsetArray where {T}
    path = (build_path ∘ project ∘ emboss)(v)
    return [ [[path[i][j][k] for k in 1:dim] for i in 1:steps ] for j in 1:N ]
end


function coeff_from_minimizer(v::Vector{T})::OffsetArray where {T}
    path = (project ∘ emboss)(v)
    return [ [[path[i][j][k] for k in 1:dim] for i in 1:F ] for j in 1:N ]
end