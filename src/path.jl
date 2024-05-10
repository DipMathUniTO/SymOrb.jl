
fourier_series(A::Coefficients)::Path  = [sum(A[k] * sin(k * h * dt) for k ∈ 1:F) for h ∈ 0:steps+1]

inverse_fourier(x::Path)::Coefficients = 1/(steps+1)*[sum(x[h] * sin(k * h * dt) for h ∈ 1:steps) for k ∈ 0:F+1]

segment(a, b) = [a + k * (b - a) / (steps + 1) for k ∈ 0:steps+1]


function build_path(A::Coefficients)::Path
    y = segment(A[0], A[F+1]) + fourier_series(A)
    return OffsetArrays.Origin(0)(y)
end


function emboss(v::Vector{T})::Coefficients where {T}
    Γ = OffsetArray([[zeros(T, dim) for i ∈ 1:N] for k ∈ 0:F+1], 0:F+1)

    for j ∈ 1:dim, i ∈ 1:N, k ∈ 0:F+1
        Γ[k][i][j] = v[k*N*dim+(i-1)*dim+j]
    end

    return Γ
end

function flatten(Γ::OffsetVector{Vector{Vector{T}}}) where {T}
    return [((Γ...)...)...]
end



function random_starting_path()::OffsetArray
    A = [[rand(Float64, dim) .- 0.5 for i ∈ 1:N] for j ∈ 0:F+1]
    return OffsetArrays.Origin(0)(A)
end


function circular_starting_path()::OffsetArray
    A::Path = OffsetArrays.Origin(0)([[zeros(dim) for i ∈ 1:N] for j ∈ 1:steps+2])

    for h ∈ 0:steps+1
        for i ∈ 1:N
            A[h][i][1] = cos(h * dt + (i-1) * 2 * π / N)
            A[h][i][2] = sin(h * dt + (i-1) * 2 * π / N)
        end
    end
    return OffsetArrays.Origin(0)(inverse_fourier(A))
end