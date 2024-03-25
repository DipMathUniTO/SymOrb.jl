

Config = Vector{Vector{Float64}}

Path = OffsetArray{Config}
Coefficients = OffsetArray{Config}

Permutation = Vector{Int64}
Rotation = Matrix{Float64}


struct G
    σ::Permutation
    M::Rotation
end

@enum ActionType begin
    Cyclic = 0
    Dihedral = 1
    Brake = 2
end

G() = G(Permutation([]), Rotation(undef, 0, 0))


N::Int64 = 0
dim::Int64 = 0
action_type::ActionType = Cyclic
cyclic_order::Int64 = 0
dim_kerT::Int64 = 0
kerT = G()
isΩ::Bool = false
g = []
H_0 = G()
H_1 = G()
m = []

Ω::Matrix{Float64} = Matrix(undef, 0, 0)
K::OffsetMatrix{Float64} = Matrix(undef, 0, 0)


fourier_series(A::Coefficients)::Path = [sum(A[k]*sin(π*k*h/(steps+1)) for k in 1:F) for h in 0:steps+1]

inverse_fourier(A::Path)::Coefficients = [sum(A[h] * sin(π * k * h / (steps + 1)) for h in 1:steps) for k in 0:F+1]
segment(a,b) = [a + k * (b - a) / (steps+1)  for k ∈ 0:steps+1]



function build_path(A::Coefficients)::Path
    y = segment(A[0], A[F+1]) + fourier_series(A)
    return OffsetArrays.Origin(0)(y)
end


function emboss(v::Vector{T}) where T
    Γ = OffsetArray([[zeros(T, dim) for i ∈ 1:N] for k ∈0:F+1], 0:F+1)

    for j in 1:dim,  i in 1:N, k in 0:F+1
        Γ[k][i][j] = v[k*N*dim + (i-1)*dim + j]
    end

    return Γ
end

function flatten(Γ::OffsetVector{Vector{Vector{T}}}) where {T}
    return [((Γ...)...)...]
end


function Base.:*(M::OffsetMatrix{T}, v::OffsetVector{Vector{Vector{T}}}) where {T}
    result = OffsetVector{Vector{Vector{T}}}(undef, axes(v,1))
    result = [[zeros(T, length(v[j][1])) for _ in 1:length(v[j])] for j in axes(v,1)]

    for j ∈ axes(M, 1), k ∈ axes(M, 2), i ∈ 1:length(v[j])
        result[j][i] += m[i] * M[j, k]  * v[k][i]
    end
    return result
end

function Base.:*(v::OffsetMatrix{Vector{Vector{T}}}, w::OffsetVector{Vector{Vector{T}}}) where {T}
    result = zero(T)

    for j ∈ axes(v, 2), k ∈ axes(w, 1),i ∈ 1:length(v[j])
        result += v[j][i] * w[k][i]
    end

    return result
end


function spatial_mult(M::Matrix{T}, v::OffsetVector{Vector{Vector{T}}}) where {T}
    return [[M * v[i][j] for j in axes(v[i],1)] for i in axes(v,1)]
end
