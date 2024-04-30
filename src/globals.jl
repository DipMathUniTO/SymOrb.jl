

Config = Vector{Vector{Float64}}

Path = OffsetVector{Config}
Coefficients = OffsetVector{Config}

Permutation = Vector{Int64}
Rotation = Matrix


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


Base.@kwdef struct SymorbConfig
    N::Int64
    dim::Int64
    action_type::Int64
    kerT::String
    rotV::Matrix
    rotS::String
    refV::Matrix
    refS::String
    Omega::Matrix
end



N::Int64 = 0
dim::Int64 = 0
action_type::ActionType = Cyclic
cyclic_order::Int64 = 0
dim_kerT::Int64 = 0
kerT::Union{G, Vector{G}} = G()
isΩ::Bool = false
g::Vector{G} = []
H_0::G = G()
H_1::G = G()
m::Vector{Float64} = []

Ω::Matrix{Float64} = Matrix(undef, 0, 0)

K::OffsetMatrix{Float64} = Matrix(undef, 0, 0)



function Base.:*(M::OffsetMatrix{T}, v::Coefficients)::Coefficients where {T}
    result = [[zeros(T, length(v[j][1])) for _ in 1:length(v[j])] for j in axes(v,1)]

    for j ∈ axes(M, 1), k ∈ axes(M, 2), i ∈ 1:length(v[j])
        result[j][i] += m[i] * M[j, k]  * v[k][i]
    end
    return result
end

function Base.:*(v::OffsetMatrix{Config}, w::Coefficients):Coefficients
    result = zero(T)

    for j ∈ axes(v, 2), k ∈ axes(w, 1),i ∈ 1:length(v[j])
        result += v[j][i] * w[k][i]
    end

    return result
end


function spatial_mult(M::Matrix{T}, v::Config) where {T}
    return [[M * v[i][j] for j in axes(v[i],1)] for i in axes(v,1)]
end
