

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

fourier_series(A::Coefficients)::Path = [sum(A[k]*sin(π*k*h/(steps+1)) for k in 1:F) for h in 0:steps+1]

inverse_fourier(A::Path)::Coefficients = [sum(A[h]*sin(π*k*h/(F+1)) for h in 1:steps    ) for k in 0:F+1]

segment(a,b) = [a + k * (b - a) / (steps+1)  for k ∈ 0:steps+1]


function build_path(A::Coefficients)::Path
    y = segment(A[0], A[F+1]) + fourier_series(A)
    return OffsetArrays.Origin(0)(y)
end


function emboss(v)
    Γ = Vector{Vector{Vector{Float64}}}(undef, F+2)
    for k in 1:F+2
        Γ[k] = [zeros(dim) for _ in 1:N]
        for i in 1:N
            for j in 1:dim
                Γ[k][i][j] = v[(k-1)*N*dim + (i-1)*dim + j]
            end
        end
    end

    return OffsetArrays.Origin(0)(Γ)
end

function flatten(Γ)
    return [((Γ...)...)...]
end


