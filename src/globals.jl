using Base: @kwdef

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



@kwdef struct MinimizationOptions
    g_tol::Float64 = 1e-8
    show_trace::Bool = false
    extended_trace::Bool = true
    callback::Function = x -> false
end

MinimizationOptions(options::Optim.Options) = MinimizationOptions(
    g_tol = options.g_abstol, 
    show_trace = options.show_trace,
    extended_trace = options.extended_trace,
    callback = options.callback)

G() = G(Permutation([]), Rotation(undef, 0, 0))


@kwdef struct MinimizationResult
   fourier_coeff::Coefficients
   trajectory::Path
   iterations::Int64
   gradient_norm::Float64
   action_value::Float64
   converged::Bool
   initial::Coefficients
   minimization_method::Symbol
   minimization_options::MinimizationOptions
end


@kwdef struct SymorbConfig
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
dt::Float64 = 0
g::Vector{G} = []
H_0::G = G()
H_1::G = G()
m::Vector{Float64} = []

Ω::Matrix{Float64} = Matrix(undef, 0, 0)

K::OffsetMatrix = Matrix(undef, 0, 0)
Id = Matrix(undef, 0,0)

dx_dAk::OffsetVector = Vector(undef, 0)


function Base.:*(M::OffsetMatrix{Matrix{Matrix{T}}}, v::Coefficients)::Coefficients where {T}
    result = [[zeros(T, length(v[j][1])) for _ ∈ 1:length(v[j])] for j ∈ axes(v,1)]

    for h ∈ axes(M, 1), k ∈ axes(M, 2), i ∈ 1:length(v[1]), j ∈ 1:length(v[1])
        result[h][j] +=  M[h, k][i, j]  * v[k][i]
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

function Base.:*(M::Matrix, v::Config)
    return [[M * v[i][j] for j ∈ axes(v[i],1)] for i ∈ axes(v,1)]
end