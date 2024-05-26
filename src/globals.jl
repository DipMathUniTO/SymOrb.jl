using Base: @kwdef

Config = Vector{Vector{Float64}}

Path = OffsetVector{Config}
Coefficients = OffsetVector{Config}

Permutation = Vector{Int64}
Rotation = Matrix

FromZero = OffsetArrays.Origin(0) 

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

global F::Int64 = 0
global N::Int64 = 0
global dim::Int64 = 0
global action_type::ActionType = Cyclic
global cyclic_order::Int64 = 0
global dim_kerT::Int64 = 0
global kerT::Union{G, Vector{G}} = G()
global isΩ::Bool = false
global dt::Float64 = 0
global g::Vector{G} = []
global H_0::G = G()
global H_1::G = G()
global m::Vector{Float64} = []
global ϵ::Float64 = 0
global Ω::Matrix{Float64} = Matrix(undef, 0, 0)

global K::OffsetMatrix = Matrix(undef, 0, 0)
global Id = Matrix(undef, 0,0)

global dx_dAk::OffsetVector = Vector(undef, 0)

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

abstract type AbstractMethod end
abstract type OptimMethod <: AbstractMethod end
abstract type NLSolveMethod <: AbstractMethod end

macro new_method(struct_name, f_name, type)
    quote
        @kwdef struct $struct_name <: $type
            f_name::Symbol = $f_name
            max_iter::Int64 = 200
            show_trace::Bool = false
            tolerance::Float64 = 1e-8
            iterations::Int64 = 0
        end
        
        # Define a constructor that accepts a single Int for max_iter
        function $(esc(struct_name))(max_iter::Int)
            $(esc(struct_name))(max_iter = max_iter)
        end
        
    end
end

@new_method Newton :trust_region NLSolveMethod
@new_method BFGS :BFGS OptimMethod
@new_method ConjugateGradient :ConjugateGradient OptimMethod

@kwdef struct Methods
    init::Union{AbstractMethod, Nothing} = nothing
    first::Union{AbstractMethod, Nothing} = BFGS()
    second::Union{AbstractMethod, Nothing} = nothing
end

Methods(one_method::AbstractMethod) = Methods(init = nothing, first = one_method, second = nothing)
Methods(init_method::AbstractMethod, method::AbstractMethod) = Methods(init = init_method, first = method)



@kwdef struct MinimizationResult
    initial::Coefficients
    fourier_coeff::Coefficients
    path::Path
    iterations::Int64
    gradient_norm::Float64
    action_value::Float64
    converged::Bool
    method::AbstractMethod
end

MinimizationResult(Γ, minimizer, method, iterations, converged) = begin
   Γ_min = (project ∘ emboss)(minimizer) 
   MinimizationResult(
    initial = Γ,
    fourier_coeff = Γ_min,
    path = (reconstruct_path ∘ build_path)(Γ_min),
    iterations = iterations,
    gradient_norm = norm(∇action(minimizer)),
    action_value = action(minimizer),
    converged = converged,
    method=method
)
end 



Base.show(io::IO, method::AbstractMethod) = begin
    print(io, typeof(method), "(max_iter=", method.max_iter, ", tolerance=", method.tolerance, ")")
end

Base.show(io::IO, ::MIME"text/plain", method::AbstractMethod) = begin
    println(io, typeof(method), " minimization method:")
    println(io, "\tIterations:\t", method.max_iter)
    println(io, "\tShow trace:\t", method.show_trace)
    println(io, "\tTolerance:\t", method.tolerance)
end

Base.show(io::IO, methods::Methods) = begin
    println(io, "Methods:")
    println(io, "  Init:\t\t", typeof(methods.init))
    println(io, "  First:\t", typeof(methods.first))
    println(io, "  Second:\t", typeof(methods.second))
end

Base.show(io::IO, result::MinimizationResult) = begin
    println("Minimization result:")
    println(io, "\tMethod: \t", result.method)
    print(io, "\tConverged: \t")
    printstyled(result.converged, color = result.converged ? :green : :red)
    println()
    println(io, "\tIterations: \t", result.iterations)
    println(io, "\tGradient norm: \t", result.gradient_norm)
    println(io, "\tAction value: \t", result.action_value)
end
