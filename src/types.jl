# ==================== CONFIGURATIONS AND PATHS ====================   

# Utility renaming for the types used in the MinPathJL package
const Config = Vector{Vector{Float64}}
const Path = OffsetVector{Config}
const Coefficients = OffsetVector{Config}

const Permutation = Vector{Int}
const Rotation = Matrix{Float64}

const FromZero = Origin(0)

# Define the products used to compute the kinetic energy
"""
    Base.:*(M::OffsetMatrix{Matrix{Matrix{T}}}, v::Coefficients)::Coefficients where {T}

Multiplication of a matrix of matrices by Coefficients. It's used to multiply the kinetic energy matrix by the Fourier coefficients.
"""
function ⊙(M::OffsetMatrix{Matrix{Matrix{T}}}, v::Coefficients)::Coefficients where {T}
    result = [[zeros(T, length(v[j][1])) for _ ∈ 1:length(v[j])] for j ∈ axes(v, 1)]
    for h ∈ axes(M, 1), k ∈ axes(M, 2), i ∈ axes(M[1,1], 1),  j ∈ axes(M[1,1], 2)
        result[h][j] += M[h, k][i, j] * v[k][i]
    end
    return result
end

"""
    Base.:*(v::OffsetMatrix{Config}, w::Coefficients)::Config

Multiplication of a matrix of configurations by Coefficients. It's used to multiply the transposed coefficients by the coefficients
"""
function ⊗(v::Coefficients, w::Coefficients)::Float64
    result = 0

    for j ∈ axes(v, 1), k ∈ axes(w, 1), i ∈ axes(v[end], 1)
        result += sum(v[j][i] .*  w[k][i])
    end

    return result
end

"""
    Base.:*(M::Matrix, v::Config)::Config

Multiplication of a matrix by a configuration. It's used to multiply the rotation matrix by the configuration.
"""
function ⊙(M::Matrix, v::Config)
    return [[M * v[i][j] for j ∈ axes(v[i], 1)] for i ∈ axes(v, 1)]
end

# ==================== SYMMETRY GROUP ====================

""" The type of the action can by Cyclic, Dihedral or Brake """
@enum ActionType begin
    Cyclic = 0
    Dihedral = 1
    Brake = 2
end

"""
A group element of the symmetry group

# Fields
- `σ::Permutation`: The permutation of the particles
- `M::Rotation`: The rotation of the particles
"""
struct GroupElement
    σ::Permutation
    M::Rotation
end
GroupElement() = GroupElement(Permutation([]), Rotation(undef, 0, 0))

"""
The symmetry group of the minimization problem

# Fields
- `action_type::ActionType`: The type of the action
- `kerT::Vector{GroupElement}`: The kernel of `τ`, the action on `O(2)`
- `g::Vector{GroupElement}`: The elements of the cyclic part of the symmetry group
- `H0::GroupElement`: If the action is dihedral or brake, it's the generator of the first isotropy subspaces
- `H1::GroupElement`: If the action is dihedral or brake, it's the generator of the second isotropy subspaces
"""
struct SymmetryGroup
    action_type::ActionType
    kerT::Vector{GroupElement}
    g::Vector{GroupElement}
    H0::GroupElement
    H1::GroupElement
end

cyclic_order(G::SymmetryGroup) = length(G.g)


# ==================== PROBLEM ====================


const DimensionsTuple = @NamedTuple{F::Int, N::Int, dim::Int}

"""
The minimization problem to be solved

# Fields
- `N::Int`: The number of particles
- `dim::Int`: The dimension of the space
- `F::Int`: The number of Fourier series terms
- `steps::Int`: The number of steps in the discretization of time [0,π]
- `G::SymmetryGroup`: The symmetry group of the minimization problem
- `m::Vector{Float64}`: The masses of the particles
- `f::A`: The denominator of the potential
- `K::OffsetMatrix{Matrix{Matrix{Float64}}}`: The kinetic energy matrix
- `dx_dAk::OffsetVector{Vector{Float64}}`: The derivative of the path w.r.t. the Fourier coefficients
"""
struct Problem{M<:Function}
    G::SymmetryGroup
    m::Vector{Float64}
    f::M
    K::OffsetMatrix{Matrix{Matrix{Float64}}}
    dx_dAk::OffsetVector{Vector{Float64}}
end


function Base.show(io::IO, P::Problem)
    println(io, "\nSymmetry group: \t", P.G)
    println(io, "Denominator: \t", P.f)
end

"""
    atomic_method(type, struct_name, f_name)

Define a new atomic minimization method with the given `type`, `struct_name` and `f_name`.  
"""
macro atomic_method(type, struct_name, f_name)
    quote
        @kwdef struct $struct_name <: $type
            f_name::Symbol = $f_name
            max_iter::Int = 200
            show_trace::Bool = false
            tolerance::Float64 = 1e-8
        end

        # Define a constructor that accepts a single Int for max_iter
        function $(esc(struct_name))(max_iter::Int)
            $(esc(struct_name))(max_iter=max_iter)
        end

    end
end

# ====================== MINIMIZATION METHODS ======================


# Define the abstract types for the minimization methods
abstract type AbstractMethod end
abstract type AbstractAtomicMethod <: AbstractMethod end
abstract type AbstractCompoundMethod <: AbstractMethod end
# Define an alias for the atomic methods or nothing
const AbstrAtomicOrNothing = Union{AbstractAtomicMethod, Nothing}

# Define the abstract atomic minimization methods
abstract type OptimMethod <: AbstractAtomicMethod end
abstract type NLSolveMethod <: AbstractAtomicMethod end

# Define the concrete atomic minimization methods
@atomic_method NLSolveMethod Newton :trust_region
@atomic_method OptimMethod BFGS :BFGS
@atomic_method OptimMethod ConjugateGradient :ConjugateGradient
@atomic_method OptimMethod GradientDescent :GradientDescent

# Define the compound minimization methods
struct CompoundMethod{M <: AbstrAtomicOrNothing, N <: AbstrAtomicOrNothing, O <: AbstrAtomicOrNothing} <: AbstractCompoundMethod
    init::M
    first::N
    second::O
end
# Define the aliases for the compound methods
const OneMethod{N} = CompoundMethod{Nothing, N, Nothing}
const InitMinimizeMethod{M, N} = CompoundMethod{M, N, Nothing}

OneMethod(method) = CompoundMethod(nothing, method, nothing)
InitMinimizeMethod(init, method) = CompoundMethod(init, method, nothing)

"""
    Base.show(io::IO, method::AbstractMethod)

Pretty-print the minimization `method` in a concise way
"""
function Base.show(io::IO, method::AbstractMethod)
    print(io, typeof(method), "(max_iter=", method.max_iter, ", tolerance=", method.tolerance, ")")
end

"""
    Base.show(io::IO, method::AbstractAtomicMethod)

Pretty-print the minimization `method` in a more detailed way
"""
function Base.show(io::IO, ::MIME"text/plain", method::AbstractMethod)
    println(io, typeof(method), " minimization method:")
    println(io, "\tIterations:\t", method.max_iter)
    println(io, "\tShow trace:\t", method.show_trace)
    println(io, "\tTolerance:\t", method.tolerance)
end

"""
    Base.show(io::IO, method::AbstractCompoundMethod)

Specialized pretty-printing for the compound minimization `method` in a concise way
"""
function Base.show(io::IO, method::AbstractCompoundMethod)
    print(io, typeof(method))
end

function Base.show(io::IO, ::MIME"text/plain", method::AbstractCompoundMethod)
    println(io, typeof(method))
    println(io, "  Init:\t\t", typeof(method.init))
    println(io, "  First:\t", typeof(method.first))
    println(io, "  Second:\t", typeof(method.second))
end

# ================== MINIMIZATION RESULT ==================

"""
    MinimizationResult{T <: AbstractMethod}(Γ::Coefficients, Γ_min::Coefficients, method, converged)

The result of a minimization

# Fields
- `initial::Coefficients`: The initial Fourier coefficients
- `fourier_coeff::Coefficients`: The Fourier coefficients of the minimum
- `path`: The path corresponding to the Fourier coefficients of the minimum
- `iterations::Int`: The number of iterations
- `gradient_norm::Float64`: The norm of the gradient at the minimum
- `action_value::Float64`: The value of the action at the minimum
- `converged::Bool`: Whether the minimization converged
- `method::T`: The minimization method
"""
@kwdef struct MinimizationResult{T <: AbstractMethod}
    initial::Coefficients
    fourier_coeff::Coefficients
    path::Path
    iterations::Union{Int, Nothing}
    gradient_norm::Float64
    action_value::Float64
    converged::Bool
    method::T
end

""" 
    MinimizationResult(Γ::Coefficients, Γ_min::Coefficients, method, converged)

Create a MinimizationResult from the initial and final Fourier coefficients, the minimization method and whether it converged.
"""
function MinimizationResult(Γ::Coefficients, Γ_min::Coefficients, method, converged)
    MinimizationResult(
        initial = Γ,
        fourier_coeff = Γ_min,
        path = (extend_to_period ∘ build_path)(Γ_min),
        iterations = nothing,
        gradient_norm = norm(∇action(Γ_min)),
        action_value = action(Γ_min),
        converged = converged,
        method = method
    )
end

"""
    MinimizationResult(Γ::Coefficients, minimizer::Vector{Float64}, method, converged)

Create a MinimizationResult from the initial Fourier coefficients, the final coefficients in a flattened form, 
the minimization method and whether it converged.
"""
function MinimizationResult(Γ::Coefficients, minimizer::Vector{Float64}, method, converged)
    Γ_min = (project ∘ emboss)(minimizer)
    MinimizationResult(
        initial = Γ,
        fourier_coeff = Γ_min,
        path = (extend_to_period ∘ build_path)(Γ_min),
        iterations = NaN,
        gradient_norm = norm(∇action(Γ_min)),
        action_value = action(Γ_min),
        method = method,
        converged = converged
    )end

"""
    MinimizationResult(Γ::Coefficients, minimizer::Vector{Float64}, method, converged)

Pretty-print the minimization result
"""
function Base.show(io::IO, result::MinimizationResult)
    println(io, "Minimization result:")
    println(io, "\tMethod: \t", result.method)
    print(io, "\tConverged: \t")
    printstyled(io, result.converged, color=result.converged ? :green : :red)
    println(io)
    println(io, "\tGradient norm: \t", result.gradient_norm)
    println(io, "\tAction value: \t", result.action_value)
end

Base.zeros(T::Coefficients, size::Int, N::Int, dim::Int) = OffsetArray([[zeros(dim) for _ ∈ 1:N] for _ ∈ 0:size], 0:size)


function dims(Γ::AbstractArray{Array{Array{Float64, T}, T}, T}) where T
    F = size(Γ, 1)-2
    N = size(Γ[end], 1)
    dim = size(Γ[end][1], 1)
    return F, N, dim
end

