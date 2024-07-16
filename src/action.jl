"""
action(Γ::Coefficients)::Float64

Compute the constrained action for a given configuration ``Γ``.
    """
function action(p::Problem, Γ::Coefficients)::Float64
    _project = x -> project(p, x)
    _action = x ->  kinetic(p, x) +  potential(p, x)
    (_action ∘ _project)(Γ)
end


""" 
    ∇action(Γ::Coefficients)::AbstractVector

Compute the gradient of the constrained action for a given configuration ``Γ``.
"""
function ∇action(p::Problem, Γ::Coefficients)::Coefficients
    _project = x -> project(p, x) 
    _∇action = x -> ∇kinetic(p, x) + ∇potential(p, x)
    (_project ∘ _∇action ∘ _project)(Γ)
end


"""
    Haction(Γ::Coefficients)::AbstractMatrix

Compute the Hessian of the constrained action for a given configuration ``Γ``.
"""
function Haction(p::Problem, Γ::Coefficients)::OffsetMatrix{Matrix{Matrix{Float64}}}
    _project = x -> project(p, x)
    _Haction = x -> Hkinetic(p, x) + Hpotential(p, x)
    (_project ∘ _Haction ∘ _project)(Γ)
end


"""
    kinetic(Γ::Coefficients)::Float64

Compute the kinetic part of the action for a given configuration ``Γ``.
"""
kinetic(p::Problem, Γ::Coefficients)::Float64 = 0.5 * Γ' * (p.K * Γ)

""" 
    ∇kinetic(Γ::Coefficients)::AbstractVector

Compute the gradient of the kinetic part of the action for a given configuration ``Γ``.
"""
∇kinetic(p::Problem, Γ::Coefficients)::Coefficients = p.K * Γ

""" 
    Hkinetic(Γ::Coefficients)::AbstractMatrix

Compute the Hessian of the kinetic part of the action for a given configuration ``Γ``.
"""
Hkinetic(p::Problem, _::Coefficients)::OffsetMatrix{Matrix{Matrix{Float64}}} = p.K 


"""
    potential(Γ::Coefficients)::Float64

Compute the potential part of the action for a given configuration ``Γ``.
"""
function potential(p::Problem, Γ::Coefficients)::Float64
    F, _, _ = dims(Γ)
    x = build_path(Γ, 2*F+1)
    V = U(x, p.m, p.f)

    # The potential part of the action is the discrete integral along the path using the
    # rectangle rule of the gradient potential energy (∇V) 

    # There are exactly steps+1 rectangles to integrate, each with width dt and height V[h].
    # We use the trapezoidal rule to integrate the potential energy along the path.
    # First, integrate between 1 and steps, thus omitting the initial and final points
    potential = sum(V[1:end-1])

    # The initial and final points contribute half of their value 
    potential += 0.5 * (V[0] + V[end])

    # The integral is the sum of the rectangles multiplied by the width of each rectangle
    return  potential * π / length(V)
end

""" 
    ∇potential(Γ::Coefficients)::AbstractVector

Compute the gradient of the potential part of the action for a given configuration ``Γ``.
"""
function ∇potential(p::Problem, Γ::Coefficients)::OffsetVector{Vector{Vector{Float64}}} 
    F, N, dim = dims(Γ)
    ∇potential = zeros(Γ, F+1, N, dim)
    x = build_path(Γ, 2*F+1)
    ∇V = ∇U(x, p.m, p.f)

    # The gradient of potential part of the action is the discrete integral along the path using the
    # rectangle rule of the gradient potential energy (∇V) times the derivative of the path x 
    # with respect to the Fourier coefficients Ak's (dx_dAk).
    # We use the trapezoidal rule to integrate the potential energy along the path.
    # There are exactly steps+1 rectangles to integrate, each with width dt and height ∇V[h].

    # First, integrate between 1 and steps, thus omitting the initial and final points
    ∇potential[0:F+1] = [sum(∇V[1:end-1] .* p.dx_dAk[k])  for k ∈ 0:F+1]

    # The initial and final points contribute half of their value
    ∇potential[0]   += 0.5 * ∇V[0]
    ∇potential[F+1] += 0.5 * ∇V[end]

    # The integral is the sum of the rectangles multiplied by the width of each rectangle
    return ∇potential * π / length(∇V)
end



"""
    Hpotential(Γ::Coefficients)::AbstractMatrix

Compute the Hessian of the potential part of the action for a given configuration ``Γ``.
"""
function Hpotential(p::Problem, Γ::Coefficients)::OffsetMatrix{Matrix{Matrix{Float64}}}
    F, N, dim = dims(Γ)
    Hpotential = OffsetMatrix([[zeros(dim, dim) for _ ∈ 1:N, _ ∈ 1:N] for _ ∈ 1:F+2, _ ∈ 1:F+2], 0:F+1, 0:F+1)
    x = build_path(Γ, 2*F+1)
    HV = HU(x, p.m, p.f)

    # The hessian of the potential part of the action is the discrete integral along the path using the
    # rectangle rule of the hessian potential energy (∇V) times the tensor product of the derivative of the path x 
    # with respect to the Fourier coefficients Ak's (dx_dAk) with itself. 
    # The second term arising from the prodduct rule, containing the second derivative of the path x 
    # with respect to the Fourier coefficients Ak's (d2x_dAk2) is zero because 
    # the path x is linear in the Fourier coefficients Ak's.

    # First, integrate between 1 and steps, thus omitting the initial and final points
    Hpotential[0:F+1, 0:F+1] =  [ sum((p.dx_dAk[k] .* p.dx_dAk[j]) .* HV[1:end-1]) for k ∈ 0:F+1, j ∈ 0:F+1]

    # The initial and final points contribute half of their value
    Hpotential[0, 0]   += 0.5 * HV[0]
    Hpotential[F+1, F+1] += 0.5 * HV[end]

    # The integral is the sum of the rectangles multiplied by the width of each rectangle
    return Hpotential * π / length(HV)
end


"""
    U(Γ::Coefficients, [n::Int64 = steps+1])::AbstractVector

Compute the potential for a given configuration ``Γ`` having an arbitrary function `f(r)` 
    at the denominator  and using  ``n`` points along the path.
"""
function U(x::Path, m::Vector{Float64}, f::Function = x -> x)::AbstractVector 
    steps, N, _ = dims(x)

    V = OffsetArray(zeros(steps+2), 0:steps+1)

    for h ∈ 0:steps+1, i ∈ 1:N-1, j ∈ (i+1):N
        V[h] += m[i] * m[j] /  f(norm(x[h][i] - x[h][j])) 
    end

    return V
end

"""
    ∇U(Γ::Coefficients, [n::Int64 = steps+1])::AbstractVector

Compute the gradient of the potential for a given configuration ``Γ``
 having an arbitrary function `f(r)` at the denominator and using  ``n`` points along the path.
"""
function ∇U(x::Coefficients, m::Vector{Float64}, f::Function = x -> x)::AbstractVector
    steps, N, dim = dims(x)

    ∇U = zeros(x, steps+1, N, dim)
    df = x -> derivative(f, x)
    

    for h ∈ 0:steps+1, i ∈ 1:N-1, j ∈ (i+1):N
        Δx = x[h][i] - x[h][j]
        r = norm(Δx)
        ∇U_ij = - m[i] * m[j] / f(r)^2 * df(r) * Δx / r
        ∇U[h][i] += ∇U_ij
        ∇U[h][j] -= ∇U_ij
    end

    return ∇U
end


""" 
    HU(Γ::Coefficients, [n::Int64 = steps+1])::AbstractMatrix

Compute the hessian for a given configuration ``Γ`` having an arbitrary function 
    `f(r)` at the denominator  and using  ``n`` points along the path.
"""
function HU(x::Path, m::Vector{Float64}, f::Function = x -> x)::AbstractVector 
    steps, N, dim = dims(x)

    HU = OffsetArray([[zeros(dim,dim) for _ ∈ 1:N, _ ∈ 1:N] for _ ∈ 0:steps+1], 0:steps+1)
    df = x -> derivative(f, x)
    d2f = x -> derivative(df, x)

    for h ∈ 0:steps+1, i ∈ 1:N-1, j ∈ (i+1):N
        Δx = x[h][i] - x[h][j]
        r = norm(Δx)

        HU_ij = - m[i] * m[j] / (f(r) * r)^2 * ( (Δx * Δx') * ( d2f(r) - df(r)/r - 2*df(r)^2 / f(r) ) + I * df(r) * r)

        HU[h][i,i] += HU_ij
        HU[h][j,j] += HU_ij
        HU[h][i,j] -= HU_ij
        HU[h][j,i] -= HU_ij
    end 

    return HU
end


"""
    K_energy(Γ::Coefficients, [n::Int64 = steps+1])::AbstractVector

Compute the kinetic energy for a given configuration ``Γ`` over ``n`` points along the path.
"""
function K_energy(Γ::Coefficients, n::Int64, m::Vector{Float64})::AbstractVector
    v = zeros(Γ, n)

    for h ∈ 0:n
        v[h] = (Γ[F+1] - Γ[0])/π + sum(k * Γ[k] * cos(k * h * π/n) for k ∈ 1:F)
    end
    
    Ek = [ 0.5 * sum( m[i] * v[h][i] ⋅ v[h][i]   for i∈1:N) for h ∈ 0:n]

    return FromZero(Ek)
end