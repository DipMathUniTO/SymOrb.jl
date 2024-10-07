"""
    action(p::Problem, Γ::Vector{T})::T
Compute the constrained action for a given configuration ``Γ``.
"""

function action(p::Problem, Γ::Vector{T})::T where {T}
    Γ .= p.Π * Γ    # Project the path
    kinetic(p, Γ) + potential(p, Γ)    
end


""" 
    ∇action(p::Problem, Γ::Vector{T})::Vector{T}

Compute the gradient of the constrained action for a given configuration ``Γ``.
"""
function ∇action(p::Problem, Γ::Vector{T})::Vector{T} where {T}
    Γ .= p.Π * Γ
    p.Π' * (∇kinetic(p, Γ) + ∇potential(p, Γ))
end



"""
    Haction(p::Problem, Γ::Vector{T})::Matrix{T}

Compute the Hessian of the constrained action for a given configuration ``Γ``.
"""
function Haction(p::Problem, Γ::Vector{T})::Matrix{T} where T
    Γ .= p.Π * Γ
    p.Π' * (Hkinetic(p, Γ) + Hpotential(p, Γ)) * p.Π
end


"""
    kinetic(p::Problem, Γ::Vector{T})::T

Compute the kinetic part of the action for a given configuration ``Γ``.
"""
function kinetic(p::Problem, Γ::Vector{T})::T where {T} 
    0.5 * Γ' * (p.K * Γ)
end

""" 
    ∇kinetic(p::Problem, Γ::Vector{T})::Vector{T}

Compute the gradient of the kinetic part of the action for a given configuration ``Γ``.
"""
function ∇kinetic(p::Problem, Γ::Vector{T})::Vector{T} where {T}
    p.K * Γ
end

""" 
    Hkinetic(Γ::Coefficients)::AbstractMatrix

Compute the Hessian of the kinetic part of the action for a given configuration ``Γ``.
"""
function Hkinetic(p::Problem, _::Vector{T})::Matrix{T} where {T} 
    p.K
end

"""
    potential(p::Problem, Γ::Vector{T})::T

Compute the potential part of the action for a given configuration ``Γ``.
"""
function potential(p::Problem, Γ::Vector{T})::T where{T}
    x = build_path(p, Γ)

    # The vector p.Z represents the linear transformation from U(t) to ∫U(t) dt
    p.Zu ⋅ U(p, x)
end

""" 
    ∇potential(p::Problem, Γ::Vector{T})::Vector{T}

Compute the gradient of the potential part of the action for a given configuration ``Γ``.
"""
function ∇potential(p::Problem, Γ::Vector{T})::Vector{T} where T 
    x = build_path(p, Γ)
    
    # The gradient of potential part of the action is the discrete integral along the path using the
    # rectangle rule of the gradient potential energy (∇V) times the derivative of the path x 
    # with respect to the Fourier coefficients Ak's (dx_dAk).
    
    p.Zg * ∇U(p, x)
end


"""
    Hpotential(p::Problem, Γ::Vector{T})

Compute the Hessian of the potential part of the action for a given configuration ``Γ``.
"""
function Hpotential(p::Problem, Γ::Vector{T}) where T
    
    x = build_path(p, Γ)


    # The hessian of the potential part of the action is the discrete integral along the path using the
    # rectangle rule of the hessian potential energy (∇V) times the tensor product of the derivative of the path x 
    # with respect to the Fourier coefficients Ak's (dx_dAk) with itself. 
    # The second term arising from the prodduct rule, containing the second derivative of the path x 
    # with respect to the Fourier coefficients Ak's (d2x_dAk2) is zero because 
    # the path x is linear in the Fourier coefficients Ak's.

    # The matrix p.Zl and p.A_to_x are the matrices that represents the linear transformation 
    # from HU(t) to ∫(∂x/∂A)' HU(t) ∂x/∂A dt

    p.Zg * HU(p, x) * p.A_to_x 

end


function U_t( x::AbstractArray{T, 2}, m::Vector, f::Function) where T
    N = size(x, 2)
    pot = 0.0
    for i ∈ 1:N-1, j ∈ (i+1):N
        pot += m[i] * m[j] / f(norm(x[:, i] - x[:, j]))
    end
    return pot
end


"""
   U(P, x::Array{T, 3})::Vector{T}

Compute the potential at every time step along the path ``x``.
"""
function U(P, x::Array{T, 3})::Vector{T} where T
    N = size(x, 2)
    V = zeros(T, size(x, 3))
    for h ∈ axes(x,3)
        @views V[h] += U_t(x[:, :, h], P.m, P.f)
    end

    return V
end



function ∇U_t!(grad::AbstractArray{T, 2}, x::AbstractArray{T, 2}, m::Vector, f::Function, df::Function) where T
    N = size(x, 2)
    for i ∈ 1:N-1, j ∈ (i+1):N
        Δx = x[:, i] - x[:, j]
        r = norm(Δx)
        ∇U_ij = - m[i] * m[j] / f(r)^2 * df(r) * Δx / r
        @views @. grad[:, i] += ∇U_ij
        @views @. grad[:, j] -= ∇U_ij
    end
end



"""
    ∇U(Γ::Coefficients, [n::Int = steps+1])::AbstractVector

Compute the gradient of the potential at every time step along the path ``x``.
"""
function ∇U(P::Problem, x::Array{T, 3})::Vector{T} where {T}
    dim, N, steps  = size(x)
    ∇U =  zeros(T, dim, N, steps)

    df = x -> derivative(P.f, x)
    for h ∈ 1:steps
        @views ∇U_t!(∇U[:, :, h], x[:, :, h], P.m, P.f, df)
    end

    reshape(∇U, dim*N*steps)
end


function HU_t!(hessian::AbstractArray{T, 4},  x::AbstractArray{T, 2}, m::Vector, f::Function, df::Function, d2f::Function) where T
    N = size(x, 2)

    for i ∈ 1:N-1, j ∈ (i+1):N
        Δx = x[:, i] - x[:, j]
        r = norm(Δx)

        HU_ij = - m[i] * m[j] / (f(r) * r)^2 * ( (Δx * Δx') * ( d2f(r) - df(r)/r - 2*df(r)^2 / f(r) ) + I * df(r) * r)

        @views @. hessian[:, i, :, i] += HU_ij
        @views @. hessian[:, j, :, j] += HU_ij
        @views @. hessian[:, i, :, j] -= HU_ij
        @views @. hessian[:, j, :, i] -= HU_ij
    end
end



""" 
    HU(Γ::Coefficients, [n::Int = steps+1])::AbstractMatrix

Compute the hessian of the potential at every time step along the path ``x``
"""
function HU(P::Problem, x::Array{T, 3})::Matrix{T} where {T} 
    
    dim, N, steps = size(x)

    HU = zeros(T, dim, N, steps, dim, N, steps)
    df = x -> derivative(P.f, x)
    d2f = x -> derivative(df, x)
    
    for h ∈ axes(x, 3)
        HU_t!((@view HU[:, :, h, :, :, h]), (@view x[:, :, h]), P.m, P.f, df, d2f)
    end 

    reshape(HU, dim*N*steps, dim*N*steps)
end



"""
    K_energy(P::Problem, Γ::Vector{T}, n::Int)::Vector{T} 

Compute the kinetic energy for a given configuration ``Γ`` over ``n`` points along the path.
"""
function K_energy(P::Problem, Γ::Vector{T}, n::Int)::Vector{T} where T
    Γ = P.R * Γ

    F = size(Γ, 1) ÷ (P.dim * P.N) - 2 

    Γ = reshape(Γ, P.dim, P.N, F+2)
    v = zeros(P.dim, P.N, n)

    for h ∈ 1:n
        v[:, :, h] = (Γ[:, :, end] - Γ[:, :, 1])/π + sum(k * Γ[:, :, k+1] * cos(k * (h-1) * π/(n-1)) for k ∈ 1:F)
    end
    
    Ek = [ 0.5 * sum( P.m[i] * v[:, i, h] ⋅ v[:, i, h] for i ∈ 1:P.N) for h ∈ 1:n]

    return Ek
end