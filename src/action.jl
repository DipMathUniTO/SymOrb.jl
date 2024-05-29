"""
    _action(Γ::Coefficients)::Float64

Compute the free action for a given configuration ``Γ``, without taking account the constraints.
"""
_action(Γ::Coefficients)::Float64  =  kinetic(Γ) +  potential(Γ)


"""
    action(Γ::Coefficients)::Float64

Compute the constrained action for a given configuration ``Γ``.
"""
action(Γ::Coefficients)::Float64 = (_action ∘ project)(Γ)

""" 
    action(v::Vector{Float64})::Float64

Compute the constrained action for a given flattened configuration ``v``. 
This function is used to interface with optimization algorithms that require a flat representation of the configuration.
"""
action(v::Vector{Float64})::Float64 =  (action ∘ emboss)(v)

"""
    _∇action(Γ::Coefficients)::AbstractVector

Compute the gradient of the free action for a given configuration ``Γ``, without taking account the constraints.
"""
_∇action(Γ::Coefficients)::AbstractVector = ∇kinetic(Γ) + ∇potential(Γ)

""" 
    ∇action(Γ::Coefficients)::AbstractVector

Compute the gradient of the constrained action for a given configuration ``Γ``.
"""
∇action(Γ::Coefficients)::AbstractVector = (project ∘ _∇action ∘ project)(Γ)

""" 
    ∇action(v::Vector{Float64})::Vector{Float64}

Compute the gradient of the constrained action for a given flattened configuration ``v``.
This function is used to interface with optimization algorithms that require a flat representation of the configuration.
"""
∇action(v::Vector{Float64})::Vector{Float64} = (flatten ∘ ∇action ∘ emboss)(v)



"""
    _Haction(Γ::Coefficients)::AbstractMatrix

Compute the Hessian of the free action for a given configuration ``Γ``, without taking account the constraints.
"""
_Haction(Γ::Coefficients)::AbstractMatrix = Hkinetic(Γ) + Hpotential(Γ)

"""
    Haction(Γ::Coefficients)::AbstractMatrix

Compute the Hessian of the constrained action for a given configuration ``Γ``.
"""
Haction(Γ::Coefficients)::AbstractMatrix = (project ∘ _Haction ∘ project)(Γ)

""" 
    Haction(v::Vector{Float64})::Matrix{Float64}

Compute the Hessian of the constrained action for a given flattened configuration ``v``.
This function is used to interface with optimization algorithms that require a flat representation of the configuration.
"""
Haction(v::Vector{Float64})::Matrix{Float64} = (flatten ∘ Haction ∘ emboss)(v)

"""
    kinetic(Γ::Coefficients)::Float64

Compute the kinetic part of the action for a given configuration ``Γ``.
"""
kinetic(Γ::Coefficients)::Float64 = 0.5 * Γ' * (K * Γ)

""" 
    ∇kinetic(Γ::Coefficients)::AbstractVector

Compute the gradient of the kinetic part of the action for a given configuration ``Γ``.
"""
∇kinetic(Γ::Coefficients)::AbstractVector = K * Γ

""" 
    Hkinetic(Γ::Coefficients)::AbstractMatrix

Compute the Hessian of the kinetic part of the action for a given configuration ``Γ``.
"""
Hkinetic(_::Coefficients)::AbstractMatrix = K 


"""
    potential(Γ::Coefficients)::Float64

Compute the potential part of the action for a given configuration ``Γ``.
"""
potential(Γ::Coefficients)::Float64 = begin
    V = U(Γ)

    # The potential part of the action is the discrete integral along the path using the
    # rectangle rule of the gradient potential energy (∇V) 

    # There are exactly steps+1 rectangles to integrate, each with width dt and height V[h].

    # First, integrate between 1 and steps, thus omitting the initial and final points
    potential = sum(V[1:steps])

    # For the remaining rectangle, use the mean of the initial and final points' contribution 
    potential += 0.5 * (V[0] + V[steps+1])

    # The integral is the sum of the rectangles multiplied by the width of each rectangle
    return  potential * dt
end

""" 
    ∇potential(Γ::Coefficients)::AbstractVector

Compute the gradient of the potential part of the action for a given configuration ``Γ``.
"""
∇potential(Γ) = begin
    ∇potential = OffsetArray([[zeros(dim) for _ ∈ 1:N] for _ ∈ 1:F+2], 0:F+1)

    ∇V = ∇U(Γ)

    # The gradient of potential part of the action is the discrete integral along the path using the
    # rectangle rule of the gradient potential energy (∇V) times the derivative of the path x 
    # with respect to the Fourier coefficients Ak's (dx_dAk).
        
    # There are exactly steps+1 rectangles to integrate, each with width dt and height ∇V[h].

    # First, integrate between 1 and steps, thus omitting the initial and final points
    ∇potential[0:F+1] = [sum(∇V[1:steps] .* dx_dAk[k])  for k ∈ 0:F+1]

    # For the remaining rectangle, use the mean of the initial and final points' contribution 
    ∇potential[0]   += 0.5 * ∇V[0]
    ∇potential[F+1] += 0.5 * ∇V[steps+1]

    # The integral is the sum of the rectangles multiplied by the width of each rectangle
    return ∇potential * dt
end



"""
    Hpotential(Γ::Coefficients)::AbstractMatrix

Compute the Hessian of the potential part of the action for a given configuration ``Γ``.
"""
Hpotential(Γ::Coefficients)::AbstractMatrix = begin

    Hpotential = OffsetMatrix([[zeros(dim, dim) for _ ∈ 1:N, _ ∈ 1:N] for _ ∈ 1:F+2, _ ∈ 1:F+2], 0:F+1, 0:F+1)

    HV = HU(Γ)

    # The hessian of the potential part of the action is the discrete integral along the path using the
    # rectangle rule of the hessian potential energy (∇V) times the tensor product of the derivative of the path x 
    # with respect to the Fourier coefficients Ak's (dx_dAk) with itself. 
    # The second term arising from the prodduct rule, containing the second derivative of the path x 
    # with respect to the Fourier coefficients Ak's (d2x_dAk2) is zero because 
    # the path x is linear in the Fourier coefficients Ak's.

    # First, integrate between 1 and steps, thus omitting the initial and final points
    Hpotential[0:F+1, 0:F+1] =  [ sum((dx_dAk[k] .* dx_dAk[j]) .* HV[1:steps]) for k ∈ 0:F+1, j ∈ 0:F+1]

    # For the remaining rectangle, use the mean of the initial and final points' contribution 
    Hpotential[0, 0]   += 0.5 * HV[0]
    Hpotential[F+1, F+1] += 0.5 * HV[steps+1]

    # The integral is the sum of the rectangles multiplied by the width of each rectangle
    return Hpotential * dt
end


"""
    U(Γ::Coefficients)::AbstractVector

Compute the potential for a given configuration ``Γ`` having an arbitrary function `f(r)` at the denominator
"""
U(Γ::Coefficients)::AbstractVector = begin
    V = OffsetArray(zeros(steps+2), 0:steps+1)

    x = build_path(Γ)

    for h ∈ 0:steps+1, i ∈ 1:N-1, j ∈ (i+1):N
        V[h] += m[i] * m[j] /  f(norm(x[h][i] - x[h][j])) 
    end

    return V
end

"""
    ∇U(Γ::Coefficients)::AbstractVector

Compute the gradient potential for a given configuration ``Γ`` having an arbitrary function `f(r)` at the denominator
"""
∇U(Γ::Coefficients)::AbstractVector = begin
    ∇U = OffsetArray([[zeros(dim) for _ ∈ 1:N] for _ ∈ 1:steps+2], 0:steps+1)
    
    x = build_path(Γ)

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
    HU(Γ::Coefficients)::AbstractMatrix

Compute the hessian for a given configuration ``Γ`` having an arbitrary function `f(r)` at the denominator
"""
HU(Γ) = begin
    HU = OffsetArray([[zeros(dim,dim) for _ ∈ 1:N, _ ∈ 1:N] for _ ∈ 1:steps+2], 0:steps+1)

    x = build_path(Γ)

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
