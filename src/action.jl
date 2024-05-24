_action(Γ)  =  kinetic(Γ) +  potential(Γ)
_∇action(Γ) = ∇kinetic(Γ) + ∇potential(Γ)
_Haction(Γ) = Hkinetic(Γ) + Hpotential(Γ)

action(v) =  (_action ∘ project ∘ emboss)(v)
∇action(v) = (flatten ∘ project ∘ _∇action ∘ project ∘ emboss)(v)
Haction(v) = (flatten ∘ project ∘ _Haction ∘ project ∘ emboss)(v)

# The kinetic part of the action and its gradient
kinetic(Γ) = 0.5 * Γ' * (K * Γ)
∇kinetic(Γ) = K * Γ
Hkinetic(_) = K 

function myfun!(F, J, v)

    if ! isnothing(F)
        F = ∇action(v)
    end
    if ! isnothing(J)
        J = Haction(v)
    end
end

# The potential part of the action
potential(Γ) = begin
    V = U(Γ)
    potential = sum(V[1:steps])
    potential += 0.5 * (V[0] + V[steps+1])
    return  potential * dt
end

# The gradient of the potential part of the action
∇potential(Γ) = begin
    ∇potential = OffsetArray([[zeros(dim) for _ ∈ 1:N] for _ ∈ 1:F+2], 0:F+1)

    ∇V = ∇U(Γ)

    ∇potential[0:F+1] = [sum(∇V[1:steps] .* dx_dAk[k])  for k ∈ 0:F+1]

    ∇potential[0]   += 0.5 * ∇V[0]
    ∇potential[F+1] += 0.5 * ∇V[steps+1]

    return ∇potential * dt
end

# The Hessian of the potential part of the action
Hpotential(Γ) = begin

    Hpotential = OffsetMatrix([[zeros(dim, dim) for _ ∈ 1:N, _ ∈ 1:N] for _ ∈ 1:F+2, _ ∈ 1:F+2], 0:F+1, 0:F+1)

    HV = HU(Γ)

    Hpotential[0:F+1, 0:F+1] =  [ sum((dx_dAk[k] .* dx_dAk[j]) .* HV[1:steps]) for k ∈ 0:F+1, j ∈ 0:F+1]

    Hpotential[0, 0]   += 0.5 * HV[0]
    Hpotential[F+1, F+1] += 0.5 * HV[steps+1]

    return Hpotential * dt
end

# The potential
U(Γ) = begin
    V = OffsetArray(zeros(steps+2), 0:steps+1)

    x = build_path(Γ)

    for h ∈ 0:steps+1, i ∈ 1:N-1, j ∈ (i+1):N
        V[h] += m[i] * m[j] /  f(norm(x[h][i] - x[h][j])) 
    end

    return V
end

# The gradient of the potential
∇U(Γ) = begin
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


# The Hessian of the potential
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


# f(x) = sqrt(x^2+ϵ)
# df(x) = x/sqrt(x^2 + ϵ)
# d2f(x) = 0.01/(x^2 + ϵ)^(3/2)

f(x) = x
df(x) = 1
d2f(x) = 0
