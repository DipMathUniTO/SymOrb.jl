_action(Γ)::Float64 = kinetic(Γ) + potential(Γ)
_∇action(Γ) = ∇kinetic(Γ) + ∇potential(Γ)

kinetic(Γ) = 0.5 * Γ' * (K * Γ)
∇kinetic(Γ) = K * Γ

action(v) = (_action ∘ project ∘ emboss)(v)
∇action(v) = (flatten ∘ project ∘ _∇action ∘ emboss)(v)

potential(Γ) = begin
    x = build_path(Γ)
    V = 0
    for i ∈ 1:N-1, j ∈ (i+1):N
        V += sum(m[i] * m[j] / norm(x[h][i] - x[h][j]) for h ∈ 0:steps+1)
    end
    return V
end


∇U(Γ) = begin
    x = build_path(Γ)

    ∇V =OffsetArray([[zeros(dim) for _ ∈ 1:N] for _ ∈ 1:steps+2], 0:steps+1)

    for h ∈ 0:steps+1, i ∈ 1:N-1, j ∈ (i+1):N
        r = x[h][i] - x[h][j]
        ∇V_ij = -m[i] * m[j] * r / norm(r)^3
        ∇V[h][i] += ∇V_ij
        ∇V[h][j] -= ∇V_ij
    end
    return ∇V
end


∇potential(Γ) = begin
    ∇V = ∇U(Γ)

    ∇potential = OffsetArray([[zeros(dim) for _ in 1:N] for _ in 1:F+2], 0:F+1)

    for h ∈ 0:steps+1
        ∇potential[0] += ∇V[h] * (1 - h / (steps + 1))
        ∇potential[F+1] += ∇V[h] * (h / (steps + 1))
        ∇potential[1:F] .+= [∇V[h] * sin(π * k * h / (steps + 1)) for k ∈ 1:F]
    end

    return ∇potential
end

