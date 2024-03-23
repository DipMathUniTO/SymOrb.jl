action(Γ)::Float64 = @timeit to "action" kinetic(Γ) + potential(Γ)

kinetic(Γ) = if (isΩ)  K_linear(Γ) + K_centrifugal(Γ) + K_coriolis(Γ) else  K_linear(Γ) end

function K_linear(Γ)
  K = 0
  for i ∈ 1:N
    Ki = 1/π * (Γ[0][i] - Γ[F+1][i])' * (Γ[0][i] - Γ[F+1][i])
    for k ∈ 1:F
        Ki += π/2. * k^2 * Γ[k][i]' * Γ[k][i]
    end
    K += 0.5*m[i]*Ki
  end
  return K
end

function K_centrifugal(Γ)
    K = 0
    for i ∈ 1:N
        Ki = -π/3. * ( Γ[0][i] +  Γ[F+1][i])' * Ω2 *(Γ[0][i] +  Γ[F+1][i])
        for k ∈ 1:F
             Ki += 2.0/k * (Γ[0][i]' * Ω2 * Γ[k][i]  + (-1)^k * Γ[F+1][i]' * Ω2 * Γ[k][i]) +  π/2. * Γ[k][i]' * Ω2 * Γ[k][i]
        end
        K += 0.5*m[i]*Ki
    end
    return K
end

function K_coriolis(Γ)
    K = 0
    for i ∈ 1:N
        Ki = Γ[0][i]' * Ω * Γ[F+1][i] - Γ[F+1][i]' * Ω * Γ[0][i]
        for k ∈ 1:F
            if isodd(k) continue end
            Ki += 4.0/(k*π) * (Γ[0][i]'*Ω*Γ[k][i] - Γ[k][i]'*Ω*Γ[0][i] - Γ[F+1][i]'*Ω*Γ[k][i] + Γ[k][i]'*Ω*Γ[F+1][i])

            for j ∈ 1:F
                if isodd(k + j) || k == j continue end
                Ki += (k^2 + j^2)/(k^2 - j^2) * Γ[k][i]' * Γ[j][i]
            end
        end
        K += 0.5*m[i]*Ki
    end
    return K
end

function potential(Γ)
    x = build_path(Γ)

    V = 0
    for i ∈ 1:N-1
        for j ∈ (i+1):N
            V += sum(m[i]*m[j] / norm(x[t][i] - x[t][j]) for t ∈ 0:steps+1)
        end
    end
    return V
end


function v_action(v)

    @timeit to "emboss" begin
        Γ = emboss(v)
    end
    Γ = project(Γ)
    @timeit to "emboss" begin
        v = flatten(Γ)
    end
    return action(Γ)
end



∫sin(k) =(1-(-1)^k) / k
∫tsin(k) = -pi*(-1)^k / k
∫tcos(k) = (-1 + (-1)^k) / k^2
∫coscos(k, h) = ifelse(h == k, pi / 2, 0)
∫sinsin(k, h) = ifelse(h == k, pi / 2, 0)
∫sincos(k, h) = ifelse(h == k, 0, k / (k^2 - h^2) * (1 - (-1)^(k + h)))

