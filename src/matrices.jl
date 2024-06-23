function K_linear()
    K = OffsetArray([[zeros(dim, dim) for _ ∈ 1:N, _ ∈ 1:N ] for _ ∈ 0:F+1, _ ∈ 0:F+1] , 0:F+1, 0:F+1)

    K[0, 0] = 1 / π * Id
    K[F+1, F+1] =  K[0, 0]
    K[0, F+1]   = -K[0, 0]
    K[F+1, 0]   = -K[0, 0]

    for k ∈ 1:F
        K[k, k] = π / 2.0 * k^2 * Id
    end
    return K
end


function K_centrifugal(Ω2)
    K = OffsetArray([[zeros(dim, dim) for _ ∈ 1:N, _ ∈ 1:N ] for _ ∈ 0:F+1, _ ∈ 0:F+1] , 0:F+1, 0:F+1)

    K[0, 0]     = -π / 3.0 * Id
    K[F+1, F+1] = K[0, 0] 
    K[0, F+1]   = K[0, 0]
    K[F+1, 0]   = K[0, 0]

    for k ∈ 1:F
        K[0, k]   = 1.0 / k * Id
        K[k, 0]   = 1.0 / k * Id 
        K[F+1, k] = (-1)^k / k * Id
        K[k, F+1] = (-1)^k / k * Id
        K[k, k]   = π / 2.0 * Id
    end

    for k ∈ 0:F+1, h ∈ 0:F+1, i ∈ 1:N, j ∈ 1:N
        K[k, h][i, j] *= Ω2
    end

    return K
end

function K_coriolis(Ω)
    K =  OffsetArray([[zeros(dim, dim) for _ ∈ 1:N, _ ∈ 1:N ] for _ ∈ 0:F+1, _ ∈ 0:F+1] , 0:F+1, 0:F+1)

    K[0, F+1] = -1.0 * Id
    K[F+1, 0] = 1.0 * Id

    for k ∈ 1:F
        K[0, k] = 2.0 * ( (-1)^k - 1) / (k * π) * Id
        K[k, 0] = - K[0, k] 
        K[F+1, k] =-  K[0, k]
        K[k, F+1] = K[0, k] 

        for j ∈ 1:F
            if k == j continue end
            K[k, j] = 2*((-1)^(j+k) -1) *j*k / (j^2 - k^2) * Id
            K[j, k] = -K[k, j]
        end
    end


    for k ∈ 0:F+1, h ∈ 0:F+1, i ∈ 1:N, j ∈ 1:N
        K[k, h][i, j] *= Ω
    end

    return K
end

function compute_dx_dAk()
    dx_dAk = OffsetArray([zeros(steps) for _ in 0:F+1], 0:F+1)

    dx_dAk[0] = [  (1 - h * dt / π) for h ∈ 1:steps]
    dx_dAk[F+1] = [  (h * dt / π) for h ∈ 1:steps]
    dx_dAk[1:F] = [  [sin(k * h * dt) for h ∈ 1:steps] for k ∈ 1:F]

    return dx_dAk       
end 
     