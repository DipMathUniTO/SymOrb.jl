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


function K_centrifugal()
    K = OffsetArray(zeros(F + 2, F + 2), 0:F+1, 0:F+1)

    K[0, 0]     = -π / 3.0 * Id
    K[F+1, F+1] = K[0, 0] 
    K[0, F+1]   = K[0, 0]
    K[F+1, 0]   = K[0, 0]

    for k ∈ 1:F
        K[0, k]   = 1.0 / k * Id
        K[k, 0]   = K[0, k]
        K[F+1, k] = K[0, k]* (-1)^k
        K[k, F+1] = K[F+1, k]
        K[k, k]   = π / 2.0 * Id
    end

    return K
end

function K_coriolis()
    K = OffsetArray(zeros(F + 2, F + 2), 0:F+1, 0:F+1)

    K[0, F+1] = 1.0 * Id
    K[F+1, 0] = -1 * Id

    for k ∈ 2:2:F #only the even ones
        K[0, k] = 2.0 / (k * π) * Id
        K[k, 0] = - K[0, k] 
        K[F+1, k] = K[k, 0]
        K[k, F+1] = K[0, k] 

        for j ∈ 2:2:F
            if k == j continue end
            K[k, j] = (k^2 + j^2) / (k^2 - j^2) * Id
            K[j, k] = -K[k, j]
        end
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
     