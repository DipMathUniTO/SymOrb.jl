function K_linear()
    K = OffsetArray(zeros(F + 2, F + 2), 0:F+1, 0:F+1)

    K[0, 0]     = 1 / π
    K[F+1, F+1] = 1 / π
    K[0, F+1]   = -1 / π
    K[F+1, 0]   = -1 / π

    for k ∈ 1:F
        K[k, k] = π / 2.0 * k^2
    end
    return K
end


function K_centrifugal()
    K = OffsetArray(zeros(F + 2, F + 2), 0:F+1, 0:F+1)

    K[0, 0]     = -π / 3.0
    K[F+1, F+1] = -π / 3.0
    K[0, F+1]   = -π / 3.0
    K[F+1, 0]   = -π / 3.0

    for k ∈ 1:F
        K[0, k], K[k, 0] = 1.0 / k
        K[F+1, k], K[k, F+1]  = 1.0 / k * (-1)^k
        K[k, k] = π / 2.0
    end

    return K
end

function K_coriolis()
    K = OffsetArray(zeros(F + 2, F + 2), 0:F+1, 0:F+1)

    K[0, F+1] = 1.0
    K[F+1, 0] = -1

    for k ∈ 2:2:F #only the even ones
        K[0, k] = 2.0 / (k * π)
        K[k, 0] = -K[0, k]
        K[F+1, k] = -2.0 / (k * π)
        K[k, F+1] = -K[k, F+1]

        for j ∈ 2:2:F
            if k == j continue end
            K[k, j] = (k^2 + j^2) / (k^2 - j^2)
            K[j, k] = -K[k, j]
        end
    end

    return K
end



∫sin(k) = (1 - (-1)^k) / k
∫tsin(k) = -pi * (-1)^k / k
∫tcos(k) = (-1 + (-1)^k) / k^2
∫coscos(k, h) = ifelse(h == k, pi / 2, 0)
∫sinsin(k, h) = ifelse(h == k, pi / 2, 0)
∫sincos(k, h) = ifelse(h == k, 0, k / (k^2 - h^2) * (1 - (-1)^(k + h)))

