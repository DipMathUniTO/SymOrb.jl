"""
    K_linear(N, F, dim, m)

Compute the linear part of the kinetic energy operator

# Arguments
- `N::Int64`: The number of particles
- `F::Int64`: The number of Fourier series terms
- `dim::Int64`: The dimension of the space
- `m::Matrix{Matrix{Float64}}`: The masses of the particles
"""
function K_linear(N, F, dim, M)
    K = OffsetArray([[zeros(dim, dim) for _ ∈ 1:N, _ ∈ 1:N ] for _ ∈ 0:F+1, _ ∈ 0:F+1] , 0:F+1, 0:F+1)
    
    K[0, 0] = 1 / π * M
    K[F+1, F+1] =  K[0, 0]
    K[0, F+1]   = -K[0, 0]
    K[F+1, 0]   = -K[0, 0]

    for k ∈ 1:F
        K[k, k] = π / 2.0 * k^2 * M
    end
    return K
end


"""
    K_centrifugal(Ω2, N, F, dim, m)

Compute the centrifugal part of the kinetic energy operator

# Arguments
- `Ω2::Float64`: The square of the generator of the rotation
- `N::Int64`: The number of particles
- `F::Int64`: The number of Fourier series terms
- `dim::Int64`: The dimension of the space
- `m::Matrix{Matrix{Float64}}`: The masses of the particles
"""
function K_centrifugal(Ω2, N, F, dim, M)
    K = OffsetArray([[zeros(dim, dim) for _ ∈ 1:N, _ ∈ 1:N ] for _ ∈ 0:F+1, _ ∈ 0:F+1] , 0:F+1, 0:F+1)

    K[0, 0]     = -π / 3.0 * M
    K[F+1, F+1] = K[0, 0] 
    K[0, F+1]   = K[0, 0]
    K[F+1, 0]   = K[0, 0]

    for k ∈ 1:F
        K[0, k]   = 1.0 / k * M
        K[k, 0]   = 1.0 / k * M 
        K[F+1, k] = (-1)^k / k * M
        K[k, F+1] = (-1)^k / k * M
        K[k, k]   = π / 2.0 * M
    end

    for k ∈ 0:F+1, h ∈ 0:F+1, i ∈ 1:N, j ∈ 1:N
        K[k, h][i, j] *= Ω2
    end

    return K
end


"""
    K_coriolis(Ω, N, F, dim, M)

Compute the Coriolis part of the kinetic energy operator

# Arguments
- `Ω::Float64`: The generator of the rotation
- `N::Int64`: The number of particles
- `F::Int64`: The number of Fourier series terms
- `dim::Int64`: The dimension of the space
- `m::Matrix{Matrix{Float64}}`: The masses of the particles
"""
function K_coriolis(Ω, N, F, dim, M)
    K =  OffsetArray([[zeros(dim, dim) for _ ∈ 1:N, _ ∈ 1:N ] for _ ∈ 0:F+1, _ ∈ 0:F+1] , 0:F+1, 0:F+1)
  
    K[0, F+1] = -1.0 * M
    K[F+1, 0] = 1.0 * M

    for k ∈ 1:F
        K[0, k] = 2.0 * ( (-1)^k - 1) / (k * π) * M
        K[k, 0] = - K[0, k] 
        K[F+1, k] =-  K[0, k]
        K[k, F+1] = K[0, k] 

        for j ∈ 1:F
            if k == j continue end
            K[k, j] = 2*((-1)^(j+k) -1) *j*k / (j^2 - k^2) * M
            K[j, k] = -K[k, j]
        end
    end


    for k ∈ 0:F+1, h ∈ 0:F+1, i ∈ 1:N, j ∈ 1:N
        K[k, h][i, j] *= Ω
    end

    return K
end


"""
    compute_dx_dAk(F, steps)

Compute the derivative of the path with respect to the Fourier coefficients

# Arguments
- `F::Int64`: The number of Fourier series terms
- `steps::Int64`: The number of steps in the path
"""
function compute_dx_dAk(F, steps)
    dx_dAk = OffsetArray([zeros(steps) for _ in 0:F+1], 0:F+1)

    dx_dAk[0] = [  (1 - h  / (steps+1) ) for h ∈ 1:steps]
    dx_dAk[F+1] = [  (h / (steps+1) ) for h ∈ 1:steps]
    dx_dAk[1:F] = [  [sin(k * h * π / (steps+1)) for h ∈ 1:steps] for k ∈ 1:F]

    return dx_dAk       
end 