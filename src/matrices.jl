"""
    K_linear(N, F, dim, m)

Compute the linear part of the kinetic energy operator

# Arguments
- `N::Int`: The number of particles
- `F::Int`: The number of Fourier series terms
- `dim::Int`: The dimension of the space
- `m::Matrix{Matrix{Float64}}`: The masses of the particles
"""
function K_linear(M, dims)
    F, N, dim = dims
    K = zeros(dim, N, F+2, dim, N, F+2)
    
    K[:, :,  1 , :, :,  1] = 1 / π * M
    K[:, :, F+2, :, :, F+2] =  K[:, :, 1, :, :, 1]
    K[:, :, F+2, :, :,  1]   = -K[:, :, 1, :, :, 1]
    K[:, :,  1 , :, :, F+2]   = -K[:, :, 1, :, :, 1]

    for k ∈ 1:F
        K[:, :, k+1, :, :, k+1] = π / 2.0 * k^2 * M
    end
    return reshape(K, (F+2)*N*dim, (F+2)*N*dim)
end


"""
    K_centrifugal(Ω2, N, F, dim, m)

Compute the centrifugal part of the kinetic energy operator

# Arguments
- `Ω2::Float64`: The square of the generator of the rotation
- `N::Int`: The number of particles
- `F::Int`: The number of Fourier series terms
- `dim::Int`: The dimension of the space
- `m::Matrix{Matrix{Float64}}`: The masses of the particles
"""
function K_centrifugal(Ω2, M, dims)
    F, N, dim = dims
    K = zeros(dim, N, F+2, dim, N, F+2)

    K[:, :,  1,  :, :,  1]  = -π / 3.0 * M
    K[:, :, F+2, :, :, F+2] = K[:, :, 1, :, :, 1] 
    K[:, :,  1,  :, :, F+2] = K[:, :, 1, :, :, 1]
    K[:, :, F+2, :, :,  1]  = K[:, :, 1, :, :, 1]

    for k ∈ 2:F+1
        K[:, :,  1, :, :, k]   = 1.0 / k * M
        K[:, :,  k, :, :, 1]   = 1.0 / k * M 
        
        K[:, :, F+2, :, :,  k]  = (-1)^k / k * M
        K[:, :,  k , :, :, F+2] = (-1)^k / k * M
        
        K[:, :, k, :, :, k] = π / 2.0 * M
    end

    for k ∈ 1:F+2, h ∈ 1:F+2, i ∈ 1:N, j ∈ 1:N
        K[:, i, k, :, j, h] *= Ω2
    end

    return reshape(K, (F+2)*N*dim, (F+2)*N*dim)
end


"""
    K_coriolis(Ω, N, F, dim, M)

Compute the Coriolis part of the kinetic energy operator

# Arguments
- `Ω::Float64`: The generator of the rotation
- `N::Int`: The number of particles
- `F::Int`: The number of Fourier series terms
- `dim::Int`: The dimension of the space
- `m::Matrix{Matrix{Float64}}`: The masses of the particles
"""
function K_coriolis(Ω, M, dims)
    F, N, dim = dims
    K = zeros(dim, N, F+2, dim, N, F+2)
  
    K[:, :,  1 , :, F+2] = -1.0 * M
    K[:, :, F+2, :,  1 ] = 1.0 * M

    for k ∈ 1:F
        K[:, :, 1, :, :, k] = 2.0 * ( (-1)^k - 1) / (k * π) * M
        K[:, :, k, :, :, 1] = - K[:, :, 1, :, :, k] 
        
        K[:, :, F+2, :, F+2, k ] = -K[:, :, 1, :, :, k]
        K[:, :,  k,  :,  k, F+2] =  K[:, :, 1, :, :, k] 

        for j ∈ 2:F+1
            if k == j continue end
            K[:, :, k, :, :, j] = 2*((-1)^(j+k) -1) *j*k / (j^2 - k^2) * M
            K[:, :, j, :, :, k] = -K[:, :, k, :, :, j]
        end
    end


    for k ∈ 1:F+2, h ∈ 1:F+2, i ∈ 1:N, j ∈ 1:N
        K[:, i, k, :, j, h] *= Ω
    end

    return reshape(K, (F+2)*N*dim, (F+2)*N*dim)
end


"""
    compute_dx_dAk(F, steps)

Compute the derivative of the path with respect to the Fourier coefficients

# Arguments
- `F::Int`: The number of Fourier series terms
- `steps::Int`: The number of steps in the path
"""
function compute_dx_dA(F, steps)
    M = zeros(steps+2, F+2)

    M[:, 1] = [  (1 - h  / (steps+1) ) for h ∈ 0:steps+1]
    M[:, F+2] = [ (h / (steps+1) ) for h ∈ 0:steps+1]

    M[ 2:steps+1, 2:F+1] = [ sin(k * h * π / (steps+1))  for h ∈ 1:steps, k ∈ 1:F] 
    return  M      
end 

function compute_dA_dx(F, steps)
    M =  zeros(F+2, steps+2)
    
    M[1, 1] = 1.0
    M[2:F+1, 1] = [- 2/(k*π) for k ∈ 1:F] 

    M[F+2, steps+2] =  1.0
    M[2:F+1, steps+2] = [ 2/(k*π)*(-1)^k for k ∈ 1:F]
    
    M[2:F+1, 2:steps+1] = [  2/(steps+1) * sin(k * h * π / (steps+1)) for k ∈ 1:F, h ∈ 1:steps]

    return M
end


function kinetic_matrix(Ω, m, dims)
    _, N, dim = dims

    M = zeros(dim, N, dim, N)

    for i ∈ 1:N
        M[:, i, :, i] = m[i] * I(dim)
    end

    K = K_linear(M, dims)          # linear part of the kinetic energy matrix 
    
    if  (!iszero(Ω))
        K .+=  K_centrifugal(Ω*Ω, M, dims) + K_coriolis(Ω, M, dims)
    end
    
    return K
end



function compute_integration_factors(steps)
    factors = ones(steps+2)
    factors[1] = 0.5
    factors[end] = 0.5
    factors * π / (steps+1)
end