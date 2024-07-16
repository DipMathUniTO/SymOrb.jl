using Test, SHA 
using MinPath, OffsetArrays


# Test path.jl

@testset "All tests" verbose=true begin

configs = readdir("config", join=true)

function test_axes(M, ax)
    @test collect.(axes(M)) == collect.(ax)
end

function test_coefficients(Γ)
    @test Γ isa MinPath.Coefficients
    test_axes(Γ, (0:P.dims.F+1,))
    test_configuration(Γ[0])
end


function test_hessian(M)
    @test M isa OffsetArrays.OffsetMatrix{Matrix{Matrix{Float64}}}
    test_axes(M, (0:P.dims.F+1, 0:P.dims.F+1))
    test_hessian_configuration(M[0,0])
end

function test_hessian_configuration(M)
    test_axes(M, (1:P.dims.N, 1:P.dims.N))
    test_axes(M[1,1], (1:P.dims.dim, 1:P.dims.dim))
end


function test_path(x)
    @test x isa MinPath.Path
    test_axes(x, (0:P.dims.F*2+1, ))
    test_configuration(x[0])
end

function test_configuration(c)
    test_axes(c, (1:P.dims.N, ))
    test_axes(c[1], (1:P.dims.dim, ))
end


@testset "Load configuration files" begin
    for file in configs
        P = initialize(file)
        @test P isa MinPath.Problem
    end
end


P = initialize(configs[1])

@testset "Random starting path" begin
    Γ = MinPath.random_starting_path(P.dims)
    test_coefficients(Γ)
end

@testset "Retrieve dimensions" begin
    Γ = MinPath.random_starting_path(P.dims)
    F, N, dim = MinPath.dims(Γ)
    @test F == P.dims.F
    @test N == P.dims.N
    @test dim == P.dims.dim
end

@testset "Circular starting path" begin
    Γ = MinPath.circular_starting_path(P.dims)
    test_coefficients(Γ)
end

@testset "Perturb path" begin
    Γ = MinPath.circular_starting_path(P.dims)
    Γ = MinPath.perturbe_path(Γ)
    test_coefficients(Γ)
end

@testset "Perturbed circular starting path" begin
    Γ = MinPath.perturbed_circular_starting_path(P.dims)
    test_coefficients(Γ)
end

@testset "Fourier series" begin
    n = 2*P.dims.F + 1
    Γ = MinPath.random_starting_path(P.dims)
    s = MinPath.fourier_series(Γ, n)
    
    MinPath.inverse_fourier(s, P.dims.F)
    segment = MinPath.segment(Γ[0], Γ[end], n)
    test_path(segment)
    
    @test Γ[0] ≈ segment[0]
    @test Γ[end] ≈ segment[end]

    x1 = MinPath.build_path(Γ, n)
    test_path(x1)

    Γ1 = MinPath.fourier_coefficients(x1, P.dims.F)
    test_coefficients(Γ1)

    @test x1[0] ≈ Γ[0]
    @test x1[end] ≈ Γ[end]
    @test Γ ≈ Γ1 
end

@testset "Read and write path from file" begin
    tempdir = mktempdir()
    out_path = joinpath(tempdir, "test.txt")
    
    Γ = MinPath.circular_starting_path(P.dims)
    MinPath.print_path_to_file(Γ, out_path)
    
    @test sha2_256(open(out_path)) == sha2_256(open("./out/circular_test_1.txt"))
    
    Γ1 = MinPath.read_path_from_file(out_path, P.dims...)
    @test Γ ≈ Γ1

    rm(tempdir, recursive=true)
end


@testset "Flatten and emboss paths" begin
    Γ = MinPath.circular_starting_path(P.dims)
    v = MinPath.flatten(Γ)
    Γ1 = MinPath.emboss(v, P.dims)
    test_coefficients(Γ1)
    @test Γ ≈ Γ1
end

@testset "Action" begin
    n = P.dims.F*2+1

    Γ = MinPath.circular_starting_path(P.dims)

    x = MinPath.build_path(Γ, n)
    kin = MinPath.kinetic(P, Γ)
    @test kin isa Float64
    @test kin > 0
    U = MinPath.U(x, P.m, P.f)
    @test U isa OffsetArrays.OffsetVector{Float64}
    test_axes(U, (0:P.dims.F*2+1, ))

    potential = MinPath.potential(P, Γ)
    @test potential isa Float64

    action = MinPath.action(P, Γ)
    @test action isa Float64
end

@testset "Action gradient" begin
    n = P.dims.F*2+1

    Γ = MinPath.circular_starting_path(P.dims)
    x = MinPath.build_path(Γ, n)
    ∇kin = MinPath.∇kinetic(P, Γ)
    @test ∇kin isa MinPath.Coefficients
    test_coefficients(∇kin)

    ∇U = MinPath.∇U(x, P.m, P.f)
    @test ∇U isa OffsetArrays.OffsetVector{Vector{Vector{Float64}}}
    test_axes(∇U, (0:P.dims.F*2+1, ))

    ∇potential = MinPath.∇potential(P, Γ)
    test_coefficients(∇potential)

    ∇action = MinPath.∇action(P, Γ)
    test_coefficients(∇action)
end


@testset "Action hessian" begin
    n = P.dims.F*2+1
    Γ = MinPath.circular_starting_path(P.dims)
    x = MinPath.build_path(Γ, n)
    Hkin = MinPath.Hkinetic(P, Γ)
    test_hessian(Hkin)

    HU = MinPath.HU(x, P.m, P.f)
    @test HU isa OffsetArrays.OffsetVector{Matrix{Matrix{Float64}}}
    test_axes(HU, (0:P.dims.F*2+1, ))
    test_hessian_configuration(HU[0])

    Hpotential = MinPath.Hpotential(P, Γ)
    test_hessian(Hpotential)

    Haction = MinPath.Haction(P, Γ)
    test_hessian(Haction)
end

end