using Test, SHA 
using SymOrb


# Test path.jl

@testset "All tests" verbose=true begin

configs = readdir("config", join=true)


function test_coefficients(Γ)
    @test Γ isa Vector
    @test size(Γ) == (dims.dim*(dims.N-1)*(dims.F+2),)
end


@testset "Load configuration files" begin
    for file in configs
        P = initialize(file)
        @test P isa SymOrb.Problem
        @test P.Π^2 ≈ P.Π atol=1e-3
    end
end


P = initialize(configs[1])
dims = (F = P.F, N = P.N, dim = P.dim)

@testset "Random starting path" begin
    Γ = SymOrb.random_starting_path(P)
    test_coefficients(Γ)
end

@testset "Circular starting path" begin
    Γ = SymOrb.circular_starting_path(P)
    test_coefficients(Γ)
end

@testset "Perturbe path" begin
    Γ = SymOrb.circular_starting_path(P)
    Γ = SymOrb.perturbe_path(P, Γ)
    test_coefficients(Γ)
end

@testset "Perturbed circular starting path" begin
    Γ = SymOrb.perturbed_circular_starting_path(P)
    test_coefficients(Γ)
end

@testset "Fourier series" begin
  
    Γ = SymOrb.random_starting_path(P)

    x1 = SymOrb.build_path(P, Γ, 500)
    @test x1 isa Array{Float64, 3}

    Γ1 = SymOrb.fourier_coefficients(P, x1, 500)
    test_coefficients(Γ1)

    @test Γ ≈ Γ1  atol = 1e-3
end

@testset "Read and write path from file" begin
    tempdir = mktempdir()
    out_path = joinpath(tempdir, "test.toml")
    
    Γ = SymOrb.circular_starting_path(P)
    SymOrb.print_path_to_file(P, Γ, Γ, out_path)
    
    P1, Γ1 = SymOrb.read_path_from_file(out_path)

    @test Γ ≈ Γ1  atol = 1e-3
    rm(tempdir, recursive=true)
end


@testset "Action" begin
    
    Γ = SymOrb.random_starting_path(P)

    x = SymOrb.build_path(P, Γ)
    kin = SymOrb.kinetic(P, Γ)
    @test kin isa Float64
    @test kin > 0
    U = SymOrb.U(P, x)
    @test U isa Vector{Float64}
    @test size(U) == (P.steps+2, )

    potential = SymOrb.potential(P, Γ)
    @test potential isa Float64

    action = SymOrb.action(P, Γ)
    @test action isa Float64
end

@testset "Action gradient" begin

    Γ = SymOrb.random_starting_path(P)
    x = SymOrb.build_path(P, Γ)

    ∇kin = SymOrb.∇kinetic(P, Γ)
    test_coefficients(∇kin)

    ∇U = SymOrb.∇U(P, x)
    @test ∇U isa Vector{Float64}
    @test size(∇U) == (P.N*P.dim*(P.steps+2), )

    ∇potential = SymOrb.∇potential(P, Γ)
    test_coefficients(∇potential)

    ∇action = SymOrb.∇action(P, Γ)
    test_coefficients(∇action)
end


@testset "Action hessian" begin

    Γ = SymOrb.random_starting_path(P)
    x = SymOrb.build_path(P, Γ)
    Hkin = SymOrb.Hkinetic(P, Γ)
    
    @test Hkin isa Matrix{Float64}
    @test size(Hkin) == ((P.N-1)*P.dim*(P.F+2), (P.N-1)*P.dim*(P.F+2))

    HU = SymOrb.HU(P, x)
    @test HU isa Matrix{Float64}
    @test size(HU) == (P.N*P.dim*(P.steps+2), P.N*P.dim*(P.steps+2))

    Hpotential = SymOrb.Hpotential(P, Γ)
    @test size(Hpotential) == ((P.N-1)*P.dim*(P.F+2), (P.N-1)*P.dim*(P.F+2))

    Haction = SymOrb.Haction(P, Γ)
    @test size(Hpotential) == ((P.N-1)*P.dim*(P.F+2), (P.N-1)*P.dim*(P.F+2))
end

end