using Test, SHA 
using SymPath


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
        @test P isa SymPath.Problem
        @test P.Π^2 ≈ P.Π atol=1e-3
    end
end


P = initialize(configs[1])
dims = (F = P.F, N = P.N, dim = P.dim)

@testset "Random starting path" begin
    Γ = SymPath.random_starting_path(P)
    test_coefficients(Γ)
end

@testset "Circular starting path" begin
    Γ = SymPath.circular_starting_path(P)
    test_coefficients(Γ)
end

@testset "Perturbe path" begin
    Γ = SymPath.circular_starting_path(P)
    Γ = SymPath.perturbe_path(P, Γ)
    test_coefficients(Γ)
end

@testset "Perturbed circular starting path" begin
    Γ = SymPath.perturbed_circular_starting_path(P)
    test_coefficients(Γ)
end

@testset "Fourier series" begin
  
    Γ = SymPath.random_starting_path(P)

    x1 = SymPath.build_path(P, Γ, 500)
    @test x1 isa Array{Float64, 3}

    Γ1 = SymPath.fourier_coefficients(P, x1, 500)
    test_coefficients(Γ1)

    @test Γ ≈ Γ1  atol = 1e-3
end

@testset "Read and write path from file" begin
    tempdir = mktempdir()
    out_path = joinpath(tempdir, "test.txt")
    
    Γ = SymPath.circular_starting_path(P)
    SymPath.print_path_to_file(P, Γ, out_path)
    
    Γ1 = SymPath.read_path_from_file(P, out_path)

    @test Γ ≈ Γ1  atol = 1e-3
    rm(tempdir, recursive=true)
end


@testset "Action" begin
    
    Γ = SymPath.random_starting_path(P)

    x = SymPath.build_path(P, Γ)
    kin = SymPath.kinetic(P, Γ)
    @test kin isa Float64
    @test kin > 0
    U = SymPath.U(P, x)
    @test U isa Vector{Float64}
    @test size(U) == (P.steps+2, )

    potential = SymPath.potential(P, Γ)
    @test potential isa Float64

    action = SymPath.action(P, Γ)
    @test action isa Float64
end

@testset "Action gradient" begin

    Γ = SymPath.random_starting_path(P)
    x = SymPath.build_path(P, Γ)

    ∇kin = SymPath.∇kinetic(P, Γ)
    test_coefficients(∇kin)

    ∇U = SymPath.∇U(P, x)
    @test ∇U isa Vector{Array{Float64, 2}}
    @test size(∇U) == (P.steps+2, )
    @test size(∇U[1]) == (P.dim, P.N)

    ∇potential = SymPath.∇potential(P, Γ)
    test_coefficients(∇potential)

    ∇action = SymPath.∇action(P, Γ)
    test_coefficients(∇action)
end


@testset "Action hessian" begin

    Γ = SymPath.random_starting_path(P)
    x = SymPath.build_path(P, Γ)
    Hkin = SymPath.Hkinetic(P, Γ)
    
    @test Hkin isa Matrix{Float64}
    @test size(Hkin) == ((P.N-1)*P.dim*(P.F+2), (P.N-1)*P.dim*(P.F+2))

    HU = SymPath.HU(P, x)
    @test HU isa Vector{Array{Float64, 4}}
    @test size(HU) == (P.steps+2, )
    @test size(HU[1]) == (P.dim, P.N, P.dim, P.N)

    Hpotential = SymPath.Hpotential(P, Γ)
    @test size(Hpotential) == ((P.N-1)*P.dim*(P.F+2), (P.N-1)*P.dim*(P.F+2))

    Haction = SymPath.Haction(P, Γ)
    @test size(Hpotential) == ((P.N-1)*P.dim*(P.F+2), (P.N-1)*P.dim*(P.F+2))
end

end