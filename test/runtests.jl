using Test, SHA 
using MinPath, OffsetArrays

macro test_no_error(expr)
    quote
    try 
        res = $(esc(expr))
        @test true
        res
    catch 
        @test false
    end
    end
  end



# Test path.jl

@testset "All tests" verbose=true begin

configs = readdir("config", join=true)

function test_axes(M, ax)
    @test collect.(axes(M)) == collect.(ax)
end

function test_coefficients(Γ)
    @test Γ isa MinPath.Coefficients
    test_axes(Γ, (0:P.F+1,))
    test_configuration(Γ[0])
end


function test_hessian(M)
    @test M isa OffsetArrays.OffsetMatrix{Matrix{Matrix{Float64}}}
    test_axes(M, (0:P.F+1, 0:P.F+1))
    test_hessian_configuration(M[0,0])
end

function test_hessian_configuration(M)
    test_axes(M, (1:P.N, 1:P.N))
    test_axes(M[1,1], (1:P.dim, 1:P.dim))
end


function test_path(x)
    @test x isa MinPath.Path
    test_axes(x, (0:P.steps+1, ))
    test_configuration(x[0])
end

function test_configuration(c)
    test_axes(c, (1:P.N, ))
    test_axes(c[1], (1:P.dim, ))
end


@testset "Load configuration files" begin
    for file in configs
        P = @test_no_error initialize(file)
        @test P isa MinPath.Problem
    end
end


 P = initialize(configs[1])

@testset "Random starting path" begin
    Γ = @test_no_error MinPath.random_starting_path()
    test_coefficients(Γ)
end

@testset "Circular starting path" begin
    Γ = @test_no_error MinPath.circular_starting_path()
    test_coefficients(Γ)
end

@testset "Perturb path" begin
    Γ = MinPath.circular_starting_path()
    Γ = @test_no_error MinPath.perturbe_path(Γ)
    test_coefficients(Γ)
end

@testset "Perturbed circular starting path" begin
    Γ = @test_no_error MinPath.perturbed_circular_starting_path()
    test_coefficients(Γ)
end


@testset "Fourier series" begin
    Γ = MinPath.random_starting_path()
    s = @test_no_error MinPath.fourier_series(Γ)
    
    @test_no_error MinPath.inverse_fourier(s)
    segment = @test_no_error MinPath.segment(Γ[0], Γ[end])
    test_path(segment)
    
    @test Γ[0] ≈ segment[0]
    @test Γ[end] ≈ segment[end]

    x1 = @test_no_error MinPath.build_path(Γ)
    test_path(x1)

    Γ1 = @test_no_error MinPath.fourier_coefficients(x1)
    test_coefficients(Γ1)

    @test x1[0] ≈ Γ[0]
    @test x1[end] ≈ Γ[end]
    @test Γ ≈ Γ1 
end

@testset "Read and write path from file" begin
    tempdir = mktempdir()
    out_path = joinpath(tempdir, "test.txt")
    
    Γ = MinPath.circular_starting_path()
    @test_no_error MinPath.print_path_to_file(Γ, out_path)
    
    @test sha2_256(open(out_path)) == sha2_256(open("./out/circular_test_1.txt"))
    
    Γ1 = @test_no_error MinPath.read_path_from_file(out_path)
    @test Γ ≈ Γ1

    rm(tempdir, recursive=true)
end


@testset "Flatten and emboss paths" begin
    Γ = MinPath.circular_starting_path()
    v = @test_no_error MinPath.flatten(Γ)
    Γ1 = @test_no_error MinPath.emboss(v)
    test_coefficients(Γ1)
    @test Γ ≈ Γ1
end

@testset "Action" begin
    Γ = MinPath.circular_starting_path()

    x = MinPath.build_path(Γ)
    kin = @test_no_error MinPath.kinetic(Γ)
    @test kin isa Float64
    @test kin > 0
    U = @test_no_error MinPath.U(MinPath.build_path(Γ))
    @test U isa OffsetArrays.OffsetVector{Float64}
    test_axes(U, (0:P.steps+1, ))

    potential = @test_no_error MinPath.potential(Γ)
    @test potential isa Float64

    action = @test_no_error MinPath._action(Γ)
    @test action isa Float64
    @test action ≈ (kin + potential)

    action = @test_no_error MinPath.action(Γ)
    @test action isa Float64
    @test MinPath.action(MinPath.flatten(Γ)) ≈ action
end

@testset "Action gradient" begin
    Γ = MinPath.circular_starting_path()
    x = MinPath.build_path(Γ)
    ∇kin = @test_no_error MinPath.∇kinetic(Γ)
    @test ∇kin isa MinPath.Coefficients
    test_coefficients(∇kin)

    ∇U = @test_no_error MinPath.∇U(MinPath.build_path(Γ))
    @test ∇U isa OffsetArrays.OffsetVector{Vector{Vector{Float64}}}
    test_axes(∇U, (0:P.steps+1, ))

    ∇potential = @test_no_error MinPath.∇potential(Γ)
    test_coefficients(∇potential)

    ∇action = @test_no_error MinPath._∇action(Γ)
    test_coefficients(∇action)
    @test ∇action ≈ (∇kin + ∇potential)

    ∇action = @test_no_error MinPath.∇action(Γ)
    test_coefficients(∇action)
    @test MinPath.∇action(MinPath.flatten(Γ)) ≈ MinPath.flatten(∇action)
end


@testset "Action hessian" begin
    Γ = MinPath.circular_starting_path()
    x = MinPath.build_path(Γ)
    Hkin = @test_no_error MinPath.Hkinetic(Γ)
    test_hessian(Hkin)

    HU = @test_no_error MinPath.HU(MinPath.build_path(Γ))
    @test HU isa OffsetArrays.OffsetVector{Matrix{Matrix{Float64}}}
    test_axes(HU, (0:P.steps+1, ))
    test_hessian_configuration(HU[0])

    Hpotential = @test_no_error MinPath.Hpotential(Γ)
    test_hessian(Hpotential)

    Haction = @test_no_error MinPath._Haction(Γ)
    test_hessian(Haction)
    @test Haction ≈ (Hkin + Hpotential)

    Haction = @test_no_error MinPath.Haction(Γ)
    test_hessian(Haction)
    @test MinPath.Haction(MinPath.flatten(Γ)) ≈ MinPath.flatten(Haction)
end

end