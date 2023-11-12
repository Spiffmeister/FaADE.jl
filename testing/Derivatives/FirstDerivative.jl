#=
UNIT TESTING FIRST DERIVATIVE
=#

using Test

using LinearAlgebra
using Pkg

Pkg.activate(".")
using FaADE


function buildgrid(n)
    x = collect(range(0,stop=1,length=n+1))
    Δx = x[2]-x[1]
    return n+1, x, Δx
end



##======##
# SECOND ORDER
##======##
# In the interior we should be able to exactly compute the soln to a quadratic

# Linear function #1
@testset "Second order Dx" begin
    @testset "Linear" begin
        n, x, Δx = buildgrid(10)
        u = x
        ∂ₓuₑ = ones(n)
        ∂ₓu = D₁(u,Δx,order=2)

        @test norm(∂ₓuₑ[2:n-1] .- ∂ₓu[2:n-1]) ≤ 1.0e-14

        # Linear function #2
        n, x, Δx = buildgrid(10)
        u = 2x
        ∂ₓuₑ = 2ones(n)
        ∂ₓu = D₁(u,Δx,order=2)

        @test norm(∂ₓuₑ[2:n-1] .- ∂ₓu[2:n-1]) ≤ 1.0e-14
    end

    @testset "Quadratic" begin
        # Quadratic function
        n, x, Δx = buildgrid(10)
        u = x.^2
        ∂ₓuₑ = 2x
        ∂ₓu = D₁(u,Δx,order=2)
        
        @test norm(∂ₓuₑ[2:n-1] .- ∂ₓu[2:n-1]) ≤ 1.0e-14  
    end

    @testset "Cubic" begin
        # Cubic function
        n, x, Δx = buildgrid(10)
        u = x.^3
        ∂ₓuₑ = 3x.^2
        ∂ₓu = D₁(u,Δx,order=2)

        lowres = norm(∂ₓuₑ[2:n-1] .- ∂ₓu[2:n-1])

        n, x, Δx = buildgrid(20)
        u = x.^3
        ∂ₓuₑ = 3x.^2
        ∂ₓu = D₁(u,Δx,order=2)

        highres = norm(∂ₓuₑ[2:n-1] .- ∂ₓu[2:n-1])

        # @test lowres - highres 
    end
end






##======##
# FOURTH ORDER
##======##
# In the interior we should be able to exactly compute the soln to a cubic
# On the boundaries we should be able to compute to third order exactly

@testset "Fourth order Dx linear" begin

    @testset "Linear" begin
        # Linear function #1
        n, x, Δx = buildgrid(15)
        u = x
        ∂ₓuₑ = ones(n)
        ∂ₓu = D₁(u,Δx,order=4)
        @test norm(∂ₓuₑ[1:4] - ∂ₓu[1:4]) ≤ 1.0e-14      #Left boundary
        @test norm(∂ₓuₑ[5:n-4] - ∂ₓu[5:n-4]) ≤ 1.0e-14  #Interior
        @test norm(∂ₓuₑ[n-3:n] - ∂ₓu[n-3:n]) ≤ 1.0e-14  #Right boundary


        # Linear function #2
        n, x, Δx = buildgrid(10)
        u = 2x
        ∂ₓuₑ = 2ones(n)
        ∂ₓu = D₁(u,Δx,order=4)
        @test norm(∂ₓuₑ[2:n-1] .- ∂ₓu[2:n-1]) ≤ 1.0e-14
    end

    @testset "Quadratic" begin
        # Quadratic function
        n, x, Δx = buildgrid(15)
        u = x.^2
        ∂ₓuₑ = 2x
        ∂ₓu = D₁(u,Δx,order=4)
        @test norm(∂ₓuₑ[5:n-4] - ∂ₓu[5:n-4]) ≤ 1.0e-14
    end

    @testset "Cubic" begin
        # Cubic
        n, x, Δx = buildgrid(15)
        u = x.^3
        ∂ₓuₑ = 3x.^2
        ∂ₓu = D₁(u,Δx,order=4)
        @test norm(∂ₓuₑ[5:n-4] - ∂ₓu[5:n-4]) ≤ 1.0e-14
    end

    @testset "Quartic" begin
        # Quartic
        n, x, Δx = buildgrid(20)
        u = x.^4
        ∂ₓuₑ = 4x.^3
        ∂ₓu = D₁(u,Δx,order=4)

        @test norm(∂ₓuₑ[5:n-4] - ∂ₓu[5:n-4]) ≤ 1.0e-14
    end

    @testset "Quintic" begin
        # Quintic - This test should return the Test Broken expression as it fails to be under the tolerance
        n, x, Δx = buildgrid(15)
        u = x.^5
        ∂ₓuₑ = 5x.^4
        ∂ₓu = D₁(u,Δx,order=4)

        @test_broken norm(∂ₓuₑ[5:n-4] - ∂ₓu[5:n-4]) ≤ 1.0e-14
    end
end



#=
##======##
# SIXTH ORDER
##======##
# In the interior we should be able to exactly compute the soln to a hextic polynomial
# 

# Linear function
n, x, Δx = buildgrid(20)
u = x
∂ₓuₑ = ones(n)
∂ₓu = D₁(u,n,Δx,order=6)

@test norm(∂ₓuₑ[7:end-6] .- ∂ₓu[7:end-6]) ≤ 1.0e-14

# Quadratic
n, x, Δx = buildgrid(20)
u = x.^2
∂ₓuₑ = 2x
∂ₓu = D₁(u,n,Δx,order=6)

@test norm(∂ₓuₑ[7:end-6] .- ∂ₓu[7:n-6]) ≤ 1.0e-14

# Cubic
n, x, Δx = buildgrid(20)
u = x.^3
∂ₓuₑ = 3x.^2
∂ₓu = D₁(u,n,Δx,order=6)

@test norm(∂ₓuₑ[7:end-6] .- ∂ₓu[7:n-6]) ≤ 1.0e-14

# Quartic
n, x, Δx = buildgrid(20)
u = x.^4
∂ₓuₑ = 4x.^3
∂ₓu = D₁(u,n,Δx,order=6)

@test norm(∂ₓuₑ[7:end-6] .- ∂ₓu[7:n-6]) ≤ 1.0e-14

# Quintic
n, x, Δx = buildgrid(20)
u = x.^5
∂ₓuₑ = 5x.^4
∂ₓu = D₁(u,n,Δx,order=6)

@test norm(∂ₓuₑ[7:n-6] - ∂ₓu[7:n-6]) ≤ 1.0e-14

# Quintic - This test should return the Test Broken expression as it fails to be under the tolerance
n, x, Δx = buildgrid(20)
u = x.^6
∂ₓuₑ = 6x.^5
∂ₓu = D₁(u,n,Δx,order=6)

@test norm(∂ₓuₑ[7:n-6] - ∂ₓu[7:n-6]) ≤ 1.0e-14

# Hextic - This test should return the Test Broken expression as it fails to be under the tolerance
n, x, Δx = buildgrid(20)
u = x.^7
∂ₓuₑ = 7x.^6
∂ₓu = D₁(u,n,Δx,order=6)

@test_broken norm(∂ₓuₑ[7:n-6] - ∂ₓu[7:n-6]) ≤ 1.0e-14
=#



