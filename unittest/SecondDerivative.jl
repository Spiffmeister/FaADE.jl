# == TESTING FIRST DERIVATIVE ==#
using Test

using LinearAlgebra
using Plots
using Pkg


using SBP_operators





function buildgrid(n)
    x = collect(range(0.0,stop=1.0,length=n+1))
    Δx = x[2]-x[1]
    return n+1, x, Δx
end






##======##
# SECOND ORDER
order = 2
##======##
@testset "Second order constant coefficient" begin
    n, x, Δx = buildgrid(10)
    c = ones(n)

    # Linear Function, ∂ₓ(1 ∂ₓx) = 0
    u = x
    ∂ₓₓuₑ = zeros(n)
    ∂ₓₓu = Dₓₓ(u,c,n,Δx,order=order)
    @test all(∂ₓₓuₑ .== ∂ₓₓu) atol=1e-14

    # Quadratic function, ∂ₓ(1 ∂ₓx²) = 2
    u = x.^2
    ∂ₓₓuₑ = 2*ones(n)
    ∂ₓₓu = Dₓₓ(u,n,Δx,c)
    @test norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1]) ≤ 1.0e-13

    # Cubic function, ∂ₓ(1 ∂ₓx³) = 3x²
    u = x.^3
    ∂ₓₓuₑ = 6x
    ∂ₓₓu = Dₓₓ(u,n,Δx,c)
    @test norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1]) ≤ 1.0e-13
end

# # Quartic function, ∂ₓ(1 ∂ₓx⁴) = 12x²
# u = x.^4
# ∂ₓₓuₑ = 12x.^2
# ∂ₓₓu = Dₓₓ(u,n,Δx,c)

# lowres = norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1])

# n, x, Δx = buildgrid(20)
# c = ones(n)
# u = x.^4
# ∂ₓₓuₑ = 12x.^2
# ∂ₓₓu = Dₓₓ(u,n,Δx,c)

# highres = norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1])

# (lowres-highres)/highres


@testset "Second order variable coefficient" begin
    # Linear Funciton, ∂ₓ(x ∂ₓx) = 1
    n, x, Δx = buildgrid(10)
    c = x
    u = x
    ∂ₓₓuₑ = ones(n)
    ∂ₓₓu = Dₓₓ(u,n,Δx,c)
    @test norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1]) ≤ 1.0e-13

    # Quadratic Function, ∂ₓ(x ∂ₓx²) = 4x
    n, x, Δx = buildgrid(10)
    c = x
    u = x.^2
    ∂ₓₓuₑ = 4x
    ∂ₓₓu = Dₓₓ(u,n,Δx,c)
    @test norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1]) ≤ 1.0e-13

    # Cubic Function, ∂ₓ(x ∂ₓx³) = 9x²
    n, x, Δx = buildgrid(10)
    c = x
    u = x.^3
    ∂ₓₓuₑ = 9x.^2
    ∂ₓₓu = Dₓₓ(u,n,Δx,c)
    @test_broken norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1]) ≤ 1.0e-13
end





##======##
# FOURTH ORDER
##======##

### Constant coefficient
# Linear Function
n, x, Δx = buildgrid(15)
c = ones(n)
u = x
∂ₓₓuₑ = zeros(n)
∂ₓₓu = Dₓₓ(u,n,Δx,c,order=4)
@test norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-13

# Quadratic function
n, x, Δx = buildgrid(15)
c = ones(n)
u = x.^2
∂ₓₓuₑ = 2*ones(n)
∂ₓₓu = Dₓₓ(u,n,Δx,c,order=4)
@test norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-12

# Cubic function
n, x, Δx = buildgrid(15)
c = ones(n)
u = x.^3
∂ₓₓuₑ = 6x
∂ₓₓu = Dₓₓ(u,n,Δx,c,order=4)
@test norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-13

# Quartic function
n, x, Δx = buildgrid(15)
c = ones(n)
u = x.^4
∂ₓₓuₑ = 12x.^2
∂ₓₓu = Dₓₓ(u,n,Δx,c,order=4)
@test norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-13

# Quintic function
n, x, Δx = buildgrid(15)
c = ones(n)
u = x.^5
∂ₓₓuₑ = 20x.^3
∂ₓₓu = Dₓₓ(u,n,Δx,c,order=4)
@test norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-13

# Hextic function
n, x, Δx = buildgrid(15)
c = ones(n)
u = x.^6
∂ₓₓuₑ = 30x.^4
∂ₓₓu = Dₓₓ(u,n,Δx,c,order=4)
@test_broken norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-12


### Variable coefficient
# Linear function
n, x, Δx = buildgrid(15)
c = x
u = x
∂ₓₓuₑ = ones(n)
∂ₓₓu = Dₓₓ(u,n,Δx,c,order=4)
@test norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-13

# Quadratic function
n, x, Δx = buildgrid(15)
c = x
u = x.^2
∂ₓₓuₑ = 4x
∂ₓₓu = Dₓₓ(u,n,Δx,c,order=4)
@test norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-13

# Cubic function
n, x, Δx = buildgrid(15)
c = x
u = x.^3
∂ₓₓuₑ = 9x.^2
∂ₓₓu = Dₓₓ(u,n,Δx,c,order=4)
@test norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-13

# BROKEN TESTS
# Quartic function
n, x, Δx = buildgrid(15)
c = x
u = x.^4
∂ₓₓuₑ = 16x.^3
∂ₓₓu = Dₓₓ(u,n,Δx,c)
@test_broken norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-14


# Quadratic coefficient
n, x, Δx = buildgrid(15)
c = x.^2
u = x
∂ₓₓuₑ = ones(n)
∂ₓₓu = Dₓₓ(u,n,Δx,c,order=4)
@test_broken norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-13


##======##
# SIXTH ORDER
##======##

### Constant coefficient
# Linear Function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x
∂ₓₓuₑ = zeros(n)
∂ₓₓu = Dₓₓ(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ[6:n-5] .- ∂ₓₓu[6:n-5]) ≤ 1.0e-12

# Quadratic function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^2
∂ₓₓuₑ = 2*ones(n)
∂ₓₓu = Dₓₓ(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ .- ∂ₓₓu) ≤ 1.0e-10

# Cubic function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^3
∂ₓₓuₑ = 6x
∂ₓₓu = Dₓₓ(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ .- ∂ₓₓu) ≤ 1.0e-10

# Quartic function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^4
∂ₓₓuₑ = 12x.^2
∂ₓₓu = Dₓₓ(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ .- ∂ₓₓu) ≤ 1.0e-12

# Quintic function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^5
∂ₓₓuₑ = 20x.^3
∂ₓₓu = Dₓₓ(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ .- ∂ₓₓu) ≤ 1.0e-12

# Hextic function - only body is 6th order
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^6
∂ₓₓuₑ = 30x.^4
∂ₓₓu = Dₓₓ(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ .- ∂ₓₓu) ≤ 1.0e-13

# Heptic coefficient
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^7
∂ₓₓuₑ = 42x.^5
∂ₓₓu = Dₓₓ(u,n,Δx,c,order=6)

@test_broken norm(∂ₓₓuₑ .- ∂ₓₓu) ≤ 1.0e-13



### VARIABLE COEFFICIENT
# Linear Function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x
∂ₓₓuₑ = zeros(n)
∂ₓₓu = Dₓₓ(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ[6:n-5] .- ∂ₓₓu[6:n-5]) ≤ 1.0e-12

# Quadratic function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^2
∂ₓₓuₑ = 2*ones(n)
∂ₓₓu = Dₓₓ(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ[6:n-5] .- ∂ₓₓu[6:n-5]) ≤ 1.0e-12

# Cubic function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^3
∂ₓₓuₑ = 6x
∂ₓₓu = Dₓₓ(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ[6:n-5] .- ∂ₓₓu[6:n-5]) ≤ 1.0e-12

# Quartic function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^4
∂ₓₓuₑ = 12x.^2
∂ₓₓu = Dₓₓ(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ[6:n-5] .- ∂ₓₓu[6:n-5]) ≤ 1.0e-12

# Quintic function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^5
∂ₓₓuₑ = 20x.^3
∂ₓₓu = Dₓₓ(u,n,Δx,c,order=6)

@test norm(∂ₓₓuₑ[7:n-6] .- ∂ₓₓu[7:n-6]) ≤ 1.0e-12

# Hextic function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^6
∂ₓₓuₑ = 30x.^4
∂ₓₓu = Dₓₓ(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ[7:n-6] .- ∂ₓₓu[7:n-6]) ≤ 1.0e-12

# Heptic coefficient
n, x, Δx = buildgrid(20)
c = x
u = x.^7
∂ₓₓuₑ = 42x.^5
∂ₓₓu = Dₓₓ(u,n,Δx,c,order=6)
@test_broken norm(∂ₓₓuₑ[7:n-6] .- ∂ₓₓu[7:n-6]) ≤ 1.0e-14
