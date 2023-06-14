# == TESTING FIRST DERIVATIVE ==#
using Test

using LinearAlgebra
using Plots
using Pkg

Pkg.activate(".")
using SPADE





function buildgrid(n)
    x = collect(range(0.0,stop=1.0,length=n+1))
    Δx = x[2]-x[1]
    return n+1, x, Δx
end


##======##
# SECOND ORDER
order = 2
##======##
@testset "1D second derivative second order" begin
    n = 21
    Dom = SPADE.Grid1D([-1.0,1.0],n)
    x = Dom.grid
    @testset "constant coefficient" begin
        k = ones(n)
    
        # Linear Function u=x, ∂ₓ(1 ∂ₓx) = 0
        u = x
        ∂ₓₓu = D2(u,k,Dom.Δ,order=order)
        @test all(zeros(n) .== ∂ₓₓu) atol=1e-14
    
        # Quadratic function, ∂ₓ(1 ∂ₓx²) = 2
        u = x.^2
        ∂ₓₓuₑ = 2*ones(n)
        ∂ₓₓu = D₂(u,k,Dom.Δ,order=order)
        @test norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1]) ≤ 1.0e-13
    
        # Cubic function, ∂ₓ(1 ∂ₓx³) = 6x
        u = x.^3
        ∂ₓₓuₑ = 6x
        ∂ₓₓu = D₂(u,k,Dom.Δ,order=order)
        @test norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1]) ≤ 1.0e-13
    end

    @testset "variable coefficient" begin
        c = x

        u = x
        ∂ₓₓuₑ = ones(n)
        ∂ₓₓu = D₂(u,k,n,Dom.Δ)
        @test norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1]) ≤ 1.0e-13
    
        # Quadratic Function, ∂ₓ(x ∂ₓx²) = 4x
        u = x.^2
        ∂ₓₓuₑ = 4x
        ∂ₓₓu = D₂(u,n,Δx,c)
        @test norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1]) ≤ 1.0e-13
    
        # Cubic Function, ∂ₓ(x ∂ₓx³) = 9x²
        u = x.^3
        ∂ₓₓuₑ = 9x.^2
        ∂ₓₓu = D₂(u,n,Δx,c)
        @test_broken norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1]) ≤ 1.0e-13
    end
end

@testset "2D second derivative second order" begin
    nx = ny = 21
end





# # Quartic function, ∂ₓ(1 ∂ₓx⁴) = 12x²
# u = x.^4
# ∂ₓₓuₑ = 12x.^2
# ∂ₓₓu = D₂(u,n,Δx,c)

# lowres = norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1])

# n, x, Δx = buildgrid(20)
# c = ones(n)
# u = x.^4
# ∂ₓₓuₑ = 12x.^2
# ∂ₓₓu = D₂(u,n,Δx,c)

# highres = norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1])

# (lowres-highres)/highres



#= 2D methods =#
@testset "Dxx second order 2D" begin
    nx = ny = 21
    Dx = Dy = [-1.0,1.0]
    Dom = SPADE.Grid2D(Dx,Dy,nx,ny)


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
∂ₓₓu = D₂(u,n,Δx,c,order=4)
@test norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-13

# Quadratic function
n, x, Δx = buildgrid(15)
c = ones(n)
u = x.^2
∂ₓₓuₑ = 2*ones(n)
∂ₓₓu = D₂(u,n,Δx,c,order=4)
@test norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-12

# Cubic function
n, x, Δx = buildgrid(15)
c = ones(n)
u = x.^3
∂ₓₓuₑ = 6x
∂ₓₓu = D₂(u,n,Δx,c,order=4)
@test norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-13

# Quartic function
n, x, Δx = buildgrid(15)
c = ones(n)
u = x.^4
∂ₓₓuₑ = 12x.^2
∂ₓₓu = D₂(u,n,Δx,c,order=4)
@test norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-13

# Quintic function
n, x, Δx = buildgrid(15)
c = ones(n)
u = x.^5
∂ₓₓuₑ = 20x.^3
∂ₓₓu = D₂(u,n,Δx,c,order=4)
@test norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-13

# Hextic function
n, x, Δx = buildgrid(15)
c = ones(n)
u = x.^6
∂ₓₓuₑ = 30x.^4
∂ₓₓu = D₂(u,n,Δx,c,order=4)
@test_broken norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-12


### Variable coefficient
# Linear function
n, x, Δx = buildgrid(15)
c = x
u = x
∂ₓₓuₑ = ones(n)
∂ₓₓu = D₂(u,n,Δx,c,order=4)
@test norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-13

# Quadratic function
n, x, Δx = buildgrid(15)
c = x
u = x.^2
∂ₓₓuₑ = 4x
∂ₓₓu = D₂(u,n,Δx,c,order=4)
@test norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-13

# Cubic function
n, x, Δx = buildgrid(15)
c = x
u = x.^3
∂ₓₓuₑ = 9x.^2
∂ₓₓu = D₂(u,n,Δx,c,order=4)
@test norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-13

# BROKEN TESTS
# Quartic function
n, x, Δx = buildgrid(15)
c = x
u = x.^4
∂ₓₓuₑ = 16x.^3
∂ₓₓu = D₂(u,n,Δx,c)
@test_broken norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-14


# Quadratic coefficient
n, x, Δx = buildgrid(15)
c = x.^2
u = x
∂ₓₓuₑ = ones(n)
∂ₓₓu = D₂(u,n,Δx,c,order=4)
@test_broken norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-13

#=
##======##
# SIXTH ORDER
##======##

### Constant coefficient
# Linear Function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x
∂ₓₓuₑ = zeros(n)
∂ₓₓu = D₂(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ[6:n-5] .- ∂ₓₓu[6:n-5]) ≤ 1.0e-12

# Quadratic function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^2
∂ₓₓuₑ = 2*ones(n)
∂ₓₓu = D₂(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ .- ∂ₓₓu) ≤ 1.0e-10

# Cubic function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^3
∂ₓₓuₑ = 6x
∂ₓₓu = D₂(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ .- ∂ₓₓu) ≤ 1.0e-10

# Quartic function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^4
∂ₓₓuₑ = 12x.^2
∂ₓₓu = D₂(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ .- ∂ₓₓu) ≤ 1.0e-12

# Quintic function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^5
∂ₓₓuₑ = 20x.^3
∂ₓₓu = D₂(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ .- ∂ₓₓu) ≤ 1.0e-12

# Hextic function - only body is 6th order
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^6
∂ₓₓuₑ = 30x.^4
∂ₓₓu = D₂(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ .- ∂ₓₓu) ≤ 1.0e-13

# Heptic coefficient
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^7
∂ₓₓuₑ = 42x.^5
∂ₓₓu = D₂(u,n,Δx,c,order=6)

@test_broken norm(∂ₓₓuₑ .- ∂ₓₓu) ≤ 1.0e-13



### VARIABLE COEFFICIENT
# Linear Function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x
∂ₓₓuₑ = zeros(n)
∂ₓₓu = D₂(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ[6:n-5] .- ∂ₓₓu[6:n-5]) ≤ 1.0e-12

# Quadratic function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^2
∂ₓₓuₑ = 2*ones(n)
∂ₓₓu = D₂(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ[6:n-5] .- ∂ₓₓu[6:n-5]) ≤ 1.0e-12

# Cubic function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^3
∂ₓₓuₑ = 6x
∂ₓₓu = D₂(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ[6:n-5] .- ∂ₓₓu[6:n-5]) ≤ 1.0e-12

# Quartic function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^4
∂ₓₓuₑ = 12x.^2
∂ₓₓu = D₂(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ[6:n-5] .- ∂ₓₓu[6:n-5]) ≤ 1.0e-12

# Quintic function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^5
∂ₓₓuₑ = 20x.^3
∂ₓₓu = D₂(u,n,Δx,c,order=6)

@test norm(∂ₓₓuₑ[7:n-6] .- ∂ₓₓu[7:n-6]) ≤ 1.0e-12

# Hextic function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^6
∂ₓₓuₑ = 30x.^4
∂ₓₓu = D₂(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ[7:n-6] .- ∂ₓₓu[7:n-6]) ≤ 1.0e-12

# Heptic coefficient
n, x, Δx = buildgrid(20)
c = x
u = x.^7
∂ₓₓuₑ = 42x.^5
∂ₓₓu = D₂(u,n,Δx,c,order=6)
@test_broken norm(∂ₓₓuₑ[7:n-6] .- ∂ₓₓu[7:n-6]) ≤ 1.0e-14
=#