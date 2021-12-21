# == TESTING FIRST DERIVATIVE ==#
using Test

using LinearAlgebra
using Plots
using Pkg

Pkg.activate(".")
using SBP_operators





function buildgrid(n)
    x = collect(range(0,stop=1,length=n+1))
    Δx = x[2]-x[1]
    return n+1, x, Δx
end






##======##
# SECOND ORDER
##======##

### CONSTANT COEFFICIENT
# Linear Function, ∂ₓ(1 ∂ₓx) = 0
n, x, Δx = buildgrid(10)
c = ones(n)
u = x
∂ₓₓuₑ = zeros(n)
∂ₓₓu = Dₓₓ(u,n,Δx,c)

@test norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1]) ≤ 1.0e-13


# Quadratic function, ∂ₓ(1 ∂ₓx²) = 2
n, x, Δx = buildgrid(10)
c = ones(n)
u = x.^2
∂ₓₓuₑ = 2*ones(n)
∂ₓₓu = Dₓₓ(u,n,Δx,c)

@test norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1]) ≤ 1.0e-13


# Cubic function, ∂ₓ(1 ∂ₓx³) = 3x²
n, x, Δx = buildgrid(10)
c = ones(n)
u = x.^3
∂ₓₓuₑ = 6x
∂ₓₓu = Dₓₓ(u,n,Δx,c)

@test norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1]) ≤ 1.0e-13


# Quartic function, ∂ₓ(1 ∂ₓx³) = 3x²
n, x, Δx = buildgrid(10)
c = ones(n)
u = x.^4
∂ₓₓuₑ = 12x.^2
∂ₓₓu = Dₓₓ(u,n,Δx,c)

lowres = norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1])

n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^4
∂ₓₓuₑ = 12x.^2
∂ₓₓu = Dₓₓ(u,n,Δx,c)

highres = norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1])


(lowres-highres)/highres


### VARIABLE COEFFICIENT

# Linear Funciton, ∂ₓ(x ∂ₓx) = 1
n, x, Δx = buildgrid(10)
c = x
u = x
∂ₓₓuₑ =  ones(n)
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







##======##
# FOURTH ORDER
##======##

# Linear Function
## Constant coefficient
n, x, Δx = buildgrid(15)
c = ones(n)
u = x
∂ₓₓuₑ = zeros(n)
∂ₓₓu = Dₓₓ(u,n,Δx,c,order=4)

norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1]) ≤ 1.0e-14

## Constant coefficient
n, x, Δx = buildgrid(15)
c = x
u = x
∂ₓₓuₑ = ones(n)
∂ₓₓu = Dₓₓ(u,n,Δx,c)

norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1]) ≤ 1.0e-14


# Quadratic function
## Constant coefficient
n, x, Δx = buildgrid(15)
c = ones(n)
u = x.^2
∂ₓₓuₑ = 2*ones(n)
∂ₓₓu = Dₓₓ(u,n,Δx,c)

norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1]) ≤ 1.0e-14

## Variable coefficient

n, x, Δx = buildgrid(15)
c = x
u = x.^2
∂ₓₓuₑ = 4x
∂ₓₓu = Dₓₓ(u,n,Δx,c)

norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1]) ≤ 1.0e-14


# Cubic function
## Constant coefficient

n, x, Δx = buildgrid(15)
c = ones(n)
u = x.^3
∂ₓₓuₑ = 3x.^2
∂ₓₓu = Dₓₓ(u,n,Δx,c)

norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1]) ≤ 1.0e-14

## Variable coefficient

n, x, Δx = buildgrid(15)
c = x
u = x.^3
∂ₓₓuₑ = 9x.^2
∂ₓₓu = Dₓₓ(u,n,Δx,c)


# Quartic function
## Constant coefficient

n, x, Δx = buildgrid(15)
c = ones(n)
u = x.^4
∂ₓₓuₑ = 4x.^3
∂ₓₓu = Dₓₓ(u,n,Δx,c)

norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1]) ≤ 1.0e-14

## Variable coefficient

n, x, Δx = buildgrid(15)
c = x
u = x.^4
∂ₓₓuₑ = 16x.^3
∂ₓₓu = Dₓₓ(u,n,Δx,c)

norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1]) ≤ 1.0e-14


# Quintic function
## Constant coefficient

n, x, Δx = buildgrid(15)
c = ones(n)
u = x.^5
∂ₓₓuₑ = 5x.^4
∂ₓₓu = Dₓₓ(u,n,Δx,c)

norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1]) ≤ 1.0e-14

## Variable coefficient

n, x, Δx = buildgrid(15)
c = x
u = x.^5
∂ₓₓuₑ = 25x.^4
∂ₓₓu = Dₓₓ(u,n,Δx,c)

norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1]) ≤ 1.0e-14


# Quintic function
## Constant coefficient

n, x, Δx = buildgrid(15)
c = ones(n)
u = x.^5
∂ₓₓuₑ = 5x.^4
∂ₓₓu = Dₓₓ(u,n,Δx,c)

norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1]) ≤ 1.0e-14

## Variable coefficient

n, x, Δx = buildgrid(15)
c = x
u = x.^5
∂ₓₓuₑ = 25x.^4
∂ₓₓu = Dₓₓ(u,n,Δx,c)

norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1]) ≤ 1.0e-14




##======##
# SIXTH ORDER
##======##

