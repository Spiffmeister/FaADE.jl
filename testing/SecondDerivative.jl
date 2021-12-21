# == TESTING FIRST DERIVATIVE ==#
# using Test

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

# Linear Function
## Constant coefficient
n, x, Δx = buildgrid(10)
c = ones(n)
u = x
∂ₓₓuₑ = zeros(n)
∂ₓₓu = Dₓₓ(u,n,Δx,c)

norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1]) ≤ 1.0e-14

## Constant coefficient
n, x, Δx = buildgrid(10)
c = x
u = x
∂ₓₓuₑ = ones(n)
∂ₓₓu = Dₓₓ(u,n,Δx,c)

norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1]) ≤ 1.0e-14


# Quadratic function
## Constant coefficient
n, x, Δx = buildgrid(10)
c = ones(n)
u = x.^2
∂ₓₓuₑ = 2*ones(n)
∂ₓₓu = Dₓₓ(u,n,Δx,c)

norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1]) ≤ 1.0e-14

## Variable coefficient

n, x, Δx = buildgrid(10)
c = x
u = x.^2
∂ₓₓuₑ = 4x
∂ₓₓu = Dₓₓ(u,n,Δx,c)

norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1]) ≤ 1.0e-14


# Cubic function
## Constant coefficient

n, x, Δx = buildgrid(10)
c = ones(n)
u = x.^3
∂ₓₓuₑ = 3x.^2
∂ₓₓu = Dₓₓ(u,n,Δx,c)

norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1]) ≤ 1.0e-14

## Variable coefficient

n, x, Δx = buildgrid(10)
c = x
u = x.^3
∂ₓₓuₑ = 9x.^2
∂ₓₓu = Dₓₓ(u,n,Δx,c)


##======##
# FOURTH ORDER
##======##

# Linear Function
## Constant coefficient
n, x, Δx = buildgrid(15)
c = ones(n)
u = x
∂ₓₓuₑ = zeros(n)
∂ₓₓu = Dₓₓ(u,n,Δx,c)

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

