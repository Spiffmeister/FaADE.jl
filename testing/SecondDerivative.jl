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

# Quadratic function
## Constant coefficient
n, x, Δx = buildgrid(10)
c = ones(n)
u = x.^2
∂ₓₓuₑ = 2*ones(n)
∂ₓₓu = Dₓₓ(u,n,Δx,c)

## Variable coefficient

n, x, Δx = buildgrid(10)
c = x
u = x.^2
∂ₓₓuₑ = 4x
∂ₓₓu = Dₓₓ(u,n,Δx,c)

# Cubic function
## Constant coefficient

n, x, Δx = buildgrid(10)
c = ones(n)
u = x.^3
∂ₓₓuₑ = 4x
∂ₓₓu = Dₓₓ(u,n,Δx,c)

## Constant coefficient

n, x, Δx = buildgrid(10)
c = x
u = x.^3
∂ₓₓuₑ = 4x
∂ₓₓu = Dₓₓ(u,n,Δx,c)



#= VARIABLE COEFFICIENT =#



n, x, Δx = buildgrid(10)
c = x
u = x.^2
∂ₓₓuₑ = 4x
∂ₓₓu = Dₓₓ(u,n,Δx,c)
