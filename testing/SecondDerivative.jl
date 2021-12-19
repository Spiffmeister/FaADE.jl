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


n, x, Δx = buildgrid(10)




##

c = ones(n)

u = x.^2
∂ₓₓuₑ = 2*ones(n)

∂ₓₓu = Dₓₓ(u,n,Δx,c)


##

n, x, Δx = buildgrid(10)


c = x
u = x.^2

# c = ones(n)
# u = x.^3

∂ₓₓuₑ = 6x


∂ₓₓu = Dₓₓ(u,n,Δx,c)

∂ₓₓuₑ .- ∂ₓₓu
