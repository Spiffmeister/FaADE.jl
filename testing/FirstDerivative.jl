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




#=

#== 2nd order method ==#
# In the interior we should be able to exactly compute the soln to a quadratic
# On the boundary there should be some error since we have first order

# Linear function #1

n, x, Δx = buildgrid(10)
u = x
∂ₓuₑ = ones(n)
∂ₓu = Dₓ(u,n,Δx)


# Linear function #2

n, x, Δx = buildgrid(10)
u = 2x
∂ₓuₑ = 2ones(n)
∂ₓu = Dₓ(u,n,Δx)


# Quadratic function

n, x, Δx = buildgrid(20)
u = x.^2
∂ₓuₑ = 2x
∂ₓu = Dₓ(u,n,Δx)


=#

#=
#== 4th order method ==#

# Linear function #1
n, x, Δx = buildgrid(15)
u = x
∂ₓuₑ = ones(n)
∂ₓu = Dₓ(u,n,Δx,order=4)

@test norm(∂ₓuₑ[5:n-4] - ∂ₓu[5:n-3]) ≤ 1.0e-14


# Quadratic function
n, x, Δx = buildgrid(15)
u = x.^2
∂ₓuₑ = 2x
∂ₓu = Dₓ(u,n,Δx,order=4)

@test norm(∂ₓuₑ[5:n-4] - ∂ₓu[5:n-3]) ≤ 1.0e-14




# Cubic
n, x, Δx = buildgrid(15)
u = x.^3
∂ₓuₑ = 3x.^2
∂ₓu = Dₓ(u,n,Δx,order=4)

@test norm(∂ₓuₑ[5:n-4] - ∂ₓu[5:n-4]) ≤ 1.0e-14



# Quartic
n, x, Δx = buildgrid(20)
u = x.^4
∂ₓuₑ = 4x.^3
∂ₓu = Dₓ(u,n,Δx,order=4)

@test norm(∂ₓuₑ[5:n-4] - ∂ₓu[5:n-4]) ≤ 1.0e-14



# Quintic - This test should return the Test Broken expression as it fails to be under the tolerance
n, x, Δx = buildgrid(15)
u = x.^5
∂ₓuₑ = 5x.^4
∂ₓu = Dₓ(u,n,Δx,order=4)

@test_broken norm(∂ₓuₑ[5:n-4] - ∂ₓu[5:n-4]) ≤ 1.0e-14

=#



#== 6th order method ==# 


# Linear function #1

n, x, Δx = buildgrid(20)
u = x
∂ₓuₑ = ones(n)
∂ₓu = Dₓ(u,n,Δx,order=6)

∂ₓuₑ[7:end-6] .≈ ∂ₓu[7:n-6]


# Quadratic function

n, x, Δx = buildgrid(20)
u = x.^2
∂ₓuₑ = 2x
∂ₓu = Dₓ(u,n,Δx)

∂ₓuₑ[7:end-6] .≈ ∂ₓu[7:n-6]


# Cubic

n, x, Δx = buildgrid(15)
u = x.^3
∂ₓuₑ = 3x.^2
∂ₓu = Dₓ(u,n,Δx)

∂ₓuₑ[7:end-6] .≈ ∂ₓu[7:n-6]


# Quartic

n, x, Δx = buildgrid(15)
u = x.^4
∂ₓuₑ = 4x.^3
∂ₓu = Dₓ(u,n,Δx)

∂ₓuₑ[7:end-6] .≈ ∂ₓu[7:n-6]


# Quintic
n, x, Δx = buildgrid(15)
u = x.^5
∂ₓuₑ = 5x.^4
∂ₓu = Dₓ(u,n,Δx,order=4)

@test norm(∂ₓuₑ[7:n-6] - ∂ₓu[7:n-6]) ≤ 1.0e-14


# Quintic - This test should return the Test Broken expression as it fails to be under the tolerance
n, x, Δx = buildgrid(15)
u = x.^5
∂ₓuₑ = 5x.^4
∂ₓu = Dₓ(u,n,Δx,order=4)

@test_broken norm(∂ₓuₑ[7:n-6] - ∂ₓu[7:n-6]) ≤ 1.0e-14


# Quintic - This test should return the Test Broken expression as it fails to be under the tolerance
n, x, Δx = buildgrid(15)
u = x.^5
∂ₓuₑ = 5x.^4
∂ₓu = Dₓ(u,n,Δx,order=4)

@test_broken norm(∂ₓuₑ[7:n-6] - ∂ₓu[7:n-6]) ≤ 1.0e-14

