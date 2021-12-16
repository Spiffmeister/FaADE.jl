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
    return n, x, Δx
end


n, x, Δx = buildgrid(10)

u = x.^2

∂ₓuₑ = 2x


∂ₓu = Dₓ!(u,n,Δx)



#=

#== 2nd order method ==#
# In the interior we should be able to exactly compute the soln to a quadratic
# On the boundary there should be some error since we have first order 
function second_order(u)
    # Compute the solution for ∂ₓu = 0 where u=x²

    
end

function second_order_soln()
    return x.^2
end



n = 100


u = second_order(x)
uₑ = second_order_soln(x)

# Test the interior
@test u[2:n-1] .- uₑ[2:n-1] .== 0

# Test the boundaries
@test u[1] - uₑ[1] < 1.0e-1
@test u[n] - uₑ[n] < 1.0e-1


#== 4th order method ==#
function forth_order()
    # Compute the solution for ∂ₓu = 0 where u=x⁴

    soln = x.^4

end

#== 6th order method ==# 
function sixth_order()
    # Compute the solution for ∂ₓu = 0 where u=x⁶

    soln = x.^6

end

=#