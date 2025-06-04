


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

"""
    1D TESTS
"""
n, x, Δx = buildgrid(1000)
u = ones(n);

k = ones(n);

∂ₓuₑ = zeros(n);
∂ₓu = zeros(n);

FaADE.Derivatives.FirstDerivativePeriodic!(∂ₓuₑ,k,u,Δx,Val(2),n,0.0)



"""

"""
n, x, Δx = buildgrid(1000)
u = sin.(12π*x)

k = ones(n);

∂ₓuₑ = zeros(n);
∂ₓu = zeros(n);

FaADE.Derivatives.FirstDerivativePeriodic!(∂ₓu,k,u,Δx,Val(2),n,0.0)




