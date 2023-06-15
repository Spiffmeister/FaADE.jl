using LinearAlgebra
using Printf
using Plots
# pyplot()
# using GLMakie
using Interpolations
using JLD2

# using BenchmarkTools

using Distributed
# addprocs(1)
@everywhere push!(LOAD_PATH,".")
@everywhere using FaADE
using SharedArrays

###
@everywhere 𝒟x = [0.5,0.68]
@everywhere 𝒟y = [-π,π]
@everywhere nx = 41
@everywhere ny = 41

@everywhere Δx = (𝒟x[2]-𝒟x[1])/(nx-1)
@everywhere Δy = (𝒟y[2]-𝒟y[1])/(ny-1)
@everywhere x = collect(range(𝒟x[1],𝒟x[2],step=Δx))
@everywhere y = collect(range(𝒟y[1],𝒟y[2],step=Δy))

@everywhere kx = zeros(Float64,nx,ny) .+ 1.0
@everywhere ky = zeros(Float64,nx,ny) .+ 1.0

# u₀(x,y) = exp(-(x-0.5)^2/0.02 - (y-π)^2/0.5)
# u₀(x,y) = 0.5sin(4*2π*x) + 0.5sin(4*y)
# @everywhere u₀(x,y) = (x-0.5)/(0.68-0.5)
@everywhere u₀(x,y) = (x-0.5)^2/(0.68-0.5)

u = SharedArray(zeros(nx,ny))
@sync @distributed for i = 1:nx
    for j = 1:ny
        u[i,j] = u₀(x[i],y[j])
    end
end
uₓₓ = D₂(u,nx,ny,Δx,Δy,kx,ky,order_x=4)

uxx = zeros(nx,ny)
for i = 1:nx
    for j = 1:ny
        uxx[i,j] = 2/(0.68-0.5)
    end
end
