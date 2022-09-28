using LinearAlgebra
using Printf
using Plots
# pyplot()

using BenchmarkTools

using Pkg
Pkg.activate(".")
using SBP_operators




###
𝒟x = [0.0,1.0]
𝒟y = [0.0,1.0]
nx = 41
ny = 41
Dom = Grid2D(𝒟x,𝒟y,nx,ny)

kx = zeros(Float64,nx,ny) .+ 1.0;
ky = zeros(Float64,nx,ny) .+ 1.0;

Δt = 0.01* min(Δx^2,Δy^2)
t_f = 2.1Δt

u₀(x,y) = exp(-((x-0.5)^2 + (y-0.5)^2) / 0.02)


BoundaryLeft = Boundary(Dirichlet,t -> 0.0,Left,1)
BoundaryRight = Boundary(Dirichlet,t -> 0.0,Right,1)
BoundaryUp = Boundary(Diricihlet,t -> 0.0,Up,2)
BoundaryDown = Boundary(Dirichlet,t -> 0.0,Down,2)

order = 2
method = :cgie

P = VariableCoefficientPDE2D(u₀,kx,ky,order,BoundaryLeft,BoundaryRight,BoundaryUp,BoundaryDown)

println(method)
println("Δx=",Δx,"      ","Δt=",Δt,"        ","final time=",t_f)


###
@benchmark solve(P,Dom,Δt,t_f,method)



