using LinearAlgebra
using Printf
# using GLMakie
# pyplot()

using BenchmarkTools
using Profile

using Pkg
Pkg.activate(".")
using SPADE




###
𝒟x = [0.0,1.0]
𝒟y = [0.0,1.0]
nx = 41
ny = 31
Dom = Grid2D(𝒟x,𝒟y,nx,ny)

kx = zeros(Float64,nx,ny) .+ 1.0;
ky = zeros(Float64,nx,ny) .+ 1.0;

Δt = 0.01* min(Dom.Δx^2,Dom.Δy^2)
t_f = 100Δt

u₀(x,y) = exp(-((x-0.5)^2 + (y-0.5)^2) / 0.02)


BoundaryLeftRight = PeriodicBoundary(1)
BoundaryUpDown = PeriodicBoundary(2)


order = 2
method = :cgie

P = VariableCoefficientPDE2D(u₀,kx,ky,order,BoundaryLeftRight,BoundaryUpDown)

println(method)
println("Δx=",Dom.Δx,"      ","Δt=",Δt,"        ","final time=",t_f)


###
@benchmark solve($P,$Dom,$Δt,$t_f,$method)



# soln = solve(P,Dom,Δt,t_f,:cgie)
# using GLMakie
# using Plots
# surface(soln.grid.gridy,soln.grid.gridx,soln.u[2],
#     xlabel="y",ylabel="x",zlabel="Temp",
#     xlims=(0.0,1.0), ylims=(0.0,1.0), zlims=(0.0,1.0))



# @time solve(P,Dom,Δt,t_f,:cgie)
# Profile.clear_malloc_data()
# @time solve(P,Dom,Δt,t_f,:cgie)