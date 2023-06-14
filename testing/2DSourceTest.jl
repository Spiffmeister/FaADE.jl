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
nx = 21
ny = 21
Dom = Grid2D(𝒟x,𝒟y,nx,ny)

kx = zeros(Float64,nx,ny) .+ 1.0;
ky = zeros(Float64,nx,ny) .+ 1.0;

Δt = 1.0* min(Dom.Δx^2,Dom.Δy^2)
t_f = 100Δt

# u₀(x,y) = exp(-((x-0.5)^2 + (y-0.5)^2) / 0.02)
u₀(x,y) = x



BoundaryLeft = Boundary(Dirichlet,(x,t) -> 0.0,Left,1)
BoundaryRight = Boundary(Dirichlet,(x,t) -> 1.0,Right,1)
BoundaryUp = Boundary(Neumann,(y,t) -> 0.0,Up,2)
BoundaryDown = Boundary(Neumann,(y,t) -> 0.0,Down,2)


# F(x,y,t) = cos(π(1.0-x))
# F(x,y,t) = 0.0

order = 2
method = :cgie

P = VariableCoefficientPDE2D(u₀,kx,ky,order,BoundaryLeft,BoundaryRight,BoundaryUp,BoundaryDown)

println(method)
println("(Δx,Δy)=",Dom.Δx,",",Dom.Δy,"      ","Δt=",Δt,"        ","final time=",t_f)


###
# @benchmark solve($P,$Dom,$Δt,$t_f,$method,source=$F)


soln = solve(P,Dom,Δt,t_f,:cgie,adaptive=false,source=F)
println("plotting")
using Plots
# surface(soln.grid.gridy,soln.grid.gridx,soln.u[2],
#     xlabel="y",ylabel="x",zlabel="Temp",
#     xlims=(0.0,1.0), ylims=(0.0,1.0), zlims=(0.0,1.0))

scatter(1:nx,soln.u[2][:,10],legend=false)

# @time solve(P,Dom,Δt,t_f,:cgie)
# Profile.clear_malloc_data()
# @time solve(P,Dom,Δt,t_f,:cgie)