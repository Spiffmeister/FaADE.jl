using LinearAlgebra
using Printf
# using GLMakie
# pyplot()

using BenchmarkTools
using Profile

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

Δt = 0.1* min(Dom.Δx^2,Dom.Δy^2)
t_f = 1000Δt


#=== MMS ===#

# Initial condition
cx = cy = 0.0

ũ(x,y,t) = cos(2π*t) * sin(2π*x + cx) * sin(2π*y + cy) #Solution

ũ₀(x,y) = sin(2π*x + cx) * sin(2π*y + cy) #Initial condition

Bxũ(y,t) = cos(2π*t) * sin(c_x) * sin(2π*y + c_y) #Boundary condition x
Byũ(x,t) = cos(2π*t) * sin(2π*x + c_x) * sin(c_y) #Boundary condition y

F(x,y,t) = -2π*sin(2π*t)*sin(2π*x+cx)*sin(2π*y+cy) - K * 16π^4 * cos(2π*t)*sin(2π*x+cx)*sin(2π*y+cy) #F = ∂ₜũ - K∇ũ 


BoundaryLeft = BoundaryRight = Boundary(Dirichlet,Bxũ,Left,1)
BoundaryUp = BoundaryDown = Boundary(Dirichlet,Byũ,Up,2)

order = 2
method = :cgie

P = VariableCoefficientPDE2D(u₀,kx,ky,order,BoundaryLeft,BoundaryRight,BoundaryUp,BoundaryDown)

println("(Δx,Δy)=(",Dom.Δx,",",Dom.Δy,")      ","Δt=",Δt,"        ","final time=",t_f,",    solver=",method,".")


###

soln = solve(P,Dom,Δt,t_f,:cgie,adaptive=false,source=F)
println("plotting")
using Plots
surface(soln.grid.gridy,soln.grid.gridx,soln.u[2],
    xlabel="y",ylabel="x",zlabel="Temp",
    xlims=(0.0,1.0), ylims=(0.0,1.0), zlims=(0.0,1.0))

scatter(1:Dom.nx,soln.u[2][:,1:end],legend=false)

# @time solve(P,Dom,Δt,t_f,:cgie)
# Profile.clear_malloc_data()
# @time solve(P,Dom,Δt,t_f,:cgie)