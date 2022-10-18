using LinearAlgebra
using Printf
# using GLMakie
# pyplot()

using BenchmarkTools
using Profile

using Pkg
Pkg.activate(".")
using SBP_operators






#=== MMS ===#
# Setting up the manufactured solution

cx = cy = 0.0
ωy = 1.0
ωx = 1.0
K = 1.0

ũ(x,y,t) = cos(2π*t) * sin(2π*x*ωx + cx) * sin(2π*y*ωy + cy) #Solution

ũ₀(x,y) = sin(2π*x*ωx + cx) * sin(2π*y + cy) #Initial condition

Bxũ(y,t) = 0.0 #cos(2π*t) * sin(cx) * sin(2π*y*ωy + cy) #Boundary condition x
Byũ(x,t) = 0.0 #cos(2π*t) * sin(2π*x*ωx + cx) * sin(cy) #Boundary condition y

F(x,y,t) = -2π*sin(2π*t)*sin(2π*x+cx)*sin(2π*y+cy) + 
    K * 4π^2 * (ωx^2 + ωy^2) * cos(2π*t)*sin(2π*x*ωx+cx)*sin(2π*y*ωy+cy) #F = ∂ₜũ - K∇ũ 


#=== Problem setup ===#
𝒟x = [0.0,1.0]
𝒟y = [0.0,1.0]
nx = 41
ny = 41
Dom = Grid2D(𝒟x,𝒟y,nx,ny)

Δt = 0.1* min(Dom.Δx^2,Dom.Δy^2)
t_f = 1_000Δt

# Diffusion coefficients
kx = ky = zeros(Float64,nx,ny) .+ 1.0;

# Boundary conditions from the MMS
BoundaryLeft = Boundary(Dirichlet,Bxũ,Left,1)
BoundaryRight = Boundary(Dirichlet,Bxũ,Right,1)
BoundaryUp = Boundary(Dirichlet,Byũ,Up,2)
BoundaryDown = Boundary(Dirichlet,Byũ,Down,2)

order = 2
method = :cgie

P = VariableCoefficientPDE2D(ũ₀,kx,ky,order,BoundaryLeft,BoundaryRight,BoundaryUp,BoundaryDown)


###

println("(Δx,Δy)=(",Dom.Δx,",",Dom.Δy,")      ","Δt=",Δt,"        ","final time=",t_f,",    solver=",method,".")
soln = solve(P,Dom,Δt,t_f,:cgie,adaptive=false,source=F)
println("plotting")
using Plots

# surface(soln.grid.gridy,soln.grid.gridx,soln.u[2],
#     xlabel="y",ylabel="x",zlabel="Temp",
#     xlims=(0.0,2π), ylims=(0.0,2π))


u_MMS = zeros(Dom.nx,Dom.ny);
for j = 1:Dom.ny
    for i = 1:Dom.nx
        u_MMS[i,j] = ũ(soln.grid.gridx[i],soln.grid.gridy[j],t_f)
    end
end
surface(soln.grid.gridy,soln.grid.gridx,u_MMS.-soln.u[2])

# surface(soln.grid.gridy,soln.grid.gridx,u_MMS)
# surface(soln.grid.gridy,soln.grid.gridx,soln.u[2])

# scatter(1:Dom.nx,soln.u[2][:,1:end],legend=false)

# @time solve(P,Dom,Δt,t_f,:cgie)
# Profile.clear_malloc_data()
# @time solve(P,Dom,Δt,t_f,:cgie)