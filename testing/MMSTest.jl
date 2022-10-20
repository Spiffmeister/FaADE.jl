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



#=== Define a function to generate the MMS solution ===#
function generate_MMS(MMS::Function,grid::SBP_operators.Helpers.Grid2D,t::Float64)
    u_MMS = zeros(grid.nx,grid.ny)
    for j = 1:grid.ny
        for i = 1:grid.nx
            u_MMS[i,j] = MMS(grid.gridx[i],grid.gridy[j],t)
        end
    end
    return u_MMS
end


#=== Problem setup ===#
𝒟x = [0.0,1.0]
𝒟y = [0.0,1.0]
# Boundary conditions from the MMS
BoundaryLeft = Boundary(Dirichlet,Bxũ,Left,1)
BoundaryRight = Boundary(Dirichlet,Bxũ,Right,1)
BoundaryUp = Boundary(Dirichlet,Byũ,Up,2)
BoundaryDown = Boundary(Dirichlet,Byũ,Down,2)

order = 2
method = :cgie

npts = [11,21,31,41,51,61,71,81]
comp_soln = []
MMS_soln = []
grids = []
relerr = []
for n in npts
    Dom = Grid2D(𝒟x,𝒟y,n,n)
    
    Δt = 0.01*Dom.Δx^2
    t_f = 0.1

    # Diffusion coefficients
    kx = ky = zeros(Float64,n,n) .+ 1.0;

    P = VariableCoefficientPDE2D(ũ₀,kx,ky,order,BoundaryLeft,BoundaryRight,BoundaryUp,BoundaryDown)

    println("Solving n=",Dom.nx," case with Δt=",Δt)
    soln = solve(P,Dom,Δt,t_f,:cgie,source=F)
    
    u_MMS = generate_MMS(ũ,Dom,t_f)

    push!(comp_soln,soln)
    push!(grids,Dom)
    push!(MMS_soln,u_MMS)
    push!(relerr, norm(u_MMS .- soln.u[2])/norm(MMS_soln))
end


log.(relerr[1:end-1]./relerr[2:end]) ./ log.( (1 ./ (npts[1:6-1].-1))./(1 ./ (npts[2:6].-1) ))


# println("plotting")
# using Plots

# surface(soln.grid.gridy,soln.grid.gridx,soln.u[2],
#     xlabel="y",ylabel="x",zlabel="Temp",
#     xlims=(0.0,2π), ylims=(0.0,2π))





# surface(soln.grid.gridy,soln.grid.gridx,u_MMS)
# surface(soln.grid.gridy,soln.grid.gridx,soln.u[2])

# scatter(1:Dom.nx,soln.u[2][:,1:end],legend=false)

# @time solve(P,Dom,Δt,t_f,:cgie)
# Profile.clear_malloc_data()
# @time solve(P,Dom,Δt,t_f,:cgie)