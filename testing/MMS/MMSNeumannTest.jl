using LinearAlgebra
using Printf
# using GLMakie
# pyplot()

using JLD2

using Pkg
Pkg.activate(".")
using SBP_operators






#=== MMS ===#
# Setting up the manufactured solution

cx = 1.0
cy = 0.0
ωx = 2.5
ωy = 2.0
K = 1.0


ũ(x,y,t) = cos(2π*t) * sin(2π*x*ωx + cx) * sin(2π*y*ωy + cy) #Solution

ũ₀(x,y) = sin(2π*ωx*x + cx) * sin(2π*ωy*y + cy) #Initial condition

BxLũ(y,t) =         2π*ωx * cos(2π*t) * cos(cx)             * sin(2π*y*ωy + cy) #Boundary condition x=0
BxRũ(y,t;Lx=1.0) =  2π*ωx * cos(2π*t) * cos(2π*Lx*ωx + cx)  * sin(2π*y*ωy + cy) #Boundary condition x=Lx
ByLũ(x,t) =         2π*ωy * cos(2π*t) * sin(2π*x*ωx + cx)   * cos(cy) #Boundary condition y=0
ByRũ(x,t;Ly=1.0) =  2π*ωy * cos(2π*t) * sin(2π*x*ωx + cx)   * cos(2π*Ly*ωy + cy) #Boundary condition y=Ly



F(x,y,t) = -2π*sin(2π*t)*sin(2π*x*ωx + cx)*sin(2π*y*ωy + cy) + 
    K * 4π^2 * (ωx^2 + ωy^2) * cos(2π*t)*sin(2π*x*ωx + cx)*sin(2π*y*ωy + cy) #F = ∂ₜũ - K∇ũ 




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
BoundaryLeft = Boundary(Neumann,BxLũ,Left,1)
BoundaryRight = Boundary(Neumann,(y,t) -> BxRũ(y,t,Lx=𝒟x[2]),Right,1)
BoundaryUp = Boundary(Neumann,ByLũ,Up,2)
BoundaryDown = Boundary(Neumann,(x,t) -> ByRũ(x,t,Ly=𝒟y[2]),Down,2)

order = 4
method = :cgie

npts = [21,31,41,51]
# npts = [41]
comp_soln = []
MMS_soln = []
grids = []
relerr = []
for n in npts
    Dom = Grid2D(𝒟x,𝒟y,n,n)
    
    Δt = 0.01*Dom.Δx^2
    t_f = 0.01
    # t_f = 200Δt

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

conv_rate = log.(relerr[1:end-1]./relerr[2:end]) ./ log.( (1 ./ (npts[1:end-1].-1))./(1 ./ (npts[2:end].-1) ))

println("The convergence rate of this MMS setup is: ",conv_rate," for order ",order," SBP operators.")


# println("plotting")
using Plots

# # l = @layout [a b c]
# p = surface(grids[end].gridy,grids[end].gridx,comp_soln[end].u[2],
#     #layout=l,
#     reuse=false,
#     xlabel="y",ylabel="x",zlabel="Solution",
#     xlims=(grids[end].gridx[1],grids[end].gridx[end]), ylims=(grids[end].gridy[1],grids[end].gridy[end]))

# surface(#p[2],
#     grids[end].gridy,grids[end].gridx,MMS_soln[end],
#     reuse=false,
#     xlabel="y",ylabel="x",zlabel="MMS Solution",
#     xlims=(grids[end].gridx[1],grids[end].gridx[end]), ylims=(grids[end].gridy[1],grids[end].gridy[end]))

surface(#p[3],
    grids[end].gridy,grids[end].gridx,(comp_soln[end].u[2].-MMS_soln[end]),
    xlabel="y",ylabel="x",zlabel="Absolute error",
    xlims=(grids[end].gridx[1],grids[end].gridx[end]), ylims=(grids[end].gridy[1],grids[end].gridy[end]))


# surface(grids[end].gridy,grids[end].gridx,MMS_soln[end])
# surface(grids[end].gridy,grids[end].gridx,comp_soln[end].u[2])
# surface(grids[1].gridy,grids[1].gridx,comp_soln[1].u[2] .- MMS_soln[end])

# scatter(1:Dom.nx,soln.u[2][:,1:end],legend=false)

# @time solve(P,Dom,Δt,t_f,:cgie)
# Profile.clear_malloc_data()
# @time solve(P,Dom,Δt,t_f,:cgie)
