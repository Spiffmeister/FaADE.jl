
using Revise
using LinearAlgebra
using FaADE


plot = false


order = 2
θ = 0.5

Ψ(x,y) = cos(π*x)*cos(π*y)


# Domain
nx = ny = 31
Dom = Grid2D([-0.5,0.5],[-0.5,0.5],nx,ny)


# Homogeneous boundary conditions
BoundaryLeft    = FaADE.SATs.SAT_Dirichlet((y,t) -> cos(0.5π)*cos(π*y)    , Dom.Δx, Left,   order)
BoundaryRight   = FaADE.SATs.SAT_Dirichlet((y,t) -> cos(-0.5π)*cos(π*y)   , Dom.Δx, Right,  order)
BoundaryUp      = FaADE.SATs.SAT_Dirichlet((x,t) -> cos(π*x)*cos(0.5π)    , Dom.Δy, Up,     order)
BoundaryDown    = FaADE.SATs.SAT_Dirichlet((x,t) -> cos(π*x)*cos(-0.5π)   , Dom.Δy, Down,   order)

BC = FaADE.Inputs.SATBoundaries(BoundaryLeft,BoundaryRight,BoundaryUp,BoundaryDown)

# Source term
F(x,t) = 2π^2*cos(π*x[1])*cos(π*x[2])
# Magnetic field
function B(X,x,p,t)
    X[1] = -π*cos(π*x[1])*sin(π*x[2])
    X[2] = π*sin(π*x[1])*cos(π*x[2])
    # X[3] = 0.0
end
# Exact solution
T(x,y,t) = (1.0 - exp(-2.0*π^2*t) )*Ψ(x,y) # k_perp = 1
# Initial condition
u₀(x,y) = 0.0
# u₀(x,y) = T(x,y,0.0)






κ_para = 1.0e10
k = 1.0 #perpendicular diffusion
println("===PARA=",κ_para,"===ORDER=",order,"===")




gdata   = construct_grid(B,Dom,[-2.0π,2.0π],ymode=:stop)

PData   = ParallelData(gdata,Dom,κ=κ_para)





# Build PDE problem
P = Problem2D(order,u₀,k,k,Dom,BC,F,PData)

# Time domain
# Δt = Dom.Δx
Δt = 0.05
t_f = 0.1
nf = round(t_f/Δt)
Δt = t_f/nf


soln = solve(P,Dom,Δt,2.1Δt,solver=:theta,  θ=θ)
soln = solve(P,Dom,Δt,t_f,solver=:theta,    θ=θ)
println("n=",nx,", Δt=",Δt,", θ=",θ)


T_exact = zeros(eltype(Dom),size(Dom));
for I in eachindex(Dom)
    # T_exact[I] = T(Dom[I]...,soln.t[end])
    T_exact[I] = T(Dom[I]...,t_f)
end


# pollution = abs(T(0.0,0.0,soln.t[2]) - soln.u[2][floor(Int,nx/2)+1,floor(Int,ny/2)+1])
pollution = abs(1/soln.u[2][floor(Int,nx/2)+1,floor(Int,ny/2)+1] - 1)
pollution_time = abs(T(0.0,0.0,soln.t[2]) - soln.u[2][floor(Int,nx/2)+1,floor(Int,ny/2)+1])
rel_error = norm(T_exact .- soln.u[2])/norm(T_exact)
println("poll=",pollution," relerr=",rel_error)



if plot
    using GLMakie
    f1 = Figure()
    ax1 = Axis3(f1[1,1])
    surface!(ax1, T_exact)

    f2 = Figure()
    ax2 = Axis3(f2[1,1])
    surface!(ax2, soln.u[2])
end

# f = Figure()
# ax = Axis(f[1,1])
# scatter!(ax,Dom.gridx[:],Dom.gridy[:])
# scatter!(ax,PData.PGrid.Fplane.x[:],PData.PGrid.Fplane.y[:])


# @. u = 1.0/(1.0 - P.κ * ττ * Δt) * 
# ( u - P.κ*Δt*ττ/2.0 * (P.w_f + P.w_b) )

# @benchmark solve($P,$Dom,$Δt,$t_f,solver=:theta,    θ=$θ)


# using ProfileView
# using Cthulhu
# @profview solve(P,Dom,Δt,t_f,solver=:theta,    θ=θ)
# @profview solve(P,Dom,1e-3,t_f,solver=:theta,    θ=θ)
