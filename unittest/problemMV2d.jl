using LinearAlgebra
using Revise
using FaADE
# using BenchmarkTools
# using ProfileView
# using Cthulhu


# Simulation parameters
order = 2
K = 1.0

nx = ny = 41

Δt = 1e-3
t = Δt

ωt = 1.0
ωx = 1.0
ωy = 1.0
cx = 0.0
cy = 0.0

Kx = 1.0
Ky = 1.0

θ = 0.5

#=
# 2D sine wave solution
# Solution
exact(x,y,t) = cos(2π*ωt*t) * sin(2π*x*ωx + cx) * sin(2π*y*ωy + cy)
u₀(x,y) = exact(x,y,0.0)
F(X,t) = begin
    x,y = X
    -2π*ωt*sin(2π*ωt*t)*sin(2π*x*ωx + cx)*sin(2π*y*ωy + cy) + 
    K * 4π^2 * ωx^2 * cos(2π*ωt*t)*sin(2π*x*ωx + cx)*sin(2π*y*ωy + cy) + 
    K * 4π^2 * ωy^2 * cos(2π*ωt*t)*sin(2π*x*ωx + cx)*sin(2π*y*ωy + cy)
end
# DIRICHLET
BxL(y,t) = cos(2π*ωt*t) * sin(cx) * sin(2π*ωy*y + cy)           #Boundary condition x=0
BxR(y,t) = cos(2π*ωt*t) * sin(2π*ωx + cx) * sin(2π*ωy*y + cy)   #Boundary condition x=Lx
ByL(x,t) = cos(2π*ωt*t) * sin(2π*ωx*x + cx) * sin(cy)           #Boundary condition y=0
ByR(x,t) = cos(2π*ωt*t) * sin(2π*ωx*x + cx) * sin(2π*ωy + cy)   #Boundary condition y=Ly
# NEUMANN
# BxLũ(y,t) = 2π*ωx * K * cos(2π*t) * cos(cx)             * sin(2π*y*ωy + cy) #Boundary condition x=0
# BxRũ(y,t) = 2π*ωx * K * cos(2π*t) * cos(2π*ωx + cx)  * sin(2π*y*ωy + cy) #Boundary condition x=Lx
# ByLũ(x,t) = 2π*ωy * K * cos(2π*t) * sin(2π*x*ωx + cx)   * cos(cy) #Boundary condition y=0
# ByRũ(x,t) = 2π*ωy * K * cos(2π*t) * sin(2π*x*ωx + cx)   * cos(2π*ωy + cy) #Boundary condition y=Ly
=#


# Plane wave solution
exact(x,y,t) = cos(2π*ωt*t) * sin(2π*x*ωx + cx)
u₀(x,y) = exact(x,y,0.0)
F(X,t) = begin
    x,y = X
    -2π*ωt*sin(2π*ωt*t)*sin(2π*x*ωx + cx) + 
    K * 4π^2 * ωx^2 * cos(2π*ωt*t)*sin(2π*x*ωx + cx)
end

# DIRICHLET
BxL(y,t) = cos(2π*ωt*t) * sin(cx)           #Boundary condition x=0
BxR(y,t) = cos(2π*ωt*t) * sin(2π*ωx + cx)   #Boundary condition x=Lx
ByL(x,t) = cos(2π*ωt*t) * sin(2π*ωx*x + cx) #Boundary condition y=0
ByR(x,t) = cos(2π*ωt*t) * sin(2π*ωx*x + cx) #Boundary condition y=Ly

#=
# Linear solution
exact(x,y,t) = x
u₀(x,y) = exact(x,y,0.0)
F = nothing

BxL(y,t) = 0.0
BxR(y,t) = 1.0
ByL(x,t) = x
ByR(x,t) = x
=#


#====== New solver 1 volume ======#
Dom = Grid2D([0.0,1.0],[0.0,1.0],nx,ny)

# New solver 1 volume
Dl = FaADE.SATs.SAT_Dirichlet(BxL,Dom.Δx,Left,order)
Dr = FaADE.SATs.SAT_Dirichlet(BxR,Dom.Δx,Right,order)
Du = FaADE.SATs.SAT_Dirichlet(ByL,Dom.Δy,Up,order)
Dd = FaADE.SATs.SAT_Dirichlet(ByR,Dom.Δy,Down,order)
BD = FaADE.Inputs.SATBoundaries(Dl,Dr,Du,Dd)

# Pl = FaADE.SATs.SAT_Periodic(Dom.Δx,1,order,Left)
# Pr = FaADE.SATs.SAT_Periodic(Dom.Δx,1,order,Right)
# Pu = FaADE.SATs.SAT_Periodic(Dom.Δy,2,order,Up)
# Pd = FaADE.SATs.SAT_Periodic(Dom.Δy,2,order,Down)
# BD = FaADE.Inputs.SATBoundaries(Pl,Pr,Pu,Pd)

P = Problem2D(order,u₀,K,K,Dom,BD,F,nothing)
println("---Solving 1 volume---")
soln = solve(P,Dom,Δt,t)


#====== New solver 1 volume ======#
println("2 volume")

# New solover 2 volume
D1 = Grid2D([0.0,0.5],[0.0,1.0],21,41)
D2 = Grid2D([0.5,1.0],[0.0,1.0],21,41)

joints = ((Joint(2,Right),),
            (Joint(1,Left),))

Dom2V = GridMultiBlock((D1,D2),joints)

Dl = FaADE.SATs.SAT_Dirichlet(BxL,D1.Δx,Left,    order)
Dr = FaADE.SATs.SAT_Dirichlet(BxR,D2.Δx,Right,   order)

Du1 = FaADE.SATs.SAT_Dirichlet(ByL,D1.Δy,Up,      order)
Dd1 = FaADE.SATs.SAT_Dirichlet(ByR,D1.Δy,Down,    order)
Du2 = FaADE.SATs.SAT_Dirichlet(ByL,D2.Δy,Up,      order)
Dd2 = FaADE.SATs.SAT_Dirichlet(ByR,D2.Δy,Down,    order)

BD2 = Dict(1 => (Dl,Du1,Dd1), 2 => (Dr,Du2,Dd2))

# Pl = FaADE.SATs.SAT_Periodic(Dom.Δx,1,order,Left)
# Pr = FaADE.SATs.SAT_Periodic(Dom.Δx,1,order,Right)
# BD2 = Dict(1 => (Pl,Du,Dd), 2 => (Pr,Du,Dd))


P2V = Problem2D(order,u₀,K,K,Dom2V,BD2,F,nothing)

println("---Solving 2 volume---")
soln2V = solve(P2V,Dom2V,Δt,t)



using GLMakie

e = zeros(size(Dom))
for I in eachindex(Dom)
    e[I] = exact(Dom[I]...,t)
end
u0 = [u₀(Dom[i]...) for i in eachindex(Dom)]


f = Figure()
ax1 = Axis3(f[1,1])
# surface!(ax1,Dom.gridx,Dom.gridy,e,colorbar=false)
surface!(ax1,Dom.gridx,Dom.gridy,soln.u[2],colorbar=false)
# surface!(ax1,Dom.gridx,Dom.gridy,soln.u[2].-e,colorbar=false)

ax2 = Axis3(f[1,2])
surface!(ax2,Dom2V.Grids[1].gridx,Dom2V.Grids[1].gridy,soln2V.u[2][1],colorbar=false)
surface!(ax2,Dom2V.Grids[2].gridx,Dom2V.Grids[2].gridy,soln2V.u[2][2],colorbar=false)

f



# p3 = surface(Dom.gridx,Dom.gridy,soln.u[2] .- e)


#=
function χ_h!(χ,x::Array{Float64},t)
    χ[2] = x[1] #p_1            qdot        θ
    χ[1] = 0.0  #q_1        pdot        ψ
end

dH(X,x,p,t) = χ_h!(X,x,t)
PGrid = FaADE.construct_grid(dH,Dom,[-2π,2π])
Pfn = FaADE.generate_parallel_penalty(PGrid,Dom,2)

P2VP = Problem2D(order,u₀,K,K,Dom2V,BD,Pfn)
soln = solve(P2VP,Dom2V,Δt,t)
@benchmark solve($P2VP,$Dom2V,$Δt,$t)
=#



#=
using GLMakie
surface(D1.gridx,D1.gridy,soln.u[2][1])
surface!(D2.gridx,D2.gridy,soln.u[2][2])
=#

#=
D1 = Grid1D([0.0,0.35],8)
D2 = Grid1D([0.35,0.65],7)
D3 = Grid1D([0.65,1.0],8)

Joints = [[(2,Right)],
            [(1,Left),(3,Right)],
            [(2,Left)]]

Dom3V = GridMultiBlock([D1,D2,D3],Joints)

Dl = FaADE.SATs.SAT_Dirichlet(t->0.0,D1.Δx,Left,1,order)
Dr = FaADE.SATs.SAT_Dirichlet(t->1.0,D3.Δx,Right,1,order)
BD = FaADE.SATs.SATBoundaries(Dl,Dr)

P3V = Problem1D(order,u₀,K,Dom3V,BD)


println("Solving")
@time soln = solve(P3V,Dom3V,Δt,t)
=#



# DBlock = FaADE.solvers.DataMultiBlock(P1,sG1,0.1,0.0)
# DBlock = FaADE.solvers.MultiDataBlock(P1,sG1)

# CGBlock = FaADE.solvers.ConjGradMultiBlock(sG1,P1.order)

# @profview soln1d = solve(P1,sG1,Δt,t)
# @profview soln1d = solve(P1,sG1,Δt,t)
# @benchmark solve($P1,$sG1,$Δt,$t)
