
using FaADE
using LinearAlgebra
using Plots
# using BenchmarkTools
# using ProfileView
# using Cthulhu


# Simulation parameters
order = 2
K = 1.0


nx = ny = 21

# Δt = 1e-3
t = 0.5
Δt = 1/(nx-1) * 1/(ny-1) * 0.1
nt = round(t/Δt)
Δt = t/nt

ωt = 1.0
ωx = 1.0
ωy = 1.0
cx = 0.0
cy = 0.0

Kx = 1.0
Ky = 1.0

θ = 1.0


# Solution
exact(t,x,y) = cos(2π*ωt*t) * sin(2π*x*ωx + cx) * sin(2π*y*ωy + cy)
u₀(x,y) = exact(0.0,x,y)
F(t,x,y) = -2π*ωt*sin(2π*ωt*t)*sin(2π*x*ωx + cx)*sin(2π*y*ωy + cy) + 
    K * 4π^2 * ωx^2 * cos(2π*ωt*t)*sin(2π*x*ωx + cx)*sin(2π*y*ωy + cy) + 
    K * 4π^2 * ωy^2 * cos(2π*ωt*t)*sin(2π*x*ωx + cx)*sin(2π*y*ωy + cy)
# DIRICHLET
BxL(t,y) = cos(2π*ωt*t) * sin(cx) * sin(2π*ωy*y + cy)           #Boundary condition x=0
BxR(t,y) = cos(2π*ωt*t) * sin(2π*ωx + cx) * sin(2π*ωy*y + cy)   #Boundary condition x=Lx
ByL(t,x) = cos(2π*ωt*t) * sin(2π*ωx*x + cx) * sin(cy)           #Boundary condition y=0
ByR(t,x) = cos(2π*ωt*t) * sin(2π*ωx*x + cx) * sin(2π*ωy + cy)   #Boundary condition y=Ly



# exact(t,x,y) = cos(2π*ωt*t) * sin(2π*x*ωx + cx)
# u₀(x,y) = exact(0.0,x,y)
# F(t,x,y) = -2π*ωt*sin(2π*ωt*t)*sin(2π*x*ωx + cx) + 
#     K * 4π^2 * ωx^2 * cos(2π*ωt*t)*sin(2π*x*ωx + cx)

#=
exact(t,x,y) = exp(-((x-0.5)^2 + (y-0.5)^2) / 0.02)
u₀(x,y) = exact(0.0,x,y)
BxL(t,y) = 0.0
BxR(t,y) = 0.0
ByL(t,x) = 0.0
ByR(t,x) = 0.0
=#



# Original solver
Dom = Grid2D([0.0,1.0],[0.0,1.0],nx,ny)

# New solver 1 volume
# Dl = FaADE.SATs.SAT_Dirichlet(BxL,Dom.Δx,Left,1,order)
# Dr = FaADE.SATs.SAT_Dirichlet(BxR,Dom.Δx,Right,1,order)
# Du = FaADE.SATs.SAT_Dirichlet(ByL,Dom.Δx,Up,2,order)
# Dd = FaADE.SATs.SAT_Dirichlet(ByR,Dom.Δx,Down,2,order)
# BD = FaADE.Inputs.SATBoundaries(Dl,Dr,Du,Dd)

Pl = FaADE.SATs.SAT_Periodic(Dom.Δx,1,order,Left)
Pr = FaADE.SATs.SAT_Periodic(Dom.Δx,1,order,Right)
Pu = FaADE.SATs.SAT_Periodic(Dom.Δy,2,order,Up)
Pd = FaADE.SATs.SAT_Periodic(Dom.Δy,2,order,Down)
BD = FaADE.Inputs.SATBoundaries(Pl,Pr,Pu,Pd)

P = newProblem2D(order,u₀,K,K,Dom,BD,F,nothing)

println("---Solving 1 volume---")
soln = solve(P,Dom,Δt,Δt)
soln = solve(P,Dom,Δt,t)

# @benchmark solve($P1V,$Dom1V,$Δt,$t)

# @profview soln1V = solve(P1V,Dom1V,Δt,t)
# @profview soln1V = solve(P1V,Dom1V,Δt,t)



#=
# New solover 2 volume
D1 = Grid2D([0.0,0.5],[0.0,1.0],21,41)
D2 = Grid2D([0.5,1.0],[0.0,1.0],21,41)

joints = (Joint(2,Up),Joint(1,Down))

Dom2V = GridMultiBlock((D1,D2),joints)

Dl = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,D1.Δx,Left,1,order)
Dr = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,D1.Δx,Right,1,order)
Du = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,D1.Δx,Up,2,order)
Dd = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,D1.Δx,Down,2,order)
BD = FaADE.Inputs.SATBoundaries(Dl,Dr,Du,Dd)

# BCs = [(1,Left,Dl),(1,Up,Du),(1,Down,Dd),(2,Right,Dr),(2,Up,Du),(2,Down,Dd)]

P2V = newProblem2D(order,u₀,K,K,Dom2V,BD)

println("---Solving 2 volume---")
soln2V = solve(P2V,Dom2V,Δt,t)
@benchmark solve($P2V,$Dom2V,$Δt,$t)
=#






e = [exact(soln.t[2],Dom.gridx[i,j],Dom.gridy[i,j]) for i in 1:Dom.nx, j in 1:Dom.ny]

E = zeros(size(Dom))
for I in eachindex(Dom)
    E[I] = exact(soln.t[2]+Δt,Dom[I]...)
end
u0 = [u₀(Dom.gridx[i,j],Dom.gridy[i,j]) for i in 1:Dom.nx, j in 1:Dom.ny]



println("n=",nx," error ",norm(e .- soln.u[2])*sqrt(Dom.Δx*Dom.Δy))
surface(Dom.gridx,Dom.gridy,soln.u[2])


l = @layout[a b]

p1 = surface(Dom.gridx,Dom.gridy,soln.u[2],colorbar=false)
p2 = surface(Dom.gridx,Dom.gridy,e,colorbar=false)
p12 = plot(p1,p2)


# p3 = surface(Dom.gridx,Dom.gridy,soln.u[2] .- e)


#=
function χ_h!(χ,x::Array{Float64},t)
    χ[2] = x[1] #p_1            qdot        θ
    χ[1] = 0.0  #q_1        pdot        ψ
end

dH(X,x,p,t) = χ_h!(X,x,t)
PGrid = FaADE.construct_grid(dH,Dom,[-2π,2π])
Pfn = FaADE.generate_parallel_penalty(PGrid,Dom,2)

P2VP = newProblem2D(order,u₀,K,K,Dom2V,BD,Pfn)
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

P3V = newProblem1D(order,u₀,K,Dom3V,BD)


println("Solving")
@time soln = solve(P3V,Dom3V,Δt,t)
=#



# DBlock = FaADE.solvers.DataMultiBlock(P1,sG1,0.1,0.0)
# DBlock = FaADE.solvers.MultiDataBlock(P1,sG1)

# CGBlock = FaADE.solvers.ConjGradMultiBlock(sG1,P1.order)

# @profview soln1d = solve(P1,sG1,Δt,t)
# @profview soln1d = solve(P1,sG1,Δt,t)
# @benchmark solve($P1,$sG1,$Δt,$t)
