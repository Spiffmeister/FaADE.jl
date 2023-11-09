
using FaADE
using BenchmarkTools
using ProfileView
using Cthulhu

#=
D1 = Grid1D([0.0,0.5],6)
D2 = Grid1D([0.5,1.0],6)
=#

order = 2
K = 1.0

Δt = 0.01
t = 10.0

sG1 = Grid1D([0.0,1.0],1001)


u₀(x) = x.^2
# u₀(x) = sin.(2π*x*2 .+ 1.0)
ũ(x,t;ωx=1.0,cx=0.0) = cos(2π*t) * sin(2π*x*ωx + cx)
ũ₀(x;ωx=1.0,cx=0.0) = sin(2π*ωx*x + cx)

F(x,t;ωx=1.0,cx=0.0,K=1.0) = 
        -2π*sin(2π*t)*sin(2π*x*ωx + cx) + 
            K * 4π^2 * ωx^2 * cos(2π*t)*sin(2π*x*ωx + cx)



Dl = FaADE.SATs.SAT_Dirichlet(x->0.0,sG1.Δx,Left,1, order)
Dr = FaADE.SATs.SAT_Dirichlet(x->1.0,sG1.Δx,Right,1,order)
# Dl = FaADE.SATs.SAT_Dirichlet(t->sin(1.0),sG1.Δx,Left,1,order)
# Dr = FaADE.SATs.SAT_Dirichlet(t->sin(4π + 1.0),sG1.Δx,Right,1,order)

B1 = FaADE.SATs.SATBoundaries(Dl,Dr)

# PL = FaADE.SATs.SAT_Periodic(sG1.Δx,1,order,Left)
# PR = FaADE.SATs.SAT_Periodic(sG1.Δx,1,order,Right)
# B1 = FaADE.SATs.SATBoundaries(PL,PR)


P1 = newProblem1D(order,u₀,K,sG1,B1)

# DBlock = FaADE.solvers.DataMultiBlock(P1,sG1,0.1,0.0)
# DBlock = FaADE.solvers.MultiDataBlock(P1,sG1)

# CGBlock = FaADE.solvers.ConjGradMultiBlock(sG1,P1.order)

println("Solving")
soln = solve(P1,sG1,Δt,t)
# @profview soln1d_tmpa = solve(P1,sG1,Δt,t)
# @profview soln1d_tmpb = solve(P1,sG1,Δt,t)
# @benchmark solve($P1,$sG1,$Δt,$t)


println("Solve 2")
BoundaryLeft = Boundary(Dirichlet,t->0.0,Left,1)
BoundaryRight = Boundary(Dirichlet,t->1.0,Right,1)
P = VariableCoefficientPDE1D(u₀,x->1.0,order,BoundaryLeft,BoundaryRight)
soln1b = solve(P,sG1,Δt,t,:cgie)
# @benchmark soln1b = solve($P,$sG1,$Δt,$t,:cgie)


#=
order = 2
u₀(x) = x.^3
K = 1.0

Δt = 0.01
t = 10.0

sG1 = Grid1D([0.0,1.0],11)
sG2 = Grid1D([1.0,2.0],6)

G = FaADE.Helpers.GridMultiBlock([sG1,sG2])

Dl = FaADE.SATs.SAT_Dirichlet(x->0.0,sG1.Δx,Left,1,2)
Dr = FaADE.SATs.SAT_Dirichlet(x->1.0,sG1.Δx,Right,1,2)
B1 = FaADE.SATs.SATBoundaries(Dl,Dr)

P1 = newProblem1D(order,u₀,K,G,B1)

DBlock = FaADE.solvers.DataMultiBlock(P1,G,0.1,0.0)



Dl = FaADE.SATs.SAT_Dirichlet(x->0.0,sG1.Δx,Left,1,2)
Dr = FaADE.SATs.SAT_Dirichlet(x->0.0,sG2.Δx,Right,1,2)

B = FaADE.SATs.SATBoundaries(Dl,Dr)



P = newProblem1D(order,u₀,K,G,B)




s2G1 = Grid2D([0.0,0.5],[0.0,1.0],6,6)
s2G2 = Grid2D([0.5,1.0],[0.0,1.0],11,6)

G2 = FaADE.Helpers.GridMultiBlock([s2G1,s2G2])

=#