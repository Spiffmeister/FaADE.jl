
using FaADE
# using BenchmarkTools
# using ProfileView
# using Cthulhu


order = 2
K = 1.0

Δt = 0.01
t = 10.0

# u₀(x) = x.^2
u₀(x,y) = exp.(-((x-0.5)^2 + (y-0.5)^2) / 0.02)


#=
Dom1V = Grid1D([0.0,1.0],21)

Dl = FaADE.SATs.SAT_Dirichlet(t->0.0,Dom1V.Δx,Left,1,order)
Dr = FaADE.SATs.SAT_Dirichlet(t->1.0,Dom1V.Δx,Right,1,order)
BD1V = FaADE.SATs.SATBoundaries(Dl,Dr)

P1V = newProblem1D(order,u₀,K,Dom1V,BD)

println("Solving")
soln = solve(P1V,Dom1V,Δt,t)
=#

#= =#


D1 = Grid2D([0.0,0.5],[0.0,1.0],11,11)
D2 = Grid2D([0.5,1.0],[0.0,1.0],11,11)

Joints = [[(1,Left),(2,Right),(1,Up),(1,Down)],
            [(1,Left),(2,Right),(2,Up),(2,Down)]]

Dom2V = GridMultiBlock([D1,D2],Joints)

Dl = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,D1.Δx,Left,1,order)
Dr = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,D1.Δx,Right,1,order)
Du = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,D1.Δx,Up,2,order)
Dd = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,D1.Δx,Down,2,order)
# Pu = FaADE.SATs.SAT_Periodic(D1.Δx,2,order,Up)
# Pd = FaADE.SATs.SAT_Periodic(D1.Δx,2,order,Down)
BD = FaADE.SATs.SATBoundaries(Dl,Dr,Du,Dd)

P2V = newProblem2D(order,u₀,K,K,Dom2V,BD)

println("Solving")
@time soln = solve(P2V,Dom2V,Δt,t)


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

#=
println("Solve 2")
BoundaryLeft = Boundary(Dirichlet,t->0.0,Left,1)
BoundaryRight = Boundary(Dirichlet,t->1.0,Right,1)
P = VariableCoefficientPDE1D(u₀,x->1.0,order,BoundaryLeft,BoundaryRight)
soln1b = solve(P,sG1,Δt,t,:cgie)
@benchmark soln1b = solve($P,$sG1,$Δt,$t,:cgie)
=#

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