
using FaADE
using BenchmarkTools
using ProfileView
using Cthulhu
using Profile


order = 2
K = 1.0

Δt = 0.01
t = 0.02
# t = 100.0

# u₀(x) = x.^2
u₀(x) = exp.(-(x-0.5)^2 / 0.02)






Dom1V = Grid1D([0.0,1.0],1001)
Dl = FaADE.SATs.SAT_Dirichlet(t->0.0,Dom1V.Δx,Left,1,order)
Dr = FaADE.SATs.SAT_Dirichlet(t->0.0,Dom1V.Δx,Right,1,order)
BD1V = FaADE.SATs.SATBoundaries(Dl,Dr)
P1V = newProblem1D(order,u₀,K,Dom1V,BD1V)
println("Solving")
solnP1V = solve(P1V,Dom1V,Δt,t)
# @benchmark solve($P1V,$Dom1V,$Δt,$t)

# DBlock = FaADE.solvers.DataMultiBlock(P1V,Dom1V,0.0,0.0)
# @code_warntype FaADE.solvers.fillBuffer(:u,DBlock,1,Left)
# @code_warntype DBlock[1].boundary[Left]



u₀(x) = exp.(-(x.-0.5).^2 ./ 0.02)
BoundaryLeft = Boundary(Dirichlet,t->0.0,Left,1)
BoundaryRight = Boundary(Dirichlet,t->0.0,Right,1)
P = VariableCoefficientPDE1D(u₀,t->K,order,BoundaryLeft,BoundaryRight)
solnO1V = solve(P,Dom1V,Δt,t,:cgie)




#= =#

#=
D1 = Grid1D([0.0,0.5],501)
D2 = Grid1D([0.5,1.0],501)

Joints = [[(2,Right)],
            [(1,Left)]]

Dom2V = GridMultiBlock([D1,D2],Joints)

Dl = FaADE.SATs.SAT_Dirichlet(t->0.0,D1.Δx,Left,1,order)
Dr = FaADE.SATs.SAT_Dirichlet(t->1.0,D2.Δx,Right,1,order)
BD = FaADE.SATs.SATBoundaries(Dl,Dr)

P2V = newProblem1D(order,u₀,K,Dom2V,BD)


println("Solving")
# Profile.clear_malloc_data()
# soln = solve(P2V,Dom2V,Δt,t)
# Profile.clear_malloc_data()
soln = solve(P2V,Dom2V,Δt,t)
# @profview soln_tmpa = solve(P2V,Dom2V,Δt,t)
# @profview soln_tmpb = solve(P2V,Dom2V,Δt,t)
@benchmark solve($P2V,$Dom2V,$Δt,$t)
=#




# DBlock = FaADE.solvers.DataMultiBlock(P2V,Dom2V,0.0,0.0)
# @code_warntype FaADE.solvers.fillBuffer(:u,DBlock,1,Left)
# @code_warntype DBlock[1].boundary[Left]

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