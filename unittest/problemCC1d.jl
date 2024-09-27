
using FaADE
using BenchmarkTools
# using ProfileView
# using Cthulhu
# using Profile


order = 2
K = 1.0

Δt = 0.01
t = 1.0
# t = 10.0

u₀(x) = x.^2
# u₀(x) = exp.(-(x-0.5)^2 / 0.02)

x = collect(LinRange(0.0,1.0,11))
x = x.^2

Grid1D()


#=
#====== Original solver ======#
Dom1V = Grid1D([0.0,1.0],101)

# u₀(x) = exp.(-(x.-0.5).^2 ./ 0.02)
BoundaryLeft = Boundary(Dirichlet,t->0.0,Left,1)
BoundaryRight = Boundary(Dirichlet,t->1.0,Right,1)
P = VariableCoefficientPDE1D(u₀,t->K,order,BoundaryLeft,BoundaryRight)
println("---Solving old---")
solnO1V = solve(P,Dom1V,Δt,t,:cgie)

# @benchmark solve($P,$Dom1V,$Δt,$t,:cgie)

# @profview solnO_tmpa = solve(P2V,Dom2V,Δt,t)
# @profview solnO_tmpb = solve(P2V,Dom2V,Δt,t)
=#




#====== New solver 1 volume ======#
Dl = FaADE.SATs.SAT_Dirichlet(t->0.0,Dom1V.Δx,Left,1,order)
Dr = FaADE.SATs.SAT_Dirichlet(t->1.0,Dom1V.Δx,Right,1,order)
BD1V = (Dl,Dr)
P1V = newProblem1D(order,u₀,K,Dom1V,BD1V)
println("---Solving 1 volume---")
soln1V = solve(P1V,Dom1V,Δt,t,solver=:theta,θ=3/4)
# soln1V = solve(P1V,Dom1V,Δt,t,solver=:cn)

# @benchmark solve($P1V,$Dom1V,$Δt,$t)



using Plots
plot(Dom1V.grid,solnO1V.u[2])

plot!(Dom1V.grid,soln1V.u[2])


