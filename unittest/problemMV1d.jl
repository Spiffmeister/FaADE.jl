
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

# u₀(x) = x.^2
# u₀(x) = exp.(-(x-0.5)^2 / 0.02)

ωx = 15.5
cx = 0.0

K = 1.0


# Solution
exact(x,t) = cos(2π*t) * sin(2π*x*ωx + cx)
# Initial condition
u₀(x) = sin(2π*ωx*x + cx)
# Source
F(x,t) = -2π*sin(2π*t)*sin(2π*x*ωx + cx) + K * 4π^2 * ωx^2 * cos(2π*t)*sin(2π*x*ωx + cx)
            
BxL(t) = cos(2π*t) * sin(cx) #Boundary condition x=0
BxR(t) = cos(2π*t) * sin(2π*ωx + cx) #Boundary condition x=Lx

#====== Original solver ======#
Dom1V = Grid1D([0.0,1.0],501)

# u₀(x) = exp.(-(x.-0.5).^2 ./ 0.02)
BoundaryLeft = Boundary(Dirichlet,BxL,Left,1)
BoundaryRight = Boundary(Dirichlet,BxR,Right,1)
P = VariableCoefficientPDE1D(u₀,t->K,order,BoundaryLeft,BoundaryRight)
println("---Solving old---")
solnO1V = solve(P,Dom1V,Δt,t,source=F,:cgie)

# @benchmark solve($P,$Dom1V,$Δt,$t,:cgie)

# @profview solnO_tmpa = solve(P2V,Dom2V,Δt,t)
# @profview solnO_tmpb = solve(P2V,Dom2V,Δt,t)





#====== New solver 1 volume ======#
Dl = FaADE.SATs.SAT_Dirichlet(BxL,Dom1V.Δx,Left,1,order)
Dr = FaADE.SATs.SAT_Dirichlet(BxR,Dom1V.Δx,Right,1,order)
BD1V = FaADE.Inputs.SATBoundaries(Dl,Dr)
P1V = newProblem1D(order,u₀,K,Dom1V,BD1V,F,nothing)
println("---Solving 1 volume---")
soln1V = solve(P1V,Dom1V,Δt,t,solver=:theta,θ=0.5)
# soln1V = solve(P1V,Dom1V,Δt,t)

# @benchmark solve($P1V,$Dom1V,$Δt,$t)


#=
#====== New solver 2 volume ======#
D1 = Grid1D([0.0,0.5],501)
D2 = Grid1D([0.5,1.0],501)

Dom2V = GridMultiBlock(D1,D2)

Dl = FaADE.SATs.SAT_Dirichlet(t->0.0,D1.Δx,Left,1,order)
Dr = FaADE.SATs.SAT_Dirichlet(t->1.0,D2.Δx,Right,1,order)
BD = FaADE.Inputs.SATBoundaries(Dl,Dr)
println("---Solving 2 volume---")
P2V = newProblem1D(order,u₀,K,Dom2V,BD)

soln2V = solve(P2V,Dom2V,Δt,t)
=#



# @benchmark solve($P2V,$Dom2V,$Δt,$t)

# solve(P2V,Dom2V,Δt,t)
# Profile.clear_malloc_data()
# solve(P2V,Dom2V,Δt,t)

# @profview soln_tmpa = solve(P2V,Dom2V,Δt,t)
# @profview soln_tmpb = solve(P2V,Dom2V,Δt,t)






#=
using LinearAlgebra
norm(solnP1V.u[2])
norm(vcat(soln.u[2][1],soln.u[2][2][2:end]))
=#



#=
D1 = Grid1D([0.0,0.35],351)
D2 = Grid1D([0.35,0.65],301)
D3 = Grid1D([0.65,1.0],351)


Dom3V = GridMultiBlock((D1,D2,D3))

Dl = FaADE.SATs.SAT_Dirichlet(t->0.0,D1.Δx,Left,1,order)
Dr = FaADE.SATs.SAT_Dirichlet(t->1.0,D3.Δx,Right,1,order)
BD = FaADE.SATs.SATBoundaries(Dl,Dr)

P3V = newProblem1D(order,u₀,K,Dom3V,BD)

println("Solving")
soln3V = solve(P3V,Dom3V,Δt,t)
=#

e = [exact(Dom1V.grid[i],solnO1V.t[2]) for i in eachindex(Dom1V)]

using Plots


l = @layout[a; b]
p1 = plot(Dom1V.grid,solnO1V.u[2],label="old")
plot!(p1, Dom1V.grid,soln1V.u[2],label="new")
plot!(p1, Dom1V.grid,e,label="exact")


p2 = plot(solnO1V.u[2] .- e,label="old err")
plot!(p2,soln1V.u[2] .- e,label="new err")
plot!(p2,soln1V.u[2] .- solnO1V.u[2],linestyle=:dash,label="new old err")

plot(p1,p2,layout=l)

# plot!(Dom2V.Grids[1].grid,soln2V.u[2][1])
# plot!(Dom2V.Grids[2].grid,soln2V.u[2][2])
#=
# plot!(Dom3V.Grids[1].grid,soln3V.u[2][1])
# plot!(Dom3V.Grids[2].grid,soln3V.u[2][2])
# plot!(Dom3V.Grids[3].grid,soln3V.u[2][3])
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