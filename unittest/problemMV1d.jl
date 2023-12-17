
using FaADE
using LinearAlgebra
using BenchmarkTools
# using ProfileView
# using Cthulhu
# using Profile


order = 2
K = 1.0

n = 321


# Δt = 1.0e-1
t = 1.7
Δt = 1/(n-1)
nt = round(t/Δt)
Δt = t/nt


ωx = 9.0
ωt = 1.0
cx = 1.0

K = 1.0

θ = 0.5


# Solution

exact(x,t) = cos(2π*ωt*t) * sin(2π*x*ωx + cx)
u₀(x) = exact(x,0.0)
F(x,t) = -2π*ωt*sin(2π*ωt*t)*sin(2π*x*ωx + cx) + K * 4π^2 * ωx^2 * cos(2π*ωt*t)*sin(2π*x*ωx + cx)
# DIRICHLET
# BxL(t) = cos(2π*ωt*t) * sin(cx) #Boundary condition x=0
# BxR(t) = cos(2π*ωt*t) * sin(2π*ωx + cx) #Boundary condition x=Lx
# NEUMANN
BxL(t) = 2π*ωx * K * cos(2π*ωt*t) * cos(cx) #Boundary condition x=0
BxR(t) = 2π*ωx * K * cos(2π*ωt*t) * cos(2π*ωx + cx) #Boundary condition x=Lx


# u₀(x) = x^3
# F(x,t) = 6*x
# exact(x,t) = x
# BxL(t) = 0.0
# BxR(t) = 1.0


# u₀(x) = x^2
# F(x,t) = 0.0
# exact(x,t) = x
# BxL(t) = 0.0
# BxR(t) = 1.0


# u₀(x) = exp.(-(x-0.5)^2 / 0.02)
# exact(x,t) = 0.0
# BxL(t) = 0.0
# BxR(t) = 0.0


# Source



#====== Original solver ======#
Dom = Grid1D([0.0,1.0],n)




#====== New solver 1 volume ======#
Dl = FaADE.SATs.SAT_Neumann(BxL,Dom.Δx,Left,1,order)
Dr = FaADE.SATs.SAT_Neumann(BxR,Dom.Δx,Right,1,order)
BD = FaADE.Inputs.SATBoundaries(Dl,Dr)

# Pl = FaADE.SATs.SAT_Periodic(Dom.Δx,1,order,Left)
# Pr = FaADE.SATs.SAT_Periodic(Dom.Δx,1,order,Right)
# BD = FaADE.Inputs.SATBoundaries(Pl,Pr)

P = newProblem1D(order,u₀,K,Dom,BD,F,nothing)
println("---Solving 1 volume---")
soln = solve(P,Dom,Δt,t,solver=:theta,θ=θ)
# soln1V = solve(P1V,Dom1V,Δt,t)

# @benchmark solve($P1V,$Dom1V,$Δt,$t)





e = [exact(Dom[i],t+Δt) for i in eachindex(Dom)]
# e = [Dom1V.grid[i] for i in eachindex(Dom1V)]
# e = zeros(Dom1V.n)
u0 = [u₀(Dom[i]) for i in eachindex(Dom)]

# println("n=",n," error ",norm(e .- soln.u[2])/norm(e))
println("n=",n," error ",norm(e .- soln.u[2])*sqrt(Dom.Δx))



using Plots

prng = 1:n

l = @layout[a; b]
p1 = plot()
plot!(p1, Dom.grid[prng],soln.u[2][prng],label="new")
plot!(p1, Dom.grid[prng],e[prng],label="exact")
# plot!(p1, Dom1V.grid[prng],[u₀(Dom1V.grid[i]) for i in prng],label="u₀",legend=true,linestyle=:dash)

p2 = plot()
plot!(p2,   abs.(soln.u[2][prng] .- e[prng]),label="err")

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


#====== New solver 2 volume ======#
#=
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
