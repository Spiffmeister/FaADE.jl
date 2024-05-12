using LinearAlgebra
using Revise

using FaADE

# using BenchmarkTools

# using ProfileView
# using Cthulhu
# using Profile


order = 2
K = 1.0

n = 101


Δt = 1.0e-3
# t = 1.0
t = 1.0
# Δt = 1/(n-1)
# nt = round(t/Δt)
# Δt = t/nt


ωx = 9.0
ωt = 1.0
cx = 0.0

K = 1.0

θ = 0.5


# Solution

exact(x,t) = cos(2π*ωt*t) * sin(2π*x*ωx + cx)
u₀(x) = exact(x,0.0)
F(x,t) = -2π*ωt*sin(2π*ωt*t)*sin(2π*x*ωx + cx) + K * 4π^2 * ωx^2 * cos(2π*ωt*t)*sin(2π*x*ωx + cx)
# DIRICHLET
BxL(t) = cos(2π*ωt*t) * sin(cx) #Boundary condition x=0
BxR(t) = cos(2π*ωt*t) * sin(2π*ωx + cx) #Boundary condition x=Lx
# NEUMANN
# BxL(t) = 2π*ωx * K * cos(2π*ωt*t) * cos(cx) #Boundary condition x=0
# BxR(t) = 2π*ωx * K * cos(2π*ωt*t) * cos(2π*ωx + cx) #Boundary condition x=Lx




#====== New solver 1 volume ======#
Dom = Grid1D([0.0,1.0],n)
Dl = FaADE.SATs.SAT_Dirichlet(BxL,Dom.Δx,Left,  order)
Dr = FaADE.SATs.SAT_Dirichlet(BxR,Dom.Δx,Right, order)
BD = FaADE.Inputs.SATBoundaries(Dl,Dr)

P = Problem1D(order,u₀,K,Dom,BD,F,nothing)
println("---Solving 1 volume---")
soln = solve(P,Dom,Δt,t,solver=:theta,θ=θ)
# soln1V = solve(P1V,Dom1V,Δt,t)

# @benchmark solve($P1V,$Dom1V,$Δt,$t)


#====== New solver 2 volume ======#
D1 = Grid1D([0.0,0.5],51)
D2 = Grid1D([0.5,1.0],51)

Dom2V = GridMultiBlock(D1,D2)

Dl = FaADE.SATs.SAT_Dirichlet(BxL,D1.Δx,Left, order)
Dr = FaADE.SATs.SAT_Dirichlet(BxR,D2.Δx,Right,order)
BD = FaADE.Inputs.SATBoundaries(Dl,Dr)
println("---Solving 2 volume---")
P2V = Problem1D(order,u₀,K,Dom2V,BD)

soln2V = solve(P2V,Dom2V,Δt,t,solver=:theta,θ=θ)









# e = [Dom1V.grid[i] for i in eachindex(Dom1V)]
# e = zeros(Dom1V.n)

# println("n=",n," error ",norm(e .- soln.u[2])/norm(e))
# println("n=",n," error ",norm(e .- soln.u[2])*sqrt(Dom.Δx))


#====== New solver 2 volume ======#
D1 = Grid1D([0.0,0.35],35)
D2 = Grid1D([0.35,0.65],31)
D3 = Grid1D([0.65,1.0],35)


Dom3V = GridMultiBlock(D1,D2,D3)

Dl = FaADE.SATs.SAT_Dirichlet(BxL,D1.Δx,Left,    order)
Dr = FaADE.SATs.SAT_Dirichlet(BxR,D3.Δx,Right,   order)
BD = FaADE.Inputs.SATBoundaries(Dl,Dr)

P3V = Problem1D(order,u₀,K,Dom3V,BD)

println("---Solving 2 volume---")
soln3V = solve(P3V,Dom3V,Δt,t)



#= Plotting =#



e = [exact(Dom[i],t) for i in eachindex(Dom)]
u0 = [u₀(Dom[i]) for i in eachindex(Dom)]

using Plots

p1 = plot()
plot!(p1, Dom.grid,e,label="exact")
plot!(p1, Dom.grid,soln.u[2],label="1 vol")
plot!(p1, Dom.grid,u0,label="u₀",linestyle=:dash)
plot!(p1, vcat([Dom.grid[1:51], Dom.grid[51:end]]...),vcat(soln2V.u[2]...),label="new")
plot!(p1, vcat([Dom.grid[1:35], Dom.grid[35:65], Dom.grid[65:end]]...),vcat(soln3V.u[2]...),label="3 vol")

# l = @layout[a; b]
# plot!(p1, Dom1V.grid[prng],[u₀(Dom1V.grid[i]) for i in prng],label="u₀",legend=true,linestyle=:dash)

# p2 = plot()
# plot!(p2,   abs.(soln.u[2][prng] .- e[prng]),label="err")

# plot(p1,p2,layout=l)


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

P1 = Problem1D(order,u₀,K,G,B1)

DBlock = FaADE.solvers.DataMultiBlock(P1,G,0.1,0.0)



Dl = FaADE.SATs.SAT_Dirichlet(x->0.0,sG1.Δx,Left,1,2)
Dr = FaADE.SATs.SAT_Dirichlet(x->0.0,sG2.Δx,Right,1,2)

B = FaADE.SATs.SATBoundaries(Dl,Dr)



P = Problem1D(order,u₀,K,G,B)




s2G1 = Grid2D([0.0,0.5],[0.0,1.0],6,6)
s2G2 = Grid2D([0.5,1.0],[0.0,1.0],11,6)

G2 = FaADE.Helpers.GridMultiBlock([s2G1,s2G2])

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

P3V = Problem1D(order,u₀,K,Dom3V,BD)

println("Solving")
soln3V = solve(P3V,Dom3V,Δt,t)
=#
