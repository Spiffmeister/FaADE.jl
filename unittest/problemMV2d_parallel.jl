
using FaADE
using LinearAlgebra
using BenchmarkTools
# using ProfileView
# using Cthulhu


# Simulation parameters
order = 2
K = 1.0

Δt = 0.01/2/2
t = 1.0
# t = 0.03

# Set initial condition
# u₀(x,y) = x.^2
u₀(x,y) = exp.(-((x-0.5)^2 + (y-0.5)^2) / 0.02)






# Original solver
Dom1V = Grid2D([-0.5,0.5],[-0.5,0.5],41,41)



# Initial condition
u₀(x,y) = 0.0
# Source term
F(x,y,t) = 2π^2*cos(π*x)*cos(π*y)
# Magnetic field
function B(X,x,p,t)
    X[1] = -π*cos(π*x[1])*sin(π*x[2])
    X[2] = π*sin(π*x[1])*cos(π*x[2])
    # X[3] = 0.0
end

Ψ(x,y) = cos(π*x)*cos(π*y)

# Exact solution for k_perp = 1
T(x,y,t) = (1.0 - exp(-2.0*π^2*t) )*Ψ(x,y)





gdata   = construct_grid(B,Dom1V,[-2.0π,2.0π],ymode=:stop)
PData   = ParallelData(gdata,Dom1V,κ=1.0e10)







BoundaryLeft    = Boundary(Dirichlet,(y,t)->0.0,Left,1)
BoundaryRight   = Boundary(Dirichlet,(y,t)->0.0,Right,1)
BoundaryUp      = Boundary(Dirichlet,(y,t)->0.0,Up,2)
BoundaryDown    = Boundary(Dirichlet,(y,t)->0.0,Down,2)
PO1V = VariableCoefficientPDE2D(u₀,(x,y)->1.0,(x,y)->1.0,order,BoundaryLeft,BoundaryRight,BoundaryUp,BoundaryDown)
println("---Solving old---")
solnO1V = solve(PO1V,Dom1V,Δt,t,:cgie,Pgrid=PData,source=F) #-Δt to ensure ends at the same time as new methods

# @benchmark solve($PO1V,$Dom1V,$Δt,100.0-$Δt,:cgie)
# @profview solnO1V = solve(PO1V,Dom1V,Δt,t,:cgie)
# @profview solnO1V = solve(PO1V,Dom1V,Δt,t,:cgie)


# New solver 1 volume
Dl = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,Dom1V.Δx,Left,1,order)
Dr = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,Dom1V.Δx,Right,1,order)
Du = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,Dom1V.Δx,Up,2,order)
Dd = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,Dom1V.Δx,Down,2,order)
BD1V = (Dl,Dr,Du,Dd)

P1V = Problem2D(order,u₀,K,K,Dom1V,BD1V,F,PData)

println("---Solving 1 volume---")
soln1V = solve(P1V,Dom1V,Δt,t,solver=:theta,θ=0.5)

# @benchmark solve($P1V,$Dom1V,$Δt,$t)
# @profview soln1V = solve(P1V,Dom1V,Δt,t)
# @profview soln1V = solve(P1V,Dom1V,Δt,t)




using Plots
surface(Dom1V.gridx,Dom1V.gridy,solnO1V.u[2])
surface(Dom1V.gridx,Dom1V.gridy,soln1V.u[2])




T_exact = zeros(size(Dom1V));
for j = 1:Dom1V.ny
    for i = 1:Dom1V.nx
        T_exact[i,j] = T(Dom1V.gridx[i],Dom1V.gridy[j],solnO1V.t[end])
    end
end
println("Old error ",norm(solnO1V.u[2]-T_exact,Inf))
for j = 1:Dom1V.ny
    for i = 1:Dom1V.nx
        T_exact[i,j] = T(Dom1V.gridx[i],Dom1V.gridy[j],soln1V.t[end])
    end
end
println("CN error ",norm(soln1V.u[2]-T_exact,Inf))


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