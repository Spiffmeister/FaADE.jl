using LinearAlgebra

using BenchmarkTools
# using ProfileView
# using Cthulhu

using FaADE




###
𝒟 = [0.0,1.0]
n = 101
Dom = Grid1D(𝒟,n)

K(x) = 1.0

Δt = 0.01
t_f = 100.0

u₀(x) = sin.(2π*x*2 .+ 1.0)


BoundaryLeft =  FaADE.SATs.SAT_Dirichlet(t->sin(1.0),       Dom.Δx, 1, order)
BoundaryRight = FaADE.SATs.SAT_Dirichlet(t->sin(4π + 1.0),  Dom.Δx, 1, order)

# BoundaryLeftRight = PeriodicBoundary(1)

order = 2
method = :cgie

P = VariableCoefficientPDE1D(u₀,K,order,BoundaryLeft,BoundaryRight)

@benchmark solve($P,$Dom,$Δt,$t_f)

# @profview solve(P,Dom,Δt,t_f,:cgie)
# @profview solve(P,Dom,Δt,t_f,:cgie)

# @time solve(P,Dom,Δt,t_f,:cgie)
# Profile.clear_malloc_data()
# @time solve(P,Dom,Δt,t_f,:cgie)


