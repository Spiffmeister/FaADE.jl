using LinearAlgebra
using Printf
using Plots

using BenchmarkTools
# using ProfileView
# using Cthulhu

push!(LOAD_PATH,".")
using FaADE




###
𝒟 = [0.0,1.0]
n = 101
Dom = Grid1D(𝒟,n)

# K = zeros(Float64,n) .+ 1.0
K(x) = 1.0

Δt = 0.01
t_f = 100.0

# u₀(x) = exp.(-(x.-0.5).^2 ./ 0.02)
u₀(x) = sin.(2π*x*2 .+ 1.0)
# u₀(x) = x
# u₀(x) = x.^2


BoundaryLeft =  FaADE.SATs.SAT_Dirichlet(t->sin(1.0),       Dom.Δx, 1, order)
BoundaryRight = FaADE.SATs.SAT_Dirichlet(t->sin(4π + 1.0),  Dom.Δx, 1, order)

# BoundaryLeftRight = PeriodicBoundary(1)

order = 2
method = :cgie

P = VariableCoefficientPDE1D(u₀,K,order,BoundaryLeft,BoundaryRight)
# P = VariableCoefficientPDE1D(u₀,K,order,BoundaryLeftRight)


println(method)
println("Δx=",Dom.Δx,"      ","Δt=",Δt,"        ","final time=",t_f,"   order=",order)


###
# @benchmark solve($P,$Dom,$Δt,$t_f,:cgie)


# soln = solve(P,Dom,Δt,t_f,:cgie)
# scatter(soln.grid.grid,soln.u[1])#,xlims=(0.0,1.0),ylims=(0.0,1.0))
# scatter!(soln.grid.grid,soln.u[2])#,xlims=(0.0,1.0),ylims=(0.0,1.0))


# @profview solve(P,Dom,Δt,t_f,:cgie)
# @profview solve(P,Dom,Δt,t_f,:cgie)
@benchmark solve($P,$Dom,$Δt,$t_f,:cgie)


# @time solve(P,Dom,Δt,t_f,:cgie)
# Profile.clear_malloc_data()
# @time solve(P,Dom,Δt,t_f,:cgie)


# ###
# anim = @animate for i=1:N
#     plot(soln.x,soln.u[i],label="t=$(@sprintf("%.5f",i*Δt))",ylims=(-0.05,1.1))
# end
# gif(anim,"yes.gif",fps=50)


