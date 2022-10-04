using LinearAlgebra
using Printf
using Plots

using BenchmarkTools
using Profile
using PProf

push!(LOAD_PATH,".")
using SBP_operators




###
𝒟 = [0.0,1.0]
n = 41
Dom = Grid1D(𝒟,n)

K = zeros(Float64,n) .+ 1.0

Δt = 0.05*Dom.Δx^2
t_f = 10.1Δt

# u₀(x) = exp.(-(x.-0.5).^2 ./ 0.02)
u₀(x) = x


BoundaryLeft = Boundary(Dirichlet,t->0.0,Left,1)
BoundaryRight = Boundary(Dirichlet,t->1.0,Right,1)

order = 2
method = :cgie

P = VariableCoefficientPDE1D(u₀,K,order,BoundaryLeft,BoundaryRight)


println(method)
println("Δx=",Dom.Δx,"      ","Δt=",Δt,"        ","final time=",t_f)


###
# @benchmark solve($P,$Dom,$Δt,$t_f,:cgie)


soln = solve(P,Dom,Δt,t_f,:cgie)
plot(soln.grid.grid,soln.u[2],xlims=(0.0,1.0),ylims=(0.0,1.0))

# solve(P,Dom,Δt,t_f,:cgie)
# @pprof solve(P,Dom,Δt,t_f,:cgie)


# @time solve(P,Dom,Δt,t_f,:cgie)
# Profile.clear_malloc_data()
# @time solve(P,Dom,Δt,t_f,:cgie)


# ###
# anim = @animate for i=1:N
#     plot(soln.x,soln.u[i],label="t=$(@sprintf("%.5f",i*Δt))",ylims=(-0.05,1.1))
# end
# gif(anim,"yes.gif",fps=50)


