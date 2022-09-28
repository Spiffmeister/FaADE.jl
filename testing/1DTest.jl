using LinearAlgebra
using Printf
using Plots

using BenchmarkTools
push!(LOAD_PATH,".")
using SBP_operators




###
𝒟 = [0.0,1.0]
n = 41
Dom = Grid1D(𝒟,n)

K = zeros(Float64,n) .+ 1.0

Δt = 0.05*Dom.Δx^2
t_f = 1000Δt

u₀(x) = exp.(-(x.-0.5).^2 ./ 0.02)

g₀(t) = 0.0
g₁(t) = 0.0

BoundaryLeft = Boundary(Neumann,t->0.0,Left,1)
BoundaryRight = Boundary(Neumann,t->0.0,Right,1)

order = 2
method = :cgie

P = VariableCoefficientPDE1D(u₀,K,order,BoundaryLeft,BoundaryRight)


println(method)
println("Δx=",Dom.Δx,"      ","Δt=",Δt,"        ","final time=",t_f)


###
# @benchmark solve($P,$Dom,$Δt,$t_f,:cgie)

soln = solve(P,Dom,Δt,t_f,:cgie)
plot(soln.grid.grid,soln.u[2],xlims=(0.0,1.0),ylims=(0.0,1.0))

# println("Plotting")

# N = length(soln.t)

# ###
# anim = @animate for i=1:N
#     plot(soln.x,soln.u[i],label="t=$(@sprintf("%.5f",i*Δt))",ylims=(-0.05,1.1))
# end
# gif(anim,"yes.gif",fps=50)



# @benchmark SBP_operators.time_solver(rate,u₀,n,x,Δx,t_f,Δt,k,g,Dirichlet,method=method,order=order)
