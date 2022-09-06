using LinearAlgebra
using Printf
using Plots

using BenchmarkTools
push!(LOAD_PATH,".")
using SBP_operators




###
𝒟 = [0.0,1.0]
n = 51
order = 2

grid = Grid1D(𝒟,n,order)



k = zeros(Float64,n) .+ 1.0

Δt = 1.0 * Δx^2
t_f = 1000Δt
# t_f = 100.0

u₀(x) = exp.(-(x.-0.5).^2 ./ 0.02)

g(t) = [1.0, 0.0]

order = 2
method = :cgie
println(method)
println("Δx=",Δx,"      ","Δt=",Δt,"        ","final time=",t_f)


###
soln = SBP_operators.time_solver(rate,u₀,n,x,Δx,t_f,Δt,k,g,Dirichlet,method=method,order=order)

println("Plotting")

N = length(soln.t)

###
anim = @animate for i=1:N
    plot(soln.x,soln.u[i],label="t=$(@sprintf("%.5f",i*Δt))",ylims=(-0.05,1.1))
end
gif(anim,"yes.gif",fps=50)



@benchmark SBP_operators.time_solver(rate,u₀,n,x,Δx,t_f,Δt,k,g,Dirichlet,method=method,order=order)
