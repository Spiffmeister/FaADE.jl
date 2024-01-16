using Revise

push!(LOAD_PATH,".")
using FaADE

_plot = true


###
𝒟 = [0.0,1.0]
n = 41
Dom = Grid1D(𝒟,n)


K = 1.0

Δt = 0.01
t_f = 100.0

u₀(x) = exp.(-(x.-0.5).^2 ./ 0.02)

order = 2

BoundaryLeft = FaADE.SATs.SAT_Periodic(Dom.Δx,1,order,Left)
BoundaryRight = FaADE.SATs.SAT_Periodic(Dom.Δx,1,order,Right)

BD = FaADE.Inputs.SATBoundaries(BoundaryLeft,BoundaryRight)

P = newProblem1D(order,u₀,K,Dom,BD)


println("Δx=",Dom.Δx,"      ","Δt=",Δt,"        ","final time=",t_f)



soln = solve(P,Dom,Δt,t_f,solver=:theta,θ=1.0)


if _plot
    using Plots
    plot(Dom.grid,soln.u[1])
    plot!(Dom.grid,soln.u[2])
end

