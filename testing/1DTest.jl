using Revise

push!(LOAD_PATH,".")
using FaADE


_plot = false

###
𝒟 = [0.0,1.0]
n = 101
Dom = Grid1D(𝒟,n)

K(x) = 1.0

Δt = 0.01
t_f = 100.0

u₀(x) = sin.(2π*x*2 .+ 1.0)

order = 2

BoundaryLeft =  SAT_Dirichlet(t->sin(1.0),       Dom.Δx, Left, order)
BoundaryRight = SAT_Dirichlet(t->sin(4π + 1.0),  Dom.Δx, Right, order)

BD = SATBoundaries(BoundaryLeft,BoundaryRight)


P = Problem1D(order,u₀,K,Dom,BD)

println("Δx=",Dom.Δx,"      ","Δt=",Δt,"        ","final time=",t_f,"   order=",order)

soln = solve(P,Dom,Δt,t_f,solver=:theta,θ=1.0)



if _plot
    using Plots
    plot(Dom.grid,soln.u[1])
    plot!(Dom.grid,soln.u[2])
end


