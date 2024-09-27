using Revise

push!(LOAD_PATH,".")
using FaADE

_plot = true


###
𝒟x = [0.0,1.0]
𝒟y = [0.0,1.0]
nx = 41
ny = 41
Dom = Grid2D(𝒟x,𝒟y,nx,ny)

kx = 1.0;
ky = 1.0;

Δt = 0.01
t_f = 100.0

u₀(x,y) = exp(-((x-0.5)^2 + (y-0.5)^2) / 0.02)

order = 2

BoundaryLeft =  FaADE.SATs.SAT_Dirichlet((y,t) -> 0.0, Dom.Δx, Left,    order)
BoundaryRight = FaADE.SATs.SAT_Dirichlet((y,t) -> 0.0, Dom.Δx, Right,   order)
BoundaryUp =    FaADE.SATs.SAT_Dirichlet((x,t) -> 0.0, Dom.Δx, Up,      order)
BoundaryDown =  FaADE.SATs.SAT_Dirichlet((x,t) -> 0.0, Dom.Δx, Down,    order)

BD = (BoundaryLeft,BoundaryRight,BoundaryUp,BoundaryDown)

P = Problem2D(order,u₀,kx,ky,Dom,BD)

soln = solve(P,Dom,Δt,t_f,solver=:theta,θ=1.0)


if _plot
    using GLMakie
    f = Figure()
    ax = Axis3(f[1,1]);
    ax2 = Axis3(f[1,2]);

    surface!(ax,    Dom.gridx,Dom.gridy,soln.u[1])
    surface!(ax2,   Dom.gridx,Dom.gridy,soln.u[2])
end


