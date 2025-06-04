
using Test
using FaADE

order = 2

nx = ny = 21

Dom = Grid2D([0.0,1.0],[0.0,2π],nx,ny)

k_perp = 1.0
k_para = 1.0e9

Δt = 0.01
t_f = 1.0

u₀(x,y) = x

BoundaryLeft    = SAT_Dirichlet((y,t) -> 0.0, Dom.Δx , Left,  order)
BoundaryRight   = SAT_Dirichlet((y,t) -> 1.0, Dom.Δx , Right, order)
BoundaryUp      = SAT_Periodic(Dom.Δy, Up  , order)
BoundaryDown    = SAT_Periodic(Dom.Δy, Down, order)

BC = (BoundaryLeft, BoundaryRight, BoundaryUp, BoundaryDown)


function B(X,x,p,t)
    X[1] = 0.0
    X[2] = 0.0
end


gridoptions = Dict("interpmode"=>:chs,"xbound"=>[0.0,1.0],"ybound"=>[0.0,2π])
interpoptions = Dict("interpolant"=>:chs,"intercept"=>nothing)

gdata   = construct_grid(B,Dom,[-2π,2π],ymode=:period,interpmode=:chs,gridoptions=gridoptions)
PData = ParallelData(gdata,Dom,order,κ=k_para,interpolant=:chs)


P = Problem2D(order,u₀,k_perp,k_perp,Dom,BC,parallel=PData)


soln = solve(P,Dom,Δt,t_f)


test_solution = repeat(collect(0.0:Dom.Δx:1.0),1,21)


@test all(soln.u[2] .== test_solution)



