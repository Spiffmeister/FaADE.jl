
using Revise
using FaADE


θ = 0.5
order = 2

k_para = 1.0e6
k_perp = 1.0

# u₀(x,y) = -((sqrt(x^2 + y^2) - 0.3)/0.3 - (0.6 - sqrt(x^2 + y^2))/0.3) * 0.5 + 0.5
u₀(x,y) = 0.0
S = nothing




δ = 0.01
xₛ = 0.45

function B(X,x,p,t)
    # bn = 1 + abs( δ*x[1]*(x[1]-1)*sin(x[2]) )^2 + abs( 2*x[1] - 2*xₛ + δ*(1-x[1])*cos(x[2]) - δ*x[1]*cos(x[2]) )^2
    # bn = sqrt(bn)
    X[1] = δ*(0.3-x[1])*(x[1]-0.6)*sin(x[2])
    X[2] = 2*x[1] - 2*xₛ - δ*(x[1]-0.6)*cos(x[2]) - δ*(x[1]-0.3)*cos(x[2])
    # X[3] = 0.0
end






coord = :Curvilinear

Rin = [3e-1]; Zin=[3e-1]
Rout = [6e-1]; Zout=[6e-1]

inner = FaADE.Grid.Torus(Rin,Zin,[1],[0])
outer = FaADE.Grid.Torus(Rout,Zout,[1],[0])

nx = ny = 21

X,Y = FaADE.Grid.meshgrid(inner,outer,0.0,nx,ny)
Dom = Grid2D(X,Y,periodicy=true)



BoundaryLeft    = SAT_Dirichlet((y,t) -> 1.0, Dom.Δx , Left,  order, Dom.Δy, coord)
BoundaryRight   = SAT_Dirichlet((y,t) -> 0.0, Dom.Δx , Right, order, Dom.Δy, coord)
BoundaryUp      = SAT_Periodic(Dom.Δy, Up  , order, Dom.Δx, coord)
BoundaryDown    = SAT_Periodic(Dom.Δy, Down, order, Dom.Δx, coord)

BC = (BoundaryLeft,BoundaryRight,BoundaryUp,BoundaryDown)




# Time setup
Δt = 1e-4
t_f = 1e-2
nf = round(t_f/Δt)
Δt = t_f/nf





# Swapping the coordinates back and forth
XtoB(x,y) = [sqrt(x^2 + y^2), atan(y,x)]
BtoX(r,θ) = [r*cos(θ), r*sin(θ)]


gridoptions = Dict("remap"=>(XtoB,BtoX),"interpmode"=>:chs,"xbound"=>[0.3,0.6],"ybound"=>[-π,π])

gdata   = construct_grid(B,Dom,[-2π,2π],ymode=:stop,remap=(XtoB,BtoX),interpmode=:chs,gridoptions=gridoptions)

periodicy=true

PData = ParallelData(gdata,Dom,order,κ=k_para,interpolant=:chs,periodicy=true)





P = Problem2D(order,u₀,k_perp,k_perp,Dom,BC,parallel=PData)

soln = solve(P, Dom, Δt, 1.1Δt, solver=:theta, θ=θ)




