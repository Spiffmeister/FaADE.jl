"""
Example of a hollow toroidal domain with an island at position specified by xₛ and a width of δ.
"""
using FaADE


order = 2

k_para = 1.0e6
k_perp = 1.0




coord = :Curvilinear

Rin = [3e-1]; Zin=[3e-1]
Rout = [6e-1]; Zout=[6e-1]

inner = FaADE.Grid.Torus(Rin,Zin,[1],[0])
outer = FaADE.Grid.Torus(Rout,Zout,[1],[0])

nx = 41
ny = 101

X,Y = FaADE.Grid.meshgrid(inner,outer,0.0,nx,ny)
Dom = Grid2D(X,Y,periodicy=true)



# Zero everywhere
# u₀(x,y) = 0.0
# S = nothing

# Linear initial condition with max at r=0.3
u₀(x,y) = -((sqrt(x^2 + y^2) - 0.3)/0.3 - (0.6 - sqrt(x^2 + y^2))/0.3) * 0.5 + 0.5
S = nothing

BoundaryLeft    = SAT_Dirichlet((y,t) -> 1.0, Dom.Δx , Left,  order, Dom.Δy, coord)
BoundaryRight   = SAT_Dirichlet((y,t) -> 0.0, Dom.Δx , Right, order, Dom.Δy, coord)
BoundaryUp      = SAT_Periodic(Dom.Δy, Up  , order, Dom.Δx, coord)
BoundaryDown    = SAT_Periodic(Dom.Δy, Down, order, Dom.Δx, coord)


BC = (BoundaryLeft,BoundaryRight,BoundaryUp,BoundaryDown)


δ = 0.1
xₛ = 0.45
function B(X,x,p,t)
    X[1] = δ*(0.3-x[1])*(x[1]-0.6)*sin(x[2])
    X[2] = 2*x[1] - 2*xₛ - δ*(x[1]-0.6)*cos(x[2]) - δ*(x[1]-0.3)*cos(x[2])
end



# Swapping the coordinates back and forth
XtoB(x,y) = [sqrt(x^2 + y^2), atan(y,x)]
BtoX(r,θ) = [r*cos(θ), r*sin(θ)]


intercept(u,x,y,t) = begin
    if (sqrt(x^2 + y^2) ≈ 0.6)
        tmp = zero(typeof(x))
        return tmp
    elseif (sqrt(x^2 + y^2) ≈ 0.3)
        tmp = one(typeof(x))
        return tmp
    else
        return u/2
    end
end


gridoptions = Dict("remap"=>(XtoB,BtoX),"interpmode"=>:chs,"xbound"=>[0.3,0.6],"ybound"=>[-π,π])
interpoptions = Dict("interpolant"=>:chs,"intercept"=>intercept)

gdata   = construct_grid(B,Dom,[-2π,2π],ymode=:stop,remap=(XtoB,BtoX),interpmode=:chs,gridoptions=gridoptions)
PData = ParallelData(gdata,Dom,order,κ=k_para,intercept=intercept,interpolant=:chs,periodicy=true)






# Time setup
Δt = 1e-4
t_f = 1e-2



P = Problem2D(order,u₀,k_perp,k_perp,Dom,BC,parallel=PData)

soln = solve(P, Dom, Δt, 1.1Δt)
soln = solve(P, Dom, Δt, t_f)





# println("Poincare")

# include("../../BADESBP_examples/FieldLines.jl")
# poindata = FieldLines.construct_poincare(B,[0.3,0.6],[0.0,π],N_trajs=200,N_orbs=400)
# poinrtheta = hcat([BtoX(poindata.ψ[I],poindata.θ[I]) for I in eachindex(poindata.ψ)]...)


# println("Plotting")

# using GLMakie
# f = Figure()
# ax = Axis3(f[1,1])
# wireframe!(ax,Dom.gridx,Dom.gridy,soln.u[2])
# scatter!(ax,poinrtheta[1,:],poinrtheta[2,:],zeros(size(poinrtheta)[2]).-1,markersize=2.5,color=:red)


