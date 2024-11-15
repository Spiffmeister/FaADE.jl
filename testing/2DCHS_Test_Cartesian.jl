
using Revise
using FaADE


θ = 0.5
order = 2

k_para = 1.0e6
k_perp = 1.0




nx = ny = 41

Dom = Grid2D([0.0,1.0],[0.0,2π],nx,ny)







# Zero everywhere
# u₀(x,y) = 0.0
# S = nothing

# Linear initial condition with max at r=0.3
u₀(x,y) = -((sqrt(x^2 + y^2) - 0.3)/0.3 - (0.6 - sqrt(x^2 + y^2))/0.3) * 0.5 + 0.5
S = nothing

BoundaryLeft    = SAT_Dirichlet((y,t) -> 1.0, Dom.Δx , Left,  order, Dom.Δy, coord)
BoundaryRight   = SAT_Dirichlet((y,t) -> 0.0, Dom.Δx , Right, order, Dom.Δy, coord)



# Wave everywhere
# u₀(x,y) = sin(2π*x) * sin(2π*y)
# F(x,y,t) = -2π*sin(2π*t)*sin(2π*x)*sin(2π*y) + 
#             4π^2 * cos(2π*t)*sin(2π*x)*sin(2π*y) + 
#             4π^2 * cos(2π*t)*sin(2π*x)*sin(2π*y)

# BxL(X,t) = cos(2π*t) * sin(2π*X[1]) * sin(2π*X[2])
# BxR(X,t) = cos(2π*t) * sin(2π*X[1]) * sin(2π*X[2])

# BoundaryLeft    = SAT_Dirichlet(BxL, Dom.Δx, Left,  order, Dom.Δy, coord)
# BoundaryRight   = SAT_Dirichlet(BxR, Dom.Δx, Right, order, Dom.Δy, coord)




BoundaryUp      = SAT_Periodic(Dom.Δy, Up  , order, Dom.Δx, coord)
BoundaryDown    = SAT_Periodic(Dom.Δy, Down, order, Dom.Δx, coord)



δ = 0.01
xₛ = 0.45

function B(X,x,p,t)
    # bn = 1 + abs( δ*x[1]*(x[1]-1)*sin(x[2]) )^2 + abs( 2*x[1] - 2*xₛ + δ*(1-x[1])*cos(x[2]) - δ*x[1]*cos(x[2]) )^2
    # bn = sqrt(bn)
    X[1] = δ*(0.3-x[1])*(x[1]-0.6)*sin(x[2])
    X[2] = 2*x[1] - 2*xₛ - δ*(x[1]-0.6)*cos(x[2]) - δ*(x[1]-0.3)*cos(x[2])
    # X[3] = 0.0
end







BC = (BoundaryLeft,BoundaryRight,BoundaryUp,BoundaryDown)




# Time setup
Δt = 1e-4
t_f = 1e-2
nf = round(t_f/Δt)
Δt = t_f/nf





# Swapping the coordinates back and forth
XtoB(x,y) = [sqrt(x^2 + y^2), atan(y,x)]
BtoX(r,θ) = [r*cos(θ), r*sin(θ)]


intercept(x,y,t) = begin
    if (sqrt(x^2 + y^2) ≈ 0.6)
        tmp = zero(typeof(x))
        # tmp = BxL((x,y),t)
        return tmp
    else
        return NaN
    end
end


gridoptions = Dict("remap"=>(XtoB,BtoX),"interpmode"=>:chs,"xbound"=>[0.3,0.6],"ybound"=>[-π,π])
interpoptions = Dict("interpolant"=>:chs,"intercept"=>intercept)

gdata   = construct_grid(B,Dom,[-2π,2π],ymode=:stop,remap=(XtoB,BtoX),interpmode=:chs,gridoptions=gridoptions)

# PData = ParallelData(gdata,Dom,order,κ=k_para,interpoptions=interpoptions)
PData = ParallelData(gdata,Dom,order,κ=k_para,intercept=intercept,interpolant=:chs,periodicy=true)


# construct_parallel(B,Dom,[-2π,2π],order,gridoptions=gridoptions,interpoptions=interpoptions)



P = Problem2D(order,u₀,k_perp,k_perp,Dom,BC,parallel=PData)

soln = solve(P, Dom, Δt, 1.1Δt)
soln = solve(P, Dom, Δt, t_f)



using GLMakie
f = Figure()
ax = Axis3(f[1,1])
surface!(ax,Dom.gridx,Dom.gridy,soln.u[2])
