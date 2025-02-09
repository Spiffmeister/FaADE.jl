
using Revise
using FaADE


θ = 0.5
order = 2

k_para = 1.0e6
k_perp = 1.0


nx = ny = 41



# Time setup
Δt = 1e-4
t_f = 1e-2
nf = round(t_f/Δt)
Δt = t_f/nf





coord = :Curvilinear

println("Curvilinear volume")


dialation = 0.2


D1 = Grid2D(
    u -> [-0.25,-0.25] + u*([0.25,-0.25] - [-0.25,-0.25]) + u*(1-u)*[0.0,-dialation],
    v -> [-0.25,-0.25] + v*([-0.25,0.25] - [-0.25,-0.25]) + v*(1-v)*[-dialation,0.0],
    v -> [0.25,-0.25] + v*([0.25,0.25] - [0.25,-0.25]) + v*(1-v)*[dialation,0.0],
    u -> [-0.25,0.25] + u*([0.25,0.25] - [-0.25,0.25]) + u*(1-u)*[0.0,dialation],
    nx,ny
)



T = FaADE.Grid.Torus([1.0],[1.0],[1],[0])

# Right domain
D2 = Grid2D(#u->[0.25, -u*0.5 + 0.25], # Bottom
            u->[0.25,0.25] + u*([0.25,-0.25] - [0.25,0.25]) + u*(1-u)*[dialation,0.0],
            v->v*(T(π/4,0.0) - [0.25,0.25]) + [0.25,0.25], # Left
            v->v*(T(7π/4,0.0) + [-0.25, 0.25]) + [0.25, -0.25], # Right
            u->T(u*(7π/4 - 9π/4) + 9π/4,0.0), # Top
            nx,ny)

# Top domain
D3 = Grid2D(#u->[u*0.5 - 0.25, 0.25], # Bottom
            u->[-0.25,0.25] + u*([0.25,0.25] - [-0.25,0.25]) + u*(1-u)*[0.0,dialation],
            v->v*(T(3π/4,0.0) + [0.25,-0.25]) + [-0.25,0.25], # Left
            v->v*(T(π/4,0.0) - [0.25,0.25]) + [0.25,0.25], # Right
            u->T(u*(π/4 - 3π/4) + 3π/4,0.0), # Top
            nx,ny)

# Left domain
D4 = Grid2D(#u->[-0.25,u*0.5 - 0.25],
            u->[-0.25,-0.25] + u*([-0.25,0.25] - [-0.25,-0.25]) + u*(1-u)*[-dialation,0.0],
            v->v*(T(5π/4,0.0) - [-0.25, -0.25]) + [-0.25, -0.25],
            v->v*(T(3π/4,0.0) - [-0.25,0.25]) + [-0.25,0.25],
            u->T(u*(3π/4 - 5π/4) + 5π/4,0.0),
            nx,ny)

# Bottom domain
D5 = Grid2D(#u->[-u*0.5 + 0.25, -0.25],
            u->[-0.25,-0.25] + u*([0.25,-0.25] - [-0.25,-0.25]) + u*(1-u)*[0.0,-dialation],
            v->v*(T(7π/4,0.0) - [0.25,-0.25]) + [0.25, -0.25],
            v->v*(T(5π/4,0.0) - [-0.25,-0.25]) + [-0.25, -0.25],
            u->T(u*(5π/4 - 7π/4) + 7π/4, 0.0),
            nx,ny)


joints = ((Joint(2,Right),Joint(3,Up),Joint(4,Left),Joint(5,Down)),
            (Joint(1,Down),Joint(3,Left),Joint(5,Right)),
            (Joint(1,Down),Joint(4,Left),Joint(2,Right)),
            (Joint(1,Down),Joint(5,Left),Joint(3,Right)),
            (Joint(1,Down),Joint(2,Left),Joint(4,Right)))


Dom = GridMultiBlock((D1,D2,D3,D4,D5),joints)




# Zero everywhere
# u₀(x,y) = 0.0
# S = nothing

u₀(x,y) = 0.0
# S = nothing

S(X,t) = begin # with forcing term from Chacón paper
    x,y = X
    4 * (1 - (x^2 + y^2))^8
end

Bxy(X,t) = 0.0


Dr = SAT_Dirichlet(Bxy, D2.Δy, Up, order, D2.Δx, :Curvilinear) # Block 2 BCs
Du = SAT_Dirichlet(Bxy, D3.Δy, Up, order, D3.Δx, :Curvilinear) # Block 3 BCs
Dl = SAT_Dirichlet(Bxy, D4.Δy, Up, order, D4.Δx, :Curvilinear) # Block 4 BCs
Dd = SAT_Dirichlet(Bxy, D5.Δy, Up, order, D5.Δx, :Curvilinear) # Block 5 BCs

BD = Dict(2 => (Dr,), 3 => (Du,), 4 => (Dl,), 5 => (Dd,))








# Magnetic field

δ = 0.05
rs = 0.7
function B(X,x::Array{Float64},params,t)
    X[1] = δ*x[1]*(1-x[1])*sin(x[2])#/bn
    X[2] = (2x[1] - 2*rs + δ*(1-x[1])*cos(x[2]) - δ*x[1]*cos(x[2]))#/bn
end
dH(X,x,params,t) = B(X,x,params,t)

XtoB(x,y) = [sqrt(x^2 + y^2), atan(y,x)]
BtoX(r,θ) = [r*cos(θ), r*sin(θ)]



intercept(u,x,y,t) = begin
    if (sqrt(x^2 + y^2) ≈ 1.0)
        tmp = zero(typeof(x))
        return tmp
    else
        return u
    end
end


gridoptions = Dict("coords"=>(XtoB,BtoX), "xbound"=>[0.0,1.0], "ybound"=>[0.0,2π])
interpoptions = Dict("interpolant"=>:chs,"intercept"=>intercept)

gdata   = construct_grid(B,Dom,[-2π,2π],gridoptions=gridoptions)

PData = FaADE.ParallelOperator.ParallelMultiBlock(gdata,Dom,order,κ=k_para,interpopts=interpoptions)



println("Solve")

P = Problem2D(order,u₀,k_perp,k_perp,Dom,BD,parallel=PData,source=S)

soln = solve(P, Dom, Δt, 1.1Δt)
soln = solve(P, Dom, Δt, t_f)





println("Poincare")

include("../../BADESBP_examples/FieldLines.jl")
poindata = FieldLines.construct_poincare(B,[0.0,1.0],[0.0,π],N_trajs=200,N_orbs=400)
poinrtheta = hcat([BtoX(poindata.ψ[I],poindata.θ[I]) for I in eachindex(poindata.ψ)]...)


println("Plotting")

using GLMakie
f = Figure()
ax = Axis3(f[1,1])
for I in eachindex(Dom.Grids)
    wireframe!(ax,Dom.Grids[I].gridx,Dom.Grids[I].gridy,soln.u[2][I])
end
scatter!(ax,poinrtheta[1,:],poinrtheta[2,:],zeros(size(poinrtheta)[2]),markersize=2.5,color=:red)





hind = Int(floor(nx/2)+1)

h = Figure()
axh = Axis(h[1,1])
lines!(axh,Dom.Grids[4].gridx[hind,:],soln.u[2][5][hind,:])
lines!(axh,Dom.Grids[1].gridx[:,hind],soln.u[2][1][:,hind])



#=


using GLMakie

g = Figure()
axg = Axis(g[1,1])
for I in eachindex(Dom.Grids)
    wireframe!(axg,Dom.Grids[I].gridx,Dom.Grids[I].gridy,soln.u[2][I])
end
scatter!(axg,[0.7071067811865477],[0.7071067811865476],color=:red)
g


=#


