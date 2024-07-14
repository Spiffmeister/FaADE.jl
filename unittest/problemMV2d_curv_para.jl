using LinearAlgebra
using Revise
using FaADE
# using BenchmarkTools
# using ProfileView
# using Cthulhu


# Simulation parameters
order = 2
K = 1.0

nx = ny = 21

Δt = 1e-3
# t = Δt
t = 5e-2

ωt = 1.0
ωx = 1.0
ωy = 1.0
cx = 0.0
cy = 0.0

Kx = 1.0
Ky = 1.0

θ = 0.5


# 2D sine wave solution
# Solution
u₀(x,y) = 0.0
F(X,t) = begin
    x,y = X
    4 * (1 - (x^2 + y^2))^8
end
# DIRICHLET
Bxy(X,t) = 0.0   #Boundary condition x=0
# NEUMANN
# BxLũ(y,t) = 2π*ωx * K * cos(2π*t) * cos(cx)             * sin(2π*y*ωy + cy) #Boundary condition x=0



#====== New solver 4 volume ======#
println("Curvilinear volume")

# D1 = Grid2D([-0.25,0.25],[-0.25,0.25],
#             nx,ny,
#             coord=CurvilinearMetric)

D1 = Grid2D(u->[u*0.5 - 0.25,-0.25],
            v->[-0.25,v*0.5 - 0.25],
            v->[0.25,v*0.5 - 0.25],
            u->[u*0.5 - 0.25, 0.25],
    nx,ny)

T = FaADE.Grid.Torus([1.0],[1.0],[1],[0])

# Right domain
D2 = Grid2D(u->[0.25, -u*0.5 + 0.25], # Bottom
            v->v*(T(π/4,0.0) - [0.25,0.25]) + [0.25,0.25], # Left
            v->v*(T(7π/4,0.0) + [-0.25, 0.25]) + [0.25, -0.25], # Right
            u->T(u*(7π/4 - 9π/4) + 9π/4,0.0), # Top
            nx,ny)

# Top domain
D3 = Grid2D(u->[u*0.5 - 0.25, 0.25], # Bottom
            v->v*(T(3π/4,0.0) + [0.25,-0.25]) + [-0.25,0.25], # Left
            v->v*(T(π/4,0.0) - [0.25,0.25]) + [0.25,0.25], # Right
            u->T(u*(π/4 - 3π/4) + 3π/4,0.0), # Top
            nx,ny)

# Left domain
D4 = Grid2D(u->[-0.25,u*0.5 - 0.25],
            v->v*(T(5π/4,0.0) - [-0.25, -0.25]) + [-0.25, -0.25],
            v->v*(T(3π/4,0.0) - [-0.25,0.25]) + [-0.25,0.25],
            u->T(u*(3π/4 - 5π/4) + 5π/4,0.0),
            nx,ny)

# Bottom domain
D5 = Grid2D(u->[-u*0.5 + 0.25, -0.25],
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

Dr = FaADE.SATs.SAT_Dirichlet(Bxy, D2.Δy, Up, order, D2.Δx, :Curvilinear) # Block 2 BCs
Du = FaADE.SATs.SAT_Dirichlet(Bxy, D3.Δy, Up, order, D3.Δx, :Curvilinear) # Block 3 BCs
Dl = FaADE.SATs.SAT_Dirichlet(Bxy, D4.Δy, Up, order, D4.Δx, :Curvilinear) # Block 4 BCs
Dd = FaADE.SATs.SAT_Dirichlet(Bxy, D5.Δy, Up, order, D5.Δx, :Curvilinear) # Block 5 BCs

BD = Dict(2 => (Dr,), 3 => (Du,), 4 => (Dl,), 5 => (Dd,))







δ = 0.005
rs = 0.5
function B(X,x::Array{Float64},params,t)
    X[1] = -x[1] * δ * (-x[1]^4 + 1) * sin(x[2])
    X[2] = -2x[1] + 2rs - 2δ * x[1] * (-x[1]^4 + 1) * cos(x[2]) + 4δ * x[1]^5 * cos(x[2])
    # time dependent
    # X[1] = -x[1] * δ * (-r^4 + 1) * sin(x[2]) * t
    # X[1] = (-2x[1] + 2rs - 2δ * x[1] * (-x[1]^4 + 1) * cos(x[2]) + 4δ * r^5 * cos(x[2])) * t
end
dH(X,x,params,t) = B(X,x,params,t)

XtoB(x,y) = [sqrt(x^2 + y^2), atan(y,x)]
BtoX(r,θ) = [r*cos(θ), r*sin(θ)]
gridoptions = Dict("coords"=>(XtoB,BtoX), "xbound"=>[0.0,1.0], "ybound"=>[0.0,2π], "remapping"=>:bilinear)

gdata = construct_grid(dH,Dom,[-2.0π,2.0π],gridoptions=gridoptions)

# gdata = remap_grid(gdata,interpmode=:idw)

PData = FaADE.ParallelOperator.ParallelMultiBlock(gdata,Dom,order,κ=1.0e6)







P = Problem2D(order,u₀,K,K,Dom,BD,F,PData)
# P = Problem2D(order,u₀,K,K,Dom,BD,F,nothing)

println("---Solving 4 volume---")
# soln = solve(P,Dom,Δt,1.1Δt)
soln = solve(P,Dom,Δt,t)






using GLMakie
gridfig = Figure()
gridfig_ax = Axis(gridfig[1,1])
scatter!(gridfig_ax,D1.gridx[:],D1.gridy[:],markersize=10.5)
scatter!(gridfig_ax,D2.gridx[:],D2.gridy[:],markersize=10.5)
scatter!(gridfig_ax,D3.gridx[:],D3.gridy[:],markersize=10.5)
scatter!(gridfig_ax,D4.gridx[:],D4.gridy[:],markersize=10.5)
scatter!(gridfig_ax,D5.gridx[:],D5.gridy[:],markersize=10.5)

# lines!(gridfig_ax,D5.gridx[:,41],D5.gridy[:,41],color=:black)
# scatter!(gridfig_ax,[-0.10203062474510068],[ -0.9947812581739387])

gridfig






colourrange = (minimum(minimum.(soln.u[2])),maximum(maximum.(soln.u[2])))

f = Figure()
ax = Axis3(f[1,1])
surface!(ax,Dom.Grids[1].gridx, Dom.Grids[1].gridy, soln.u[2][1],colorbar=false, colorrange=colourrange)
surface!(ax,Dom.Grids[2].gridx, Dom.Grids[2].gridy, soln.u[2][2],colorbar=false, colorrange=colourrange)
surface!(ax,Dom.Grids[3].gridx, Dom.Grids[3].gridy, soln.u[2][3],colorbar=false, colorrange=colourrange)
surface!(ax,Dom.Grids[4].gridx, Dom.Grids[4].gridy, soln.u[2][4],colorbar=false, colorrange=colourrange)
surface!(ax,Dom.Grids[5].gridx, Dom.Grids[5].gridy, soln.u[2][5],colorbar=false, colorrange=colourrange)

# scatter!(ax,D1.gridx[:],D1.gridy[:],-ones(length(D1)),markersize=1.5)
# scatter!(ax,D2.gridx[:],D2.gridy[:],-ones(length(D2)),markersize=1.5)
# scatter!(ax,D3.gridx[:],D3.gridy[:],-ones(length(D3)),markersize=1.5)
# scatter!(ax,D4.gridx[:],D4.gridy[:],-ones(length(D4)),markersize=1.5)
# scatter!(ax,D5.gridx[:],D5.gridy[:],-ones(length(D5)),markersize=1.5)
f

