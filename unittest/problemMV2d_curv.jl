using LinearAlgebra
using Revise
using FaADE
# using BenchmarkTools
# using ProfileView
# using Cthulhu


# Simulation parameters
order = 2
K = 1.0

nx = ny = 11

Δt = 1e-3
t = Δt
# t = 0.76

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
# exact(x,y,t) = cos(2π*ωt*t) * sin(2π*x*ωx + cx) * sin(2π*y*ωy + cy)
# u₀(x,y) = exact(x,y,0.0)
# F(X,t) = begin
#     x,y = X
#     -2π*ωt*sin(2π*ωt*t)*sin(2π*x*ωx + cx)*sin(2π*y*ωy + cy) + 
#     K * 4π^2 * ωx^2 * cos(2π*ωt*t)*sin(2π*x*ωx + cx)*sin(2π*y*ωy + cy) + 
#     K * 4π^2 * ωy^2 * cos(2π*ωt*t)*sin(2π*x*ωx + cx)*sin(2π*y*ωy + cy)
# end
# DIRICHLET
# Bxy(X,t) = cos(2π*ωt*t) * sin(2π*X[1]*ωx + cx) * sin(2π*ωy*X[2] + cy)   #Boundary condition x=0
# NEUMANN
# BxLũ(y,t) = 2π*ωx * K * cos(2π*t) * cos(cx)             * sin(2π*y*ωy + cy) #Boundary condition x=0

exact(R,Z,t) = exp.( -(R.^2 + Z.^2) / 0.1 )
u₀(R,Z) = exact(R,Z,0.0)
F(X,t) = 0.0
Bxy(X,t) = 0.0



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
D2 = Grid2D(u->u*(T(7π/4,0.0) + [-0.25, 0.25]) + [0.25, -0.25],
            v->[0.25, v*0.5 - 0.25],
            v->T(v*(9π/4 - 7π/4) + 7π/4,0.0),
            u->u*(T(π/4,0.0) - [0.25,0.25]) + [0.25,0.25],
            nx,ny)

# Top domain
D3 = Grid2D(u->[u*0.5 - 0.25, 0.25],
            v->v*(T(3π/4,0.0) + [0.25,-0.25]) + [-0.25,0.25],
            v->v*(T(π/4,0.0) - [0.25,0.25]) + [0.25,0.25],
            u->T(u*(π/4 - 3π/4) + 3π/4,0.0),
            nx,ny)

# Left domain
D4 = Grid2D(u->u*(T(5π/4,0.0) + [0.25, 0.25]) + [-0.25,0.25],
            v->T(v*(3π/4 - 5π/4) + 5π/4,0.0),
            v->[-0.25,v*0.5 - 0.25],
            u->u*([0.25,0.25] + T(3π/4,0.0)) + T(3π/4,0.0),
            nx,ny)

# Bottom domain
D5 = Grid2D(u->T(u*(7π/4 - 5π/4) + 5π/4,0.0),
            v->v*([-0.25,-0.25] - T(5π/4,0.0)) + T(5π/4,0.0),
            v->v*([0.25,-0.25] - T(7π/4,0.0)) + T(7π/4,0.0),
            u->[u*0.5 - 0.25, -0.25],
            nx,ny)


joints = ((Joint(2,Right),Joint(3,Up),Joint(4,Left),Joint(5,Down)),
            (Joint(1,Left),Joint(3,Up),Joint(5,Down)),
            (Joint(1,Down),Joint(2,Right),Joint(4,Left)),
            (Joint(1,Right),Joint(3,Up),Joint(5,Down)),
            (Joint(1,Up),Joint(4,Left),Joint(2,Right)))

# joints = (1=>((2,Right),(3,Up),(4,Left),(5,Down)),
#           2=>((1,Left),(3,Up),(5,Down)),
#           3=>((1,Down),(2,Right),(4,Left)),
#           4=>((1,Right),(3,Up),(5,Down)),
#           5=>((1,Up),(4,Left),(2,Right)))

Dom = GridMultiBlock((D1,D2,D3,D4,D5),joints)

Dr = FaADE.SATs.SAT_Dirichlet(Bxy, D2.Δy, Right,order, D2.Δx, :Curvilinear) # Block 2 BCs
Du = FaADE.SATs.SAT_Dirichlet(Bxy, D3.Δx, Up,   order, D2.Δy, :Curvilinear) # Block 3 BCs
Dl = FaADE.SATs.SAT_Dirichlet(Bxy, D4.Δy, Left, order, D2.Δx, :Curvilinear) # Block 4 BCs
Dd = FaADE.SATs.SAT_Dirichlet(Bxy, D5.Δx, Down, order, D2.Δy, :Curvilinear) # Block 5 BCs

BD = Dict(2 => (Dr,), 3 => (Du,), 4 => (Dl,), 5 => (Dd,))



using GLMakie
gridfig = Figure()
gridfix_ax = Axis(gridfig[1,1])
scatter!(gridfix_ax,D1.gridx[:],D1.gridy[:],markersize=10.5)
scatter!(gridfix_ax,D2.gridx[:],D2.gridy[:],markersize=10.5)
scatter!(gridfix_ax,D3.gridx[:],D3.gridy[:],markersize=10.5)
scatter!(gridfix_ax,D4.gridx[:],D4.gridy[:],markersize=10.5)
scatter!(gridfix_ax,D5.gridx[:],D5.gridy[:],markersize=10.5)
gridfig



P = Problem2D(order,u₀,K,K,Dom,BD,F,nothing)

println("---Solving 4 volume---")
soln = solve(P,Dom,Δt,t)


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

