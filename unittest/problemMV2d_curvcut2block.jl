using LinearAlgebra
using Revise
using FaADE
# using BenchmarkTools
# using ProfileView
# using Cthulhu


# Simulation parameters
order = 2
K = 1.0

nx = ny = 41

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

# exact(R,Z,t) = exp.( -(R.^2 + Z.^2) / 0.1 )
exact(R,Z,t) = sin(2π*R*ωx + cx)
u₀(R,Z) = exact(R,Z,0.0)
F(X,t) = -4π^2*sin(2π)
# Bxy(X,t) = exp.( -(X[1].^2 + X[2].^2) / 0.1 )
Bxy(X,t) = sin(2π*X[1]*ωx + cx)



D1 = Grid2D([-0.25,0.25],[-0.25,0.25],nx,ny)
D2 = Grid2D([0.25,0.75],[-0.25,0.25],nx,ny)

joints = ((Joint(2,Right),),
            (Joint(1,Left),),)

Domcart = GridMultiBlock((D1,D2),joints)

Dd1 = SAT_Dirichlet(Bxy, D1.Δy, Down, order) # Block 5 BCs
Du1 = SAT_Dirichlet(Bxy, D1.Δy, Up,   order) # Block 3 BCs
Dl1 = SAT_Dirichlet((X,t)->-1.0, D1.Δx, Left, order) # Block 1 BCs

Dr2 = SAT_Dirichlet((X,t)->-1.0, D2.Δy, Right,order) # Block 2 BCs
Dd2 = SAT_Dirichlet(Bxy, D2.Δx, Down, order) # Block 5 BCs
Du2 = SAT_Dirichlet(Bxy, D2.Δx, Up,   order) # Block 3 BCs

BD = Dict(1 => (Dd1,Du1,Dl1), 2 => (Dr2,Dd2,Du2))


Pcart = Problem2D(order,u₀,K,K,Domcart,BD,F,nothing)

solncart = solve(Pcart,Domcart,Δt,t)



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



joints = ((Joint(2,Right),),
            (Joint(1,Left),),
            )


Dom = GridMultiBlock((D1,D2),joints)

Dd1 = FaADE.SATs.SAT_Dirichlet(Bxy, D1.Δy, Down, order, D1.Δx, :Curvilinear) # Block 5 BCs
Du1 = FaADE.SATs.SAT_Dirichlet(Bxy, D1.Δy, Up,   order, D1.Δx, :Curvilinear) # Block 3 BCs
Dl1 = SAT_Dirichlet((X,t)->-1.0, D1.Δx, Left, order, D1.Δy, :Curvilinear) # Block 1 BCs

Dr2 = FaADE.SATs.SAT_Dirichlet(Bxy, D2.Δy, Right,order, D2.Δx, :Curvilinear) # Block 2 BCs
Dd2 = FaADE.SATs.SAT_Dirichlet(Bxy, D2.Δx, Down, order, D2.Δy, :Curvilinear) # Block 5 BCs
Du2 = FaADE.SATs.SAT_Dirichlet(Bxy, D2.Δx, Up,   order, D2.Δy, :Curvilinear) # Block 3 BCs


BD = Dict(1 => (Dd1,Du1,Dl1), 2 => (Dr2,Dd2,Du2))



using GLMakie
gridfig = Figure()
gridfix_ax = Axis(gridfig[1,1])
scatter!(gridfix_ax,D1.gridx[:],D1.gridy[:],markersize=1.5)
scatter!(gridfix_ax,D2.gridx[:],D2.gridy[:],markersize=1.5)
gridfig



P = Problem2D(order,u₀,K,K,Dom,BD,F,nothing)

println("---Solving 4 volume---")
soln = solve(P,Dom,Δt,t)


colourrange = (minimum(minimum.(soln.u[2])),maximum(maximum.(soln.u[2])))

f = Figure()

ax1 = Axis3(f[1,1])
surface!(ax1,Domcart.Grids[1].gridx, Domcart.Grids[1].gridy, solncart.u[2][1],colorbar=false, colorrange=colourrange)
surface!(ax1,Domcart.Grids[2].gridx, Domcart.Grids[2].gridy, solncart.u[2][2],colorbar=false, colorrange=colourrange)

# scatter!(ax1,D1.gridx[:],D1.gridy[:],-ones(length(D1)),markersize=1.5)
# scatter!(ax1,D2.gridx[:],D2.gridy[:],-ones(length(D2)),markersize=1.5)
# scatter!(ax1,D3.gridx[:],D3.gridy[:],-ones(length(D3)),markersize=1.5)

ax2 = Axis3(f[1,2])
surface!(ax2,Dom.Grids[1].gridx, Dom.Grids[1].gridy, soln.u[2][1],colorbar=false, colorrange=colourrange)
surface!(ax2,Dom.Grids[2].gridx, Dom.Grids[2].gridy, soln.u[2][2],colorbar=false, colorrange=colourrange)

scatter!(ax2,D1.gridx[:],D1.gridy[:],-ones(length(D1)),markersize=1.5)
scatter!(ax2,D2.gridx[:],D2.gridy[:],-ones(length(D2)),markersize=1.5)

f




h = Figure()
axh1 = Axis(h[1,1])
# lines!(axh1,Domcart.Grids[1].gridx[:,21],solncart.u[2][1][:,21])
# lines!(axh1,Domcart.Grids[2].gridx[:,21],solncart.u[2][2][:,21])
lines!(axh1,Domcart.Grids[1].gridy[1,:],solncart.u[2][1][end,:])
scatter!(axh1,Domcart.Grids[2].gridy[1,:],solncart.u[2][2][1,:])
axh2 = Axis(h[1,2])
# lines!(axh2,Dom.Grids[1].gridx[:,21],soln.u[2][1][:,21])
# lines!(axh2,Dom.Grids[2].gridx[:,21],soln.u[2][2][:,21])
lines!(axh2,Dom.Grids[1].gridy[1,:],soln.u[2][1][end,:])
scatter!(axh2,Dom.Grids[2].gridy[1,:],soln.u[2][2][1,:])
