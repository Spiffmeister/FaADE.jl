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
# t = Δt
t = 0.76

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
F(X,t) = 4π^2*ωx^2*sin(2π*X[1]*ωx + cx)
Bxy(X,t) = sin(2π*X[1]*ωx + cx)


#====== New solver 4 volume ======#
println("Curvilinear volume")

# Right domain
Dx = Grid2D(u->u*[0.5, 0.0] - [0.0,0.0],
            v->v*[cos(3π/4),sin(3π/4)],
            v->v*[cos(π/4),sin(π/4)] + [0.5,0.0],
            u->[cos(u*(π/4 - 3π/4) + 3π/4), sin(u*(π/4 - 3π/4) + 3π/4)] + 
                u*[0.5,0.0], # top of D1
            nx,ny)

Dy = Grid2D(u->u*[cos(7π/4), sin(7π/4)] - [0.0,0.25],
            v->[0.0, v/2 - 0.25],
            v->[cos(v*(9π/4 - 7π/4) + 7π/4), sin(v*(9π/4 - 7π/4) + 7π/4)] - 
                (1-v)*[0.0,0.25] + v*[0.0,0.25],
            u->u*[cos(π/4), sin(π/4)] + [0.0, 0.25], # top of D1
            nx,ny)



using GLMakie
gridfig = Figure()
gridfig_axx = Axis(gridfig[1,1])
scatter!(gridfig_axx,Dx.gridx[:],Dx.gridy[:],markersize=2.5)
gridfig_axy = Axis(gridfig[1,2])
scatter!(gridfig_axy,Dy.gridx[:],Dy.gridy[:],markersize=2.5)

gridfig
            



println("---Solving x orientation---")
Dlx = SAT_Dirichlet(Bxy, Dx.Δx, FaADE.Left,   order, Dx.Δy, :Curvilinear) #
Drx = SAT_Dirichlet(Bxy, Dx.Δx, FaADE.Right,  order, Dx.Δy, :Curvilinear) #
Ddx = SAT_Dirichlet(Bxy, Dx.Δy, FaADE.Down,   order, Dx.Δx, :Curvilinear) #
Dux = SAT_Dirichlet(Bxy, Dx.Δy, FaADE.Up,     order, Dx.Δx, :Curvilinear) #

BDx = (Dlx,Drx,Ddx,Dux)

Px = Problem2D(order,u₀,K,K,Dx,BDx,F,nothing)

solnx = solve(Px,Dx,Δt,t)




println("---Solving y orientation---")
Dly = SAT_Dirichlet(Bxy, Dy.Δx, FaADE.Left,   order, Dy.Δy, :Curvilinear) #
Dry = SAT_Dirichlet(Bxy, Dy.Δx, FaADE.Right,  order, Dy.Δy, :Curvilinear) #
Ddy = SAT_Dirichlet(Bxy, Dy.Δy, FaADE.Down,   order, Dy.Δx, :Curvilinear) #
Duy = SAT_Dirichlet(Bxy, Dy.Δy, FaADE.Up,     order, Dy.Δx, :Curvilinear) #

BDy = (Dly,Dry,Ddy,Duy)

Py = Problem2D(order,u₀,K,K,Dy,BDy,F,nothing)

solny = solve(Py,Dy,Δt,t)





f = Figure()

ax1 = Axis3(f[1,1])
surface!(ax1,Dx.gridx, Dx.gridy, solnx.u[2],colorbar=false)
ax2 = Axis3(f[1,2])
surface!(ax2,Dy.gridx, Dy.gridy, solny.u[2],colorbar=false)

f




