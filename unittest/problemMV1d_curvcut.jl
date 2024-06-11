using LinearAlgebra
using Revise
using FaADE
# using BenchmarkTools
# using ProfileView
# using Cthulhu


# Simulation parameters
order = 2
K = 1.0

nx = ny = 81

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
D1 = Grid2D(u->u*[cos(7π/4), sin(7π/4)] - [0.0,0.25],
            v->[0.0, v/2 - 0.25],
            v->[cos(v*(9π/4 - 7π/4) + 7π/4), sin(v*(9π/4 - 7π/4) + 7π/4)] - 
                (1-v)*[0.0,0.25] + v*[0.0,0.25],
            u->u*[cos(π/4), sin(π/4)] + [0.0, 0.25], # top of D1
            nx,ny)


Dl = SAT_Dirichlet(Bxy, D1.Δx, Left, order, D1.Δy, :Curvilinear) # Block 3 BCs
Dr = SAT_Dirichlet(Bxy, D1.Δx, Right,order, D1.Δy, :Curvilinear) # 
Dd = SAT_Dirichlet(Bxy, D1.Δy, Down, order, D1.Δx, :Curvilinear) # Block 5 BCs
Du = SAT_Dirichlet(Bxy, D1.Δy, Up,   order, D1.Δx, :Curvilinear) # Block 3 BCs

BD = (Dl,Dr,Dd,Du)



using GLMakie
gridfig = Figure()
gridfix_ax = Axis(gridfig[1,1])
scatter!(gridfix_ax,D1.gridx[:],D1.gridy[:],markersize=2.5)
gridfig



P = Problem2D(order,u₀,K,K,D1,BD,F,nothing)

println("---Solving 4 volume---")
soln = solve(P,D1,Δt,t)

e = zeros(size(D1))
for I in eachindex(D1)
    e[I] = exact(D1.gridx[I],D1.gridy[I],t)
end



f = Figure()

ax1 = Axis3(f[1,1])
surface!(ax1, D1.gridx, D1.gridy, e, colorbar=false)
ax2 = Axis3(f[1,2])
surface!(ax2, D1.gridx, D1.gridy, soln.u[2], colorbar=false)

f



