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
t = 0.76
# t = 10Δt

ωt = 1.0
ωx = 1.0
ωy = 1.0
cx = 1.0
cy = 0.5

Kx = 1.0
Ky = 1.0

θ = 0.5


# 2D sine wave solution
# Solution
exact(x,y,t) = cos(2π*ωt*t) * sin(2π*x*ωx + cx) * sin(2π*y*ωy + cy)
u₀(x,y) = exact(x,y,0.0)
F(X,t) = begin
    x,y = X
    -2π*ωt*sin(2π*ωt*t)*sin(2π*x*ωx + cx)*sin(2π*y*ωy + cy) + 
    K * 4π^2 * ωx^2 * cos(2π*ωt*t)*sin(2π*x*ωx + cx)*sin(2π*y*ωy + cy) + 
    K * 4π^2 * ωy^2 * cos(2π*ωt*t)*sin(2π*x*ωx + cx)*sin(2π*y*ωy + cy)
end

# DIRICHLET
Bx0(y,t) = cos(2π*ωt*t) * sin(cx)               * sin(2π*ωy*y + cy)     #Boundary condition x=0
Bxl(y,t) = cos(2π*ωt*t) * sin(2π*1.5*ωx + cx)   * sin(2π*ωy*y + cy)     #Boundary condition x=Lx=1.5
By0(x,t) = cos(2π*ωt*t) * sin(2π*ωx*x + cx)     * sin(cy)               #Boundary condition y=0
Byl(x,t) = cos(2π*ωt*t) * sin(2π*ωx*x + cx)     * sin(2π*1.5*ωy + cy)   #Boundary condition y=Ly=1.5

Bxh(y,t) = cos(2π*ωt*t) * sin(2π*ωx + cx)       * sin(2π*y*ωy + cy) #Boundary condition x=1.0
Byh(x,t) = cos(2π*ωt*t) * sin(2π*x*ωx + cx)     * sin(2π*ωy + cy)   #Boundary condition y=1.0

#====== 8 Block wedge ======#
#   ---------
#   | 7 | 8 |
#   -------------
#   | 4 | 5 | 6 |
#   -------------
#   | 1 | 2 | 3 |
#   -------------   
#

# Bottom row
D1 = Grid2D([0.0,0.5],[0.0,0.5],21,21)
D2 = Grid2D([0.5,1.0],[0.0,0.5],21,21)
D3 = Grid2D([1.0,1.5],[0.0,0.5],21,21)
# Second row
D4 = Grid2D([0.0,0.5],[0.5,1.0],21,21)
D5 = Grid2D([0.5,1.0],[0.5,1.0],21,21)
D6 = Grid2D([1.0,1.5],[0.5,1.0],21,21)
# Third row
D7 = Grid2D([0.0,0.5],[1.0,1.5],21,21)
D8 = Grid2D([0.5,1.0],[1.0,1.5],21,21)


joints = ((Joint(2,Right),Joint(4,Up)),                                 # Block 1
            (Joint(3,Right),Joint(5,Up),Joint(1,Left)),                 # Block 2
            (Joint(6,Up),Joint(2,Left)),                                # Block 3
            (Joint(5,Right),Joint(7,Up),Joint(1,Down)),                 # Block 4
            (Joint(6,Right),Joint(8,Up),Joint(4,Left),Joint(2,Down)),   # Block 5
            (Joint(5,Left),Joint(3,Down)),                              # Block 6
            (Joint(8,Right),Joint(4,Down)),                             # Block 7
            (Joint(7,Left),Joint(5,Down))                               # Block 8
            )

# Build the domain
Dom = GridMultiBlock((D1,D2,D3,D4,D5,D6,D7,D8),joints)


# Block 1 BCs
Dl1 = FaADE.SATs.SAT_Dirichlet(Bx0,D1.Δx,Left,  order)
Dd1 = FaADE.SATs.SAT_Dirichlet(By0,D1.Δy,Down,  order)
# Block 2 BCs
Dd2 = FaADE.SATs.SAT_Dirichlet(By0,D2.Δx,Down,  order)
# Block 3 BCs
Dd3 = FaADE.SATs.SAT_Dirichlet(By0,D2.Δx,Down,  order)
Dr3 = FaADE.SATs.SAT_Dirichlet(Bxl,D2.Δy,Right, order)
# Block 4 BCs
Dl4 = FaADE.SATs.SAT_Dirichlet(Bx0,D4.Δx,Left,  order)
# Block 6 BCs
Dr6 = FaADE.SATs.SAT_Dirichlet(Bxl,D6.Δx,Right, order)
Du6 = FaADE.SATs.SAT_Dirichlet(Byh,D6.Δy,Up,    order) # Bottom of wedge
# Block 7 BCs
Dl7 = FaADE.SATs.SAT_Dirichlet(Bx0,D7.Δx,Left,  order)
Du7 = FaADE.SATs.SAT_Dirichlet(Byl,D7.Δx,Up,    order)
# Block 8 BCs
Du8 = FaADE.SATs.SAT_Dirichlet(Byl,D8.Δy,Up,    order)
Dr8 = FaADE.SATs.SAT_Dirichlet(Bxh,D8.Δx,Right, order) # Left of wedge


BD= Dict(1 => (Dl1,Dd1), 2=> (Dd2,), 3 => (Dd3,Dr3), 4 => (Dl4,), 6 => (Dr6,Du6), 7 => (Dl7,Du7), 8 => (Du8,Dr8))



P = Problem2D(order,u₀,K,K,Dom,BD,F,nothing)

println("---Solving---")
soln = solve(P,Dom,Δt,t)







println("Plotting")
using GLMakie

# e = zeros(size(Dom))
# for I in eachindex(Dom)
#     e[I] = exact(Dom[I]...,t)
# end
# u0 = [u₀(Dom[i]...) for i in eachindex(Dom)]

colourrange = (minimum(minimum.(soln.u[2])),maximum(maximum.(soln.u[2])))


f = Figure()
# ax1 = Axis3(f[1,1])
# surface!(ax1,Dom.gridx,Dom.gridy,e,colorbar=false)

ax = Axis3(f[1,1])
for I in eachindex(Dom.Grids)
    surface!(ax,Dom.Grids[I].gridx,Dom.Grids[I].gridy,soln.u[2][I],colorbar=false,colorrange=colourrange)
end

f

# p3 = surface(Dom.gridx,Dom.gridy,soln.u[2] .- e)




g = Figure()
gax = Axis(g[1,1])
lines!(gax,Dom.Grids[1].gridx[:,end],soln.u[2][1][end,:],label="Block 1")
lines!(gax,Dom.Grids[1].gridx[:,1],soln.u[2][2][1,:],label="Block 2")