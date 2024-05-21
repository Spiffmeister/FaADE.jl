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
ωx = 5.0
ωy = 3.0
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
Bx0(y,t) = cos(2π*ωt*t) * sin(cx) * sin(2π*ωy*y + cy)           #Boundary condition x=0
Bxl(y,t) = cos(2π*ωt*t) * sin(2π*ωx + cx) * sin(2π*ωy*y + cy)   #Boundary condition x=Lx
By0(x,t) = cos(2π*ωt*t) * sin(2π*ωx*x + cx) * sin(cy)           #Boundary condition y=0
Byl(x,t) = cos(2π*ωt*t) * sin(2π*ωx*x + cx) * sin(2π*ωy + cy)   #Boundary condition y=Ly

Bxh(y,t) = cos(2π*ωt*t) * sin(π*ωx + cx) * sin(2π*y*ωy + cy) #Boundary condition x=0.5
Byh(x,t) = cos(2π*ωt*t) * sin(2π*x*ωx + cx) * sin(π*ωy + cy) #Boundary condition y=0.5

#====== New solver 4 volume ======#
#   -----
#   | 3 |
#   ---------
#   | 1 | 2 |
#   ---------   
#
D1 = Grid2D([0.0,0.5],[0.0,0.5],21,21)
D2 = Grid2D([0.5,1.0],[0.0,0.5],21,21)
D3 = Grid2D([0.0,0.5],[0.5,1.0],21,21)

joints = ((Joint(2,Right),Joint(3,Up)),
            (Joint(1,Left),),
            (Joint(1,Down),))

Dom = GridMultiBlock((D1,D2,D3),joints)


# Block 1 BCs
Dl1 = FaADE.SATs.SAT_Dirichlet(Bx0,D1.Δx,Left,    order)
Dd1 = FaADE.SATs.SAT_Dirichlet(By0,D1.Δy,Down,    order)
# Block 2 BCs
Dr2 = FaADE.SATs.SAT_Dirichlet(Bxl,D2.Δx,Right,   order)
Dd2 = FaADE.SATs.SAT_Dirichlet(By0,D2.Δy,Down,    order)
Du2 = FaADE.SATs.SAT_Dirichlet(Byh,D2.Δy,Up,      order)
# Block 3 BCs
Dl3 = FaADE.SATs.SAT_Dirichlet(Bx0,D3.Δx,Left,    order)
Du3 = FaADE.SATs.SAT_Dirichlet(Byl,D3.Δy,Up,      order)
Dr3 = FaADE.SATs.SAT_Dirichlet(Bxh,D3.Δx,Right,   order)

BD= Dict(1 => (Dl1,Dd1), 2=> (Dr2,Dd2,Du2), 3 => (Dl3,Du3,Dr3))



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
surface!(ax,Dom.Grids[1].gridx,Dom.Grids[1].gridy,soln.u[2][1],colorbar=false,colorrange=colourrange)
surface!(ax,Dom.Grids[2].gridx,Dom.Grids[2].gridy,soln.u[2][2],colorbar=false,colorrange=colourrange)
surface!(ax,Dom.Grids[3].gridx,Dom.Grids[3].gridy,soln.u[2][3],colorbar=false,colorrange=colourrange)

f

# p3 = surface(Dom.gridx,Dom.gridy,soln.u[2] .- e)




g = Figure()
gax = Axis(g[1,1])
lines!(gax,Dom.Grids[1].gridx[:,end],soln.u[2][1][end,:],label="Block 1")
lines!(gax,Dom.Grids[1].gridx[:,1],soln.u[2][2][1,:],label="Block 2")