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
# t = Δt

ωt = 1.0
ωx = 4.0
ωy = 3.0
cx = 1.0
cy = 0.5

Kx = 1.0
Ky = 1.0

θ = 0.5




# Plane wave solution
exact(x,y,t) = cos(2π*ωt*t) * sin(2π*x*ωx + cx)
u₀(x,y) = exact(x,y,0.0)
F(X,t) = begin
    x,y = X
    -2π*ωt*sin(2π*ωt*t)*sin(2π*x*ωx + cx) + 
    K * 4π^2 * ωx^2 * cos(2π*ωt*t)*sin(2π*x*ωx + cx)
end

# DIRICHLET
BxL(y,t) = cos(2π*ωt*t) * sin(cx)           #Boundary condition x=0
BxR(y,t) = cos(2π*ωt*t) * sin(2π*ωx + cx)   #Boundary condition x=Lx
ByL(x,t) = cos(2π*ωt*t) * sin(2π*ωx*x + cx) #Boundary condition y=0
ByR(x,t) = cos(2π*ωt*t) * sin(2π*ωx*x + cx) #Boundary condition y=Ly


#=
# Linear solution
exact(x,y,t) = x
u₀(x,y) = exact(x,y,0.0)
F = nothing

BxL(y,t) = 0.0
BxR(y,t) = 1.0
ByL(x,t) = x
ByR(x,t) = x
=#


#====== New solver 1 volume ======#
Dom = Grid2D([0.0,1.0],[0.0,1.0],nx,ny)

# New solver 1 volume
Dl = FaADE.SATs.SAT_Dirichlet(BxL,Dom.Δx,Left,order)
Dr = FaADE.SATs.SAT_Dirichlet(BxR,Dom.Δx,Right,order)
Dd = FaADE.SATs.SAT_Dirichlet(ByL,Dom.Δy,Down,order)
Du = FaADE.SATs.SAT_Dirichlet(ByR,Dom.Δy,Up,order)
BD = (Dl,Dr,Du,Dd)

# Pl = FaADE.SATs.SAT_Periodic(Dom.Δx,1,order,Left)
# Pr = FaADE.SATs.SAT_Periodic(Dom.Δx,1,order,Right)
# Pu = FaADE.SATs.SAT_Periodic(Dom.Δy,2,order,Up)
# Pd = FaADE.SATs.SAT_Periodic(Dom.Δy,2,order,Down)
# BD = (Pl,Pr,Pu,Pd)

P = Problem2D(order,u₀,K,K,Dom,BD,F,nothing)
println("---Solving 1 volume---")
soln = solve(P,Dom,Δt,t)


#====== New solver 2 volume ======#
println("2 volume")
D1 = Grid2D([0.0,0.5],[0.0,1.0],21,41)
D2 = Grid2D([0.5,1.0],[0.0,1.0],21,41)

joints = ((Joint(2,Right),),
            (Joint(1,Left),))
# 
Dom2V = GridMultiBlock((D1,D2),joints)

# Block 1 boundary
Dl1 = FaADE.SATs.SAT_Dirichlet(BxL,D1.Δx,Left,  order)
Dd1 = FaADE.SATs.SAT_Dirichlet(ByL,D2.Δy,Down,  order)
Du1 = FaADE.SATs.SAT_Dirichlet(ByR,D2.Δy,Up,    order)
# Block 2 boundary
Dr2 = FaADE.SATs.SAT_Dirichlet(BxR,D2.Δx,Right, order)
Dd2 = FaADE.SATs.SAT_Dirichlet(ByL,D2.Δy,Down,  order)
Du2 = FaADE.SATs.SAT_Dirichlet(ByR,D2.Δy,Up,    order)

BD2 = Dict(1 => (Dl1,Du1,Dd1), 2 => (Dr2,Du2,Dd2))

P2V = Problem2D(order,u₀,K,K,Dom2V,BD2,F,nothing)

println("---Solving 2 volume---")
soln2V = solve(P2V,Dom2V,Δt,t)




#====== New solver 3 volume ======#
println("3 volume")

D1 = Grid2D([0.0,0.4],[0.0,1.0],21,41)
D2 = Grid2D([0.4,0.6],[0.0,1.0],21,41)
D3 = Grid2D([0.6,1.0],[0.0,1.0],21,41)

joints = ((Joint(2,Right),),
            (Joint(1,Left),Joint(3,Right)),
            (Joint(2,Left),))

Dom3V = GridMultiBlock((D1,D2,D3),joints)

# Block 1 boundary
Dl1 = FaADE.SATs.SAT_Dirichlet(BxL,D1.Δx,Left,  order)
Dd1 = FaADE.SATs.SAT_Dirichlet(ByL,D1.Δy,Down,  order)
Du1 = FaADE.SATs.SAT_Dirichlet(ByR,D1.Δy,Up,    order)
# Block 2 boundary
Dd2 = FaADE.SATs.SAT_Dirichlet(ByL,D2.Δy,Down,  order)
Du2 = FaADE.SATs.SAT_Dirichlet(ByR,D2.Δy,Up,    order)
# Block 3 boundary
Dr3 = FaADE.SATs.SAT_Dirichlet(BxR,D2.Δx,Right, order)
Dd3 = FaADE.SATs.SAT_Dirichlet(ByL,D3.Δy,Down,  order)
Du3 = FaADE.SATs.SAT_Dirichlet(ByR,D3.Δy,Up,    order)


BD3 = Dict(1 => (Dl1,Du1,Dd1), 2=> (Du2, Dd2), 3 => (Dr3,Du3,Dd3))

P3V = Problem2D(order,u₀,K,K,Dom3V,BD3,F,nothing)

println("---Solving 3 volume---")
soln3V = solve(P3V,Dom3V,Δt,t)



#====== New solver 4 volume ======#
println("4 volume")
D1 = Grid2D([0.0,0.5],[0.0,0.5],21,21)
D2 = Grid2D([0.5,1.0],[0.0,0.5],21,21)
D3 = Grid2D([0.0,0.5],[0.5,1.0],21,21)
D4 = Grid2D([0.5,1.0],[0.5,1.0],21,21)

joints = ((Joint(2,Right),Joint(3,Up)),
            (Joint(1,Left),Joint(4,Up)),
            (Joint(4,Right),Joint(1,Down)),
            (Joint(3,Left),Joint(2,Down)))

Dom4V = GridMultiBlock((D1,D2,D3,D4),joints)

# Block 1 BCs
Dl1 = FaADE.SATs.SAT_Dirichlet(BxL,D1.Δx,Left,  order)
Dd1 = FaADE.SATs.SAT_Dirichlet(ByL,D1.Δy,Down,  order)
# Block 2 BCs
Dr2 = FaADE.SATs.SAT_Dirichlet(BxR,D2.Δx,Right, order)
Dd2 = FaADE.SATs.SAT_Dirichlet(ByL,D2.Δy,Down,  order)
# Block 3 BCs
Dl3 = FaADE.SATs.SAT_Dirichlet(BxL,D3.Δx,Left,  order)
Du3 = FaADE.SATs.SAT_Dirichlet(ByR,D3.Δy,Up,    order)
# Block 4 BCs
Dr4 = FaADE.SATs.SAT_Dirichlet(BxL,D4.Δx,Right, order)
Du4 = FaADE.SATs.SAT_Dirichlet(ByR,D4.Δy,Up,    order)

BD4= Dict(1 => (Dl1,Dd1), 2=> (Dd2,Dr2), 3 => (Dl3,Du3), 4 => (Dr4,Du4))



P4V = Problem2D(order,u₀,K,K,Dom4V,BD4,F,nothing)

println("---Solving 4 volume---")
soln4V = solve(P4V,Dom4V,Δt,t)







println("Plotting")
using GLMakie

e = zeros(size(Dom))
for I in eachindex(Dom)
    e[I] = exact(Dom[I]...,t)
end
u0 = [u₀(Dom[i]...) for i in eachindex(Dom)]



colourrange = (minimum(soln.u[2]),maximum(soln.u[2]))

f = Figure()
ax1 = Axis3(f[1,1])
# surface!(ax1,Dom.gridx,Dom.gridy,e,colorbar=false)
surface!(ax1,Dom.gridx,Dom.gridy,soln.u[2],colorbar=false,colorrange=colourrange)
# surface!(ax1,Dom.gridx,Dom.gridy,soln.u[2].-e,colorbar=false)

ax2 = Axis3(f[1,2])
surface!(ax2,Dom2V.Grids[1].gridx,Dom2V.Grids[1].gridy,soln2V.u[2][1],colorbar=false,colorrange=colourrange)
surface!(ax2,Dom2V.Grids[2].gridx,Dom2V.Grids[2].gridy,soln2V.u[2][2],colorbar=false,colorrange=colourrange)

ax3 = Axis3(f[2,1])
surface!(ax3,Dom3V.Grids[1].gridx,Dom3V.Grids[1].gridy,soln3V.u[2][1],colorbar=false,colorrange=colourrange)
surface!(ax3,Dom3V.Grids[2].gridx,Dom3V.Grids[2].gridy,soln3V.u[2][2],colorbar=false,colorrange=colourrange)
surface!(ax3,Dom3V.Grids[3].gridx,Dom3V.Grids[3].gridy,soln3V.u[2][3],colorbar=false,colorrange=colourrange)

ax4 = Axis3(f[2,2])
surface!(ax4,Dom4V.Grids[1].gridx,Dom4V.Grids[1].gridy,soln4V.u[2][1],colorbar=false,colorrange=colourrange)
surface!(ax4,Dom4V.Grids[2].gridx,Dom4V.Grids[2].gridy,soln4V.u[2][2],colorbar=false,colorrange=colourrange)
surface!(ax4,Dom4V.Grids[3].gridx,Dom4V.Grids[3].gridy,soln4V.u[2][3],colorbar=false,colorrange=colourrange)
surface!(ax4,Dom4V.Grids[4].gridx,Dom4V.Grids[4].gridy,soln4V.u[2][4],colorbar=false,colorrange=colourrange)

f

# p3 = surface(Dom.gridx,Dom.gridy,soln.u[2] .- e)





g = Figure()
gax1 = Axis3(g[1,1])

wireframe!(gax1,Dom3V.Grids[1].gridx,Dom3V.Grids[1].gridy,soln3V.u[2][1])
wireframe!(gax1,Dom3V.Grids[2].gridx,Dom3V.Grids[2].gridy,soln3V.u[2][2])
wireframe!(gax1,Dom3V.Grids[3].gridx,Dom3V.Grids[3].gridy,soln3V.u[2][3])




