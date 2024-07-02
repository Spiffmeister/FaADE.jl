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
exact(x,y,t) = cos(2π*ωt*t) * sin(2π*x*ωx + cx) * sin(2π*y*ωy + cy)
u₀(x,y) = exact(x,y,0.0)
F(X,t) = begin
    x,y = X
    -2π*ωt*sin(2π*ωt*t)*sin(2π*x*ωx + cx)*sin(2π*y*ωy + cy) + 
    K * 4π^2 * ωx^2 * cos(2π*ωt*t)*sin(2π*x*ωx + cx)*sin(2π*y*ωy + cy) + 
    K * 4π^2 * ωy^2 * cos(2π*ωt*t)*sin(2π*x*ωx + cx)*sin(2π*y*ωy + cy)
end
# DIRICHLET
Bxy(X,t) = cos(2π*ωt*t) * sin(2π*X[1]*ωx + cx) * sin(2π*ωy*X[2] + cy)   #Boundary condition x=0
# NEUMANN
# BxLũ(y,t) = 2π*ωx * K * cos(2π*t) * cos(cx)             * sin(2π*y*ωy + cy) #Boundary condition x=0

# exact(R,Z,t) = exp.( -(R.^2 + Z.^2) / 0.1 )

# exact(R,Z,t) = sin(2π*R*ωx + cx)
# u₀(R,Z) = exact(R,Z,0.0)
# F(X,t) = 0.0
# Bxy(X,t) = exp.( -(X[1].^2 + X[2].^2) / 0.1 )
# Bxy(X,t) = sin(2π*X[1]*ωx + cx)


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
D2 = Grid2D(u->[1.0u + 0.25,-0.25],
            v->[0.25,       v*0.5 - 0.25],
            v->[1.0,       v*0.5 - 0.25],
            u->[1.0u + 0.25,0.25],
            nx,ny)



joints = ((Joint(2,Right),),
            (Joint(1,Left),),)


Dom = GridMultiBlock((D1,D2),joints)

Byd(X,t) = cos(2π*ωt*t) * sin(2π*ωx*X[1] * ωx + cx) * sin(-π/2 * ωy + cy)
Byu(X,t) = cos(2π*ωt*t) * sin(2π*ωx*X[1] + cx) * sin(π/2 * ωy + cy)

Bxy(X,t) = cos(2π*ωt*t) * sin(2π*X[1]*ωx + cx) * sin(2π*ωy*X[2] + cy)   #Boundary condition x=0

Dd1 = SAT_Dirichlet(Byd, D1.Δy, Down, order, D1.Δx, :Curvilinear) # Block 5 BCs
Du1 = SAT_Dirichlet(Byu, D1.Δy, Up,   order, D1.Δx, :Curvilinear) # Block 3 BCs
Dl1 = SAT_Dirichlet(Bxy, D1.Δx, Left, order, D1.Δy, :Curvilinear)

Dr2 = SAT_Dirichlet(Bxy, D2.Δy, Right,order, D2.Δx, :Curvilinear) # Block 2 BCs
Dd2 = SAT_Dirichlet(Byd, D2.Δx, Down, order, D2.Δy, :Curvilinear) # Block 5 BCs
Du2 = SAT_Dirichlet(Byu, D2.Δx, Up,   order, D2.Δy, :Curvilinear) # Block 3 BCs

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




e = [zeros(size(D1)),zeros(size(D2))]
for i in eachindex(D1)
    e[1][i] = exact(D1.gridx[i],D1.gridy[i],t)
end
for i in eachindex(D2)
    e[2][i] = exact(D2.gridx[i],D2.gridy[i],t)
end


g = Figure()
axg1 = Axis3(g[1,1])
surface!(axg1,Dom.Grids[1].gridx, Dom.Grids[1].gridy, soln.u[2][1] .- e[1],colorbar=false)
surface!(axg1,Dom.Grids[2].gridx, Dom.Grids[2].gridy, soln.u[2][2] .- e[2],colorbar=false)







f = Figure()

ax1 = Axis3(f[1,1])
# surface!(ax1,Dom.Grids[1].gridx, Dom.Grids[1].gridy, soln.u[1][1],colorbar=false, colorrange=colourrange)
# surface!(ax1,Dom.Grids[2].gridx, Dom.Grids[2].gridy, soln.u[1][2],colorbar=false, colorrange=colourrange)
surface!(ax1,Dom.Grids[1].gridx, Dom.Grids[1].gridy, e[1],colorbar=false, colorrange=colourrange)
surface!(ax1,Dom.Grids[2].gridx, Dom.Grids[2].gridy, e[2],colorbar=false, colorrange=colourrange)

# scatter!(ax1,D1.gridx[:],D1.gridy[:],-ones(length(D1)),markersize=1.5)
# scatter!(ax1,D2.gridx[:],D2.gridy[:],-ones(length(D2)),markersize=1.5)

ax2 = Axis3(f[1,2])
surface!(ax2,Dom.Grids[1].gridx, Dom.Grids[1].gridy, soln.u[2][1],colorbar=false, colorrange=colourrange)
surface!(ax2,Dom.Grids[2].gridx, Dom.Grids[2].gridy, soln.u[2][2],colorbar=false, colorrange=colourrange)

# scatter!(ax2,D1.gridx[:],D1.gridy[:],-ones(length(D1)),markersize=1.5)
# scatter!(ax2,D2.gridx[:],D2.gridy[:],-ones(length(D2)),markersize=1.5)

f




h = Figure()
axh1 = Axis(h[1,1])
lines!(axh1,Dom.Grids[1].gridy[end,:],  soln.u[1][1][end,:])
scatter!(axh1,Dom.Grids[2].gridy[1,:],    soln.u[1][2][1,:])
axh2 = Axis(h[1,2])
lines!(axh2,Dom.Grids[1].gridy[end,:],soln.u[2][1][end,:])
scatter!(axh2,Dom.Grids[2].gridy[1,:],soln.u[2][2][1,:])





j = Figure()
axj1 = Axis(j[1,1])
lines!(axj1,Dom.Grids[1].gridx[:,21],soln.u[1][1][:,21])
lines!(axj1,Dom.Grids[2].gridx[:,21],soln.u[1][2][:,21])
axj2 = Axis(j[1,2])
lines!(axj2,Dom.Grids[1].gridx[:,21],soln.u[2][1][:,21])
lines!(axj2,Dom.Grids[2].gridx[:,21],soln.u[2][2][:,21])
