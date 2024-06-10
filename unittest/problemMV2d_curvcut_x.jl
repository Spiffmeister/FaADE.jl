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
F(X,t) = 4π^2*ωx^2*sin(2π*X[1]*ωx + cx)
Bxy(X,t) = sin(2π*X[1]*ωx + cx)


#====== New solver 4 volume ======#
println("Curvilinear volume")

# Right domain
D1 = Grid2D(u->u*[0.5, 0.0] - [0.0,0.0],
            v->v*[cos(3π/4),sin(3π/4)] + [0.0,0.0],
            v->v*[cos(π/4),sin(π/4)] + [0.5,0.0],
            u->[cos(u*(π/4 - 3π/4) + 3π/4), sin(u*(π/4 - 3π/4) + 3π/4)] + 
                u*[0.5,0.0], # top of D1
            nx,ny)


# Left domain
D2 = Grid2D(u->[cos(u*(-5π/4 + 7π/4) + 5π/4), sin(u*(-5π/4 + 7π/4) + 5π/4)] - 
                [cos(7π/4), sin(7π/4)] + u*[0.5,0.0] - [0.5,0.0], # bottom - same as D1 top -- done
            v->v*[cos(3π/4),sin(3π/4)] + [0.0,0.0], # left - shifted up 0.25 -- done
            v->v*[cos(π/4),sin(π/4)] + [cos(5π/4)-cos(7π/4)-0.5,sin(5π/4)-sin(7π/4)], # shifted to top of D1 -- done
            u->u*[0.5, 0.0] + [cos(3π/4) - 0.5, sin(3π/4)], # top - shifted up 0.25
            nx,ny)


joints = ((Joint(2,Left),),
            (Joint(1,Right),),)

Dom = GridMultiBlock((D1,D2),joints)

Dr1 = SAT_Dirichlet(Bxy, D1.Δx, Right,  order, D1.Δy, :Curvilinear) #
Dd1 = SAT_Dirichlet(Bxy, D1.Δy, Down,   order, D1.Δx, :Curvilinear) #
Du1 = SAT_Dirichlet(Bxy, D1.Δy, Up,     order, D1.Δx, :Curvilinear) #

Dl2 = SAT_Dirichlet(Bxy, D2.Δx, Left, order, D2.Δy, :Curvilinear) #
Du2 = SAT_Dirichlet(Bxy, D2.Δy, Up,   order, D2.Δx, :Curvilinear) #
Dd2 = SAT_Dirichlet(Bxy, D2.Δy, Down, order, D2.Δx, :Curvilinear) #


BD = Dict(1 => (Dr1,Dd1,Du1), 2 => (Dl2,Du2,Dd2))



using GLMakie
gridfig = Figure()
gridfix_ax = Axis(gridfig[1,1])
scatter!(gridfix_ax,D1.gridx[:],D1.gridy[:],markersize=2.5)
scatter!(gridfix_ax,D2.gridx[:],D2.gridy[:],markersize=2.5)
gridfig


P = Problem2D(order,u₀,K,K,Dom,BD,F,nothing)

println("---Solving 4 volume---")
soln = solve(P,Dom,Δt,t)


colourrange = (minimum(minimum.(soln.u[2])),maximum(maximum.(soln.u[2])))

e = [zeros(size(D1)),zeros(size(D2))]
for I in eachindex(D1)
    # e[1][I] = exact(D1.gridx[I],D1.gridy[I],t)
    e[1][I] = exact(D1.gridx[I],D1.gridy[I],0.0)
end
for I in eachindex(D2)
    e[2][I] = exact(D2.gridx[I],D2.gridy[I],0.0)
end


f = Figure()

ax1 = Axis3(f[1,1])
# surface!(ax1,Dom.Grids[1].gridx, Dom.Grids[1].gridy, soln.u[1][1],colorbar=false, colorrange=colourrange)
# surface!(ax1,Dom.Grids[2].gridx, Dom.Grids[2].gridy, soln.u[1][2],colorbar=false, colorrange=colourrange)
surface!(ax1,Dom.Grids[1].gridx, Dom.Grids[1].gridy, e[1],colorbar=false)#, colorrange=colourrange)
surface!(ax1,Dom.Grids[2].gridx, Dom.Grids[2].gridy, e[2],colorbar=false)#, colorrange=colourrange)

scatter!(ax1,D1.gridx[:],D1.gridy[:],-ones(length(D1)),markersize=1.5)
scatter!(ax1,D2.gridx[:],D2.gridy[:],-ones(length(D2)),markersize=1.5)

ax2 = Axis3(f[1,2])
surface!(ax2,Dom.Grids[1].gridx, Dom.Grids[1].gridy, soln.u[2][1],colorbar=false, colorrange=colourrange)
surface!(ax2,Dom.Grids[2].gridx, Dom.Grids[2].gridy, soln.u[2][2],colorbar=false, colorrange=colourrange)

scatter!(ax2,D1.gridx[:],D1.gridy[:],-ones(length(D1)),markersize=1.5)
scatter!(ax2,D2.gridx[:],D2.gridy[:],-ones(length(D2)),markersize=1.5)

f




h = Figure()
axh1 = Axis(h[1,1])
lines!(axh1,Dom.Grids[1].gridx[:,21],soln.u[1][1][:,21])
lines!(axh1,Dom.Grids[2].gridx[:,21],soln.u[1][2][:,21])
axh2 = Axis(h[1,2])
lines!(axh2,Dom.Grids[1].gridx[:,21],soln.u[2][1][:,21])
lines!(axh2,Dom.Grids[2].gridx[:,21],soln.u[2][2][:,21])

