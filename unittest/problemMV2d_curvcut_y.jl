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
t = 0.26

ωt = 3.0
ωx = 3.5 
ωy = 2.5
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
Bxy(X,t) = cos(2π*ωt*t) * sin(2π*X[1]*ωx + cx) * sin(2π*ωy*X[2] + cy)   #Boundary condition x=0
# NEUMANN
# BxLũ(y,t) = 2π*ωx * K * cos(2π*t) * cos(cx)             * sin(2π*y*ωy + cy) #Boundary condition x=0

# exact(R,Z,t) = exp.( -(R.^2 + Z.^2) / 0.1 )
# exact(R,Z,t) = sin(2π*R*ωx + cx)
# u₀(R,Z) = exact(R,Z,0.0)
# F(X,t) = 4π^2*ωx^2*sin(2π*X[1]*ωx + cx)
# F(X,t) = 0.0
# Bxy(X,t) = sin(2π*X[1]*ωx + cx)


#====== New solver 4 volume ======#
println("Curvilinear volume")

# Right domain
D1 = Grid2D(u->u*[cos(7π/4), sin(7π/4)] - [0.0,0.25],
            v->[0.0, v/2 - 0.25],
            v->[cos(v*(9π/4 - 7π/4) + 7π/4), sin(v*(9π/4 - 7π/4) + 7π/4)] - 
                (1-v)*[0.0,0.25] + v*[0.0,0.25],
            u->u*[cos(π/4), sin(π/4)] + [0.0, 0.25], # top of D1
            nx,ny)


# Left domain
D2 = Grid2D(u->u*[cos(π/4), sin(π/4)] + [0.0, 0.25], # bottom - same as D1 top -- done
            v->[cos(v*(3π/4 - 5π/4) + 5π/4), sin(v*(3π/4 - 5π/4) + 5π/4)] - 
                [cos(5π/4), sin(5π/4) - 0.25], # left - shifted up 0.25 -- done
            v->v*[0.0, 0.5] + [cos(π/4), sin(π/4) + 0.25], # shifted to top of D1 -- done
            u->u*([cos(π/4), sin(π/4) + 3/4] - [cos(3π/4) - cos(5π/4), sin(3π/4) - sin(5π/4) + 0.25]) +
                [cos(3π/4) - cos(5π/4), sin(3π/4) - sin(5π/4) + 0.25], # top - shifted up 0.25
            nx,ny)


joints = ((Joint(2,Up),),
            (Joint(1,Down),),)


Dom = GridMultiBlock((D1,D2),joints)

Dl1 = SAT_Dirichlet(Bxy, D1.Δx, Left, order, D1.Δy, :Curvilinear) # Block 3 BCs
Dr1 = SAT_Dirichlet(Bxy, D1.Δx, Right,order, D1.Δy, :Curvilinear) # 
Dd1 = SAT_Dirichlet(Bxy, D1.Δy, Down, order, D1.Δx, :Curvilinear) # Block 5 BCs

Dr2 = SAT_Dirichlet(Bxy, D2.Δx, Right,order, D2.Δy, :Curvilinear) # Block 2 BCs
Dl2 = SAT_Dirichlet(Bxy, D2.Δx, Left, order, D2.Δy, :Curvilinear) # Block 5 BCs
Du2 = SAT_Dirichlet(Bxy, D2.Δy, Up,   order, D2.Δx, :Curvilinear) # Block 3 BCs


BD = Dict(1 => (Dl1,Dr1,Dd1), 2 => (Dr2,Dl2,Du2))



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
    e[1][I] = exact(D1.gridx[I],D1.gridy[I],t)
end
for I in eachindex(D2)
    e[2][I] = exact(D2.gridx[I],D2.gridy[I],t)
end

@show norm(e[1] .- soln.u[2][1])/norm(e[1])
@show norm(e[2] .- soln.u[2][2])/norm(e[2])


f = Figure()

ax1 = Axis3(f[1,1])
# surface!(ax1,Dom.Grids[1].gridx, Dom.Grids[1].gridy, soln.u[1][1],colorbar=false, colorrange=colourrange)
# surface!(ax1,Dom.Grids[2].gridx, Dom.Grids[2].gridy, soln.u[1][2],colorbar=false, colorrange=colourrange)
surface!(ax1,Dom.Grids[1].gridx, Dom.Grids[1].gridy, e[1] .-  soln.u[2][1],colorbar=false, colorrange=colourrange)
surface!(ax1,Dom.Grids[2].gridx, Dom.Grids[2].gridy, e[2] .-  soln.u[2][2],colorbar=false, colorrange=colourrange)

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
lines!(axh1,Dom.Grids[1].gridx[:,21],soln.u[1][1][:,21])
lines!(axh1,Dom.Grids[2].gridx[:,21],soln.u[1][2][:,21])
axh2 = Axis(h[1,2])
lines!(axh2,Dom.Grids[1].gridx[:,21],soln.u[2][1][:,21])
lines!(axh2,Dom.Grids[2].gridx[:,21],soln.u[2][2][:,21])

