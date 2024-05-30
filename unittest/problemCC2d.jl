
using Revise
using FaADE
# using BenchmarkTools
# using ProfileView
# using Cthulhu


order = 2
K = 1.0

nx = ny = 41

Δt = 1e-3
# t = Δt
t = 0.76

# u₀(x,y) = x.^2
u₀(x,y) = exp.(-((x-0.5)^2 + (y-0.5)^2) / 0.02)


# New solver 1 volume
cbottom(u) = [u,0.0]
cleft(v) = [0.0,v]
cright(v) = [1.0,v]
ctop(u) = [u,1.0]

Dom1V = Grid2D(cbottom,cleft,cright,ctop,ceil(Int,nx/2),ny)

Dl = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,Dom1V.Δx,Left, order)
Dr = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,Dom1V.Δx,Right,order)
Du = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,Dom1V.Δy,Up,   order)
Dd = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,Dom1V.Δy,Down, order)
BD1V = (Dl,Dr,Du,Dd)

P1V = Problem2D(order,u₀,K,K,Dom1V,BD1V)

println("---Solving 1 volume---")
soln1V = solve(P1V,Dom1V,Δt,t)
# @benchmark solve($P1V,$Dom1V,$Δt,$t)



D1 = Grid2D([0.0,0.5],[0.0,1.0],ceil(Int,nx/2),ny, coord=CurvilinearMetric)
D2 = Grid2D([0.5,1.0],[0.0,1.0],ceil(Int,nx/2),ny, coord=CurvilinearMetric)

joints = (
    (Joint(2,Right),),
    (Joint(1,Left),)
)



Dom2V = GridMultiBlock((D1,D2),joints)

Dl1 = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,D1.Δx,Left,    order, D1.Δy, :Curvilinear)
Du1 = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,D1.Δy,Up,      order, D1.Δx, :Curvilinear)
Dd1 = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,D1.Δy,Down,    order, D1.Δx, :Curvilinear)

Dr2 = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,D2.Δx,Right,   order, D2.Δy, :Curvilinear)
Du2 = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,D2.Δy,Up,      order, D2.Δx, :Curvilinear)
Dd2 = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,D2.Δy,Down,    order, D2.Δx, :Curvilinear)

BD2V = Dict(1=>(Dl1,Du1,Dd1), 2=>(Dr2,Du2,Dd2))

P2V = Problem2D(order,u₀,K,K,Dom2V,BD2V)

println("---Solving 2 volume---")
soln2V = solve(P2V,Dom2V,Δt,t)
# @benchmark solve($P2V,$Dom2V,$Δt,$t)





#====== 2 volume double arc ======#
u₀(x,y) = exp.(-(x^2 + y^2) / 0.02)

T = FaADE.Grid.Torus([2.0],[2.0],[1],[0])

# Left domain
D1 = Grid2D(u->u*(T(5π/4,0.0) + [0.0, 1.0]) + [0.0, 1.0],
            v->T(v*(3π/4 - 5π/4) + 5π/4,0.0),
            v->[0.0,2v - 1.0],
            u->u*([0.0, 1.0] + T(3π/4,0.0)) + T(3π/4,0.0),
            nx,nx)

# Right domain
D2 = Grid2D(u->u*(T(7π/4,0.0) + [0.0, 1.0]) + [0.0, -1.0],
            v->[0.0, 2v - 1.0],
            v->T(v*(9π/4 - 7π/4) + 7π/4,0.0),
            u->u*(T(π/4,0.0) - [0.0, 1.0]) + [0.0, 1.0],
            nx,nx)

joints = (
    (Joint(2,Right),),
    (Joint(1,Left),)
)

Dom2VC = GridMultiBlock((D1,D2),joints)

Dl1C = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,D1.Δx,Left,    order, D1.Δy, :Curvilinear)
Du1C = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,D1.Δy,Up,      order, D1.Δx, :Curvilinear)
Dd1C = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,D1.Δy,Down,    order, D1.Δx, :Curvilinear)

Dr2C = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,D2.Δx,Right,   order, D2.Δy, :Curvilinear)
Du2C = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,D2.Δy,Up,      order, D2.Δx, :Curvilinear)
Dd2C = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,D2.Δy,Down,    order, D2.Δx, :Curvilinear)

BD2VC = Dict(1=>(Dl1C,Du1C,Dd1C), 2=>(Dr2C,Du2C,Dd2C))

P2VC = Problem2D(order,u₀,K,K,Dom2VC,BD2VC)

println("---Solving 2 volume curvilinear---")
soln3 = solve(P2VC,Dom2VC,Δt,t)


# g = Figure(); ax = Axis(g[1,1]);
# scatter!(ax,Dom2VC.Grids[1].gridx[:],Dom2VC.Grids[1].gridy[:],marker=:circle)
# scatter!(ax,Dom2VC.Grids[2].gridx[:],Dom2VC.Grids[2].gridy[:],marker=:circle)




#====== 2 volume square and arc ======#
u₀(x,y) = exp.(-(x^2 + y^2) / 0.02)

T = FaADE.Grid.Torus([2.0],[2.0],[1],[0])

# Left domain
D1 = Grid2D([-1.0,0.0],[-1.0,1.0],nx,ny,coord=CurvilinearMetric)

# Right domain
D2 = Grid2D(u->u*(T(7π/4,0.0) + [0.0, 1.0]) + [0.0, -1.0],
            v->[0.0, 2v - 1.0],
            v->T(v*(9π/4 - 7π/4) + 7π/4, 0.0),
            u->u*(T(π/4,0.0) - [0.0, 1.0]) + [0.0, 1.0],
            nx,nx)

joints = (
    (Joint(2,Right),),
    (Joint(1,Left),)
)

Dom4 = GridMultiBlock((D1,D2),joints)

Dl1 = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,D1.Δx,Left,    order, D1.Δy, :Curvilinear)
Du1 = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,D1.Δy,Up,      order, D1.Δx, :Curvilinear)
Dd1 = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,D1.Δy,Down,    order, D1.Δx, :Curvilinear)

Dr2 = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,D2.Δx,Right,   order, D2.Δy, :Curvilinear)
Du2 = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,D2.Δy,Up,      order, D2.Δx, :Curvilinear)
Dd2 = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,D2.Δy,Down,    order, D2.Δx, :Curvilinear)

BD4 = Dict(1=>(Dl1,Du1,Dd1), 2=>(Dr2,Du2,Dd2))

P4 = Problem2D(order,u₀,K,K,Dom4,BD4)

println("---Solving 2 volume curvilinear---")
soln4 = solve(P4,Dom4,Δt,t)








println("---Plotting---")
using GLMakie


g = Figure(); 
gax1 = Axis(g[1,1]);
scatter!(gax1,Dom2V.Grids[1].gridx[:],Dom2V.Grids[1].gridy[:],marker=:circle)
scatter!(gax1,Dom2V.Grids[2].gridx[:],Dom2V.Grids[2].gridy[:],marker=:circle)
gax2 = Axis(g[1,2])
scatter!(gax2,Dom2VC.Grids[1].gridx[:],Dom2VC.Grids[1].gridy[:],marker=:circle)
scatter!(gax2,Dom2VC.Grids[2].gridx[:],Dom2VC.Grids[2].gridy[:],marker=:circle)
gax3 = Axis(g[1,3])
scatter!(gax3,Dom4.Grids[1].gridx[:],Dom4.Grids[1].gridy[:],marker=:circle)
scatter!(gax3,Dom4.Grids[2].gridx[:],Dom4.Grids[2].gridy[:],marker=:circle)





f = Figure()
Ax1=Axis3(f[1,1])

surface!(Ax1,Dom1V.gridx,Dom1V.gridy,soln1V.u[2])

Ax2=Axis3(f[1,2])
surface!(Ax2,Dom2V.Grids[1].gridx,Dom2V.Grids[1].gridy,soln2V.u[2][1])
surface!(Ax2,Dom2V.Grids[2].gridx,Dom2V.Grids[2].gridy,soln2V.u[2][2])
lines!(Ax2,Dom2V.Grids[1].gridx[end,:],Dom2V.Grids[1].gridy[end,:],soln2V.u[2][1][end,:],color=:red)

Ax3=Axis3(f[2,1])
surface!(Ax3,Dom2VC.Grids[1].gridx,Dom2VC.Grids[1].gridy,soln3.u[2][1])
surface!(Ax3,Dom2VC.Grids[2].gridx,Dom2VC.Grids[2].gridy,soln3.u[2][2])
lines!(Ax3,Dom2VC.Grids[1].gridx[end,:],Dom2VC.Grids[1].gridy[end,:],soln3.u[2][1][end,:],color=:red)

Ax4=Axis3(f[2,2])
surface!(Ax4,Dom4.Grids[1].gridx,Dom4.Grids[1].gridy,soln4.u[2][1])
surface!(Ax4,Dom4.Grids[2].gridx,Dom4.Grids[2].gridy,soln4.u[2][2])
lines!(Ax4,Dom4.Grids[1].gridx[end,:],Dom4.Grids[1].gridy[end,:],soln4.u[2][1][end,:],color=:red)

f




