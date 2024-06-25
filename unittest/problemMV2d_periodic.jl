using Revise
using FaADE
# using BenchmarkTools
# using ProfileView
# using Cthulhu

# Simulation parameters
order = 2
K = 1.0

Δt = 0.01
# t = 100.0
t = 1.0

# Set initial condition
u₀(x,y) = x.^2
# u₀(x,y) = exp.(-((x-0.5)^2 + (y-0.5)^2) / 0.02)





# Original solver
Dom1V = Grid2D([0.0,1.0],[-π,π],21,21)

BoundaryLeft    = SAT_Dirichlet((y,t)->0.0, Dom1V.Δx, Left, order)
BoundaryRight   = SAT_Dirichlet((y,t)->1.0, Dom1V.Δx, Right, order)
BoundaryUp      = SAT_Periodic(Dom1V.Δy, order, Up)
BoundaryDown    = SAT_Periodic(Dom1V.Δy, order, Down)

BC = (BoundaryLeft, BoundaryRight, BoundaryUp, BoundaryDown)

P1V = Problem2D(order,u₀,K,K,Dom1V,BC)
println("---Solving old---")
soln1V = solve(P1V,Dom1V,Δt,t) #-Δt to ensure ends at the same time as new methods




# New solver 1 volume
D1 = Grid2D([0.0,0.5],[-π,π],11,21)
D2 = Grid2D([0.5,1.0],[-π,π],11,21)
Dom2V = GridMultiBlock((D1,D2),(
        (Joint(2,Right),),(Joint(1,Left),)
        ))

Dl1 = SAT_Dirichlet((x,t)->0.0, D1.Δx, Left, order)
Dr2 = SAT_Dirichlet((x,t)->1.0, D2.Δx, Right, order)

Pu1 = SAT_Periodic(D1.Δy, order, Up)
Pd1 = SAT_Periodic(D1.Δy, order, Down)

Pu2 = SAT_Periodic(D2.Δy, order, Up)
Pd2 = SAT_Periodic(D2.Δy, order, Down)

BC2V = Dict(1 => (Dl1, Pu1, Pd1), 2 => (Dr2, Pu2, Pd2))

P2V = Problem2D(order,u₀,K,K,Dom2V,BC2V)
println("---Solving 1 volume---")
soln2V = solve(P2V,Dom2V,Δt,t)






#= Plotting =#


using GLMakie

f = Figure()

Ax = Axis3(f[1,1])
surface!(Ax,Dom1V.gridx,Dom1V.gridy,soln1V.u[2])

Ax2  = Axis3(f[1,2])
colourrange = (minimum(minimum.(soln2V.u[2])),maximum(maximum.(soln2V.u[2])))
surface!(Ax2,Dom2V.Grids[1].gridx,Dom2V.Grids[1].gridy,soln2V.u[2][1],colorrange=colourrange)
surface!(Ax2,Dom2V.Grids[2].gridx,Dom2V.Grids[2].gridy,soln2V.u[2][2],colorrange=colourrange)


f
