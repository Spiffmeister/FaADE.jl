

#====
Testing the boundary operator construction

For Problem setup:
    1. Select PDE type (only handles variable coefficient diffusion problem at the moment)
        a. Order
        b. Domain
        c. Grid

For Dirichlet, Neumann, Robin the workflow is:
    1. Set 
====#

push!(LOAD_PATH,"."); using SPADE

𝒟x = [0.0,1.0]
𝒟y = [0.0,1.0]
nx = 11
ny = 16

Dom = Grid2D(𝒟x,𝒟y, nx,ny)

order = 2
kx = zeros(Float64,nx,ny) .+ 1.0;
ky = zeros(Float64,nx,ny) .+ 1.0;
Δt = 0.1Dom.Δx

# Define initial condition
u₀(x,y) = exp(-((x-0.5)^2 + (y-0.5)^2) / 0.02)


# Define some boundary conditions
BoundaryDirichletLeft   = Boundary(Dirichlet,t->0.0,Left,1)
BoundaryDirichletRight  = Boundary(Dirichlet,t->1.0,Right,1)
BoundaryDirichletUp     = Boundary(Dirichlet,t->0.0,Up,1)
BoundaryDirichletDown   = Boundary(Dirichlet,t->1.0,Down,1)

# Build PDE problem
PD = VariableCoefficientPDE2D(u₀,kx,ky,order,BoundaryDirichletLeft,BoundaryDirichletRight,BoundaryDirichletUp,BoundaryDirichletDown)

# Testing internal data storage construction
BStor = SPADE.Helpers.BoundaryData2D{Float64}(PD.BoundaryConditions,Dom,order)
DStor = SPADE.Helpers.DataBlock{Float64}(PD.BoundaryConditions,Dom,Δt,2,PD.Kx,PD.Ky)
CGStor = SPADE.Helpers.ConjGradBlock{Float64}(Dom,order)






SATDL = SPADE.SATs.SAT_Dirichlet(PD.BoundaryConditions[1].RHS,Dom.Δx,PD.BoundaryConditions[1].side,PD.BoundaryConditions[1].axis,order)
SATDR = SPADE.SATs.SAT_Dirichlet(PD.BoundaryConditions[2].RHS,Dom.Δx,PD.BoundaryConditions[2].side,PD.BoundaryConditions[2].axis,order)
SATDU = SPADE.SATs.SAT_Dirichlet(PD.BoundaryConditions[3].RHS,Dom.Δy,PD.BoundaryConditions[3].side,PD.BoundaryConditions[3].axis,order)
SATDD = SPADE.SATs.SAT_Dirichlet(PD.BoundaryConditions[4].RHS,Dom.Δy,PD.BoundaryConditions[4].side,PD.BoundaryConditions[4].axis,order)

SATD1, SATFn1 = SPADE.SATs.SAT(PD.BoundaryConditions[1],Dom,order,:cgie)
SATD2, SATFn2 = SPADE.SATs.SAT(PD.BoundaryConditions[2],Dom,order,:cgie)
SATD3, SATFn3 = SPADE.SATs.SAT(PD.BoundaryConditions[3],Dom,order,:cgie)
SATD4, SATFn4 = SPADE.SATs.SAT(PD.BoundaryConditions[4],Dom,order,:cgie)



#= PERIODIC BOUNDARIES =#


BoundaryPeriodicLeftRight = PeriodicBoundary(1)
BoundaryPeriodicUpDown = PeriodicBoundary(2)

PP = VariableCoefficientPDE2D(u₀,kx,ky,order,BoundaryPeriodicLeftRight,BoundaryPeriodicUpDown)

DStor = SPADE.Helpers.DataBlock{Float64}(PP.BoundaryConditions,Dom,Δt,2,PP.Kx,PP.Ky)
BStor = SPADE.Helpers.BoundaryData2D{Float64}(PP.BoundaryConditions,Dom,order)



