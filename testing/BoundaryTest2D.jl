

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

push!(LOAD_PATH,"."); using SBP_operators



𝒟x = [0.0,1.0]
𝒟y = [0.0,1.0]
nx = 41
ny = 41

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
BStor = SBP_operators.Helpers.BoundaryData2D{Float64}(PD.BoundaryConditions,Dom,order)
DStor = SBP_operators.Helpers.DataBlock{Float64}(PD.BoundaryConditions,Dom,Δt,2,PD.Kx,PD.Ky)
CGStor = SBP_operators.Helpers.ConjGradBlock{Float64}(nx,ny)


