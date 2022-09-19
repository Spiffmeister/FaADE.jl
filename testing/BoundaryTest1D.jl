

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

push!(LOAD_PATH,".")

using SBP_operators





𝒟 = [0.0,1.0]
n = 41
Dom = Grid1D(𝒟, n)


g₀(t) = 0.0
g₁(t) = 1.0
order = 2
K = ones(Float64,n)
Δt = 0.1Dom.Δx

# Define some boundary conditions
BoundaryLeft = Boundary(Dirichlet,g₀,Left,1)
BoundaryRight = Boundary(Dirichlet,g₁,Right,1)

# Define initial condition
u₀(x) = exp.(-(x.-0.5).^2 ./ 0.02)
# Build PDE problem
P = VariableCoefficientPDE1D(u₀,K,order,BoundaryLeft,BoundaryRight)


# Testing internal data storage construction
BStor = SBP_operators.Helpers.BoundaryData1D{Float64}(P.BoundaryConditions,order)
DStor = SBP_operators.Helpers.DataBlock{Float64}(P.BoundaryConditions,Dom,Δt,2,P.K)
CGStor = SBP_operators.Helpers.ConjGradBlock{Float64}(n)



# Testing internal boundary function construction
SATD = SBP_operators.SATs.SATDirichlet(P.BoundaryConditions[1].RHS,Dom.Δx,P.BoundaryConditions[1].side,P.BoundaryConditions[1].axis,order)
SATD2, SATFn2 = SBP_operators.SATs.SAT(P.BoundaryConditions[1],Dom,order,:cgie)

# Solution storage
IC = u₀(Dom.grid)
solution = SBP_operators.solvers.solution{Float64}(IC,Dom,0.0,Δt,P)

