

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




g₀(t) = 0.0
g₁(t) = 1.0

𝒟 = [0.0,1.0]
n = 41
order = 2
K = ones(Float64,n)

Dom = Grid1D(𝒟, n)

# Define some boundary conditions
BoundaryLeft = Boundary(Dirichlet,g₀,Left,1)
BoundaryRight = Boundary(Dirichlet,g₁,Right,1)
# At this stage we have only told the solver what our boundary conditions are, the SATs are constructed after we build the PDE
P = VariableCoefficientPDE1D(Dom,K,order,BoundaryLeft,BoundaryRight)

Δt = 0.1Dom.Δx


BStor = SBP_operators.Helpers.BoundaryData1D{Float64}((Dirichlet,Dirichlet),n,order)
DStor = SBP_operators.Helpers.DataBlock{Float64}(P.grid,Δt,2,BoundaryLeft,BoundaryRight)






# SBP_operators.Helpers.DataBlock(P.grid,Δt,2,BoundaryLeft,BoundaryRight)

