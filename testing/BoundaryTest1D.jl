

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



#### 1D
𝒟 = [0.0,1.0]
n = 41
order = 2
K = ones(Float64,n)


D = Grid1D(𝒟, n)


# Define some boundary conditions
BoundaryLeft = Boundary(Dirichlet,g₀,Left,1)
BoundaryRight = Boundary(Dirichlet,g₁,Right,1)
# At this stage we have only told the solver what our boundary conditions are, the SATs are constructed after we build the PDE
P = VariableCoefficientPDE1D(D,K,order,(BoundaryLeft,BoundaryRight))







