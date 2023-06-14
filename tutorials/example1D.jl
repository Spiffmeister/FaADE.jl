# # 1D example
# 
#

using SPADE

#
# As an example we'll solve the 1D heat equation with no parallel mapping. The PDE is
# ```math
#       \frac{\partial u}{\partial t} = k\frac{\partial}{\partial x}\frac{\partial u}{\partial x}
# ```
# with boundary conditions
# ```math
#       u(0,t) = 0, \qquad \partial_x u(1,t) = 0
# ```
# and initial condition
# ```math
#       u(x,0) = \exp\left(\frac{-(x-0.5)^2}{0.02}\right)
# ```
#
# 
# We first create a domain to solve the PDE using [`Grid1D`](@ref SPADE.Helpers.Grid1D),
# 

𝒟 = [0.0,1.0]
n = 41
grid = Grid1D(𝒟,n)

# The initial condition is a simple Gaussian
u₀(x) = exp.(-(x.-0.5).^2 ./ 0.02)

# The boundary conditions are defined by creating [`Boundary`](@ref SPADE.Helpers.Boundary) objects, which will then be fed to the PDE structure
BoundaryLeft = Boundary(Dirichlet,t->0.0,Left)
BoundaryRight = Boundary(Neumann,t->0.0,Right,1)

# Create a few more things we'll need for the PDE and the solver such as the order (2) and the solver (conjugate gradient implicit euler)

order = 2;
method = :cgie;

# We will set the diffusion coefficient to 1 eveywhere in the domain

K(x) = 1.0

# Now we can create a PDE object to pass to the solver, in this case a [`VariableCoefficientPDE1D`](@ref SPADE.Helpers.VariableCoefficientPDE1D),

P = VariableCoefficientPDE1D(u₀,K,order,BoundaryLeft,BoundaryRight)

# Lastly before solving we define our time step and simulation time,

Δt = 0.01grid.Δx;
t_f = 100Δt;

# Finally we call the solver (currently not working with `Documenter.jl`)
# 
soln = solve(P,grid,Δt,t_f,method)

#
# The solver outputs a [`solution`](@ref SPADE.solvers.solution) data structure, with everything packaged in that we would need to reconstruct
# the problem from the final state if we wanted to restart.
# 
# No visualisation routines are written at the moment, coming soon.
