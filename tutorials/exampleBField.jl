# # Magnetic field example
#
# For this example we will solve the head equation with a mangetic field aligned with the grid 
#
# ```math
#   \mathbf{B} = (0,0,1)
# ```
# 
# ```math
#   \frac{\partial u}{\partial t} = K\nabla_\perp u
# ```
#


using SBP_operators

#
# For this we'll solve the 2D heat equation
# ```math
#       \frac{\partial u}{\partial t} = K\frac\Delta u
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
# We first need to create a domain to solve the PDE,
# `Grid1D` ([link](@ref Grid1D))
#

𝒟x = [0.0,1.0]
𝒟y = [0.0,1.0]
nx = ny = 41
grid = Grid2D(𝒟x,𝒟y,nx,ny)

# The initial condition is a simple function
u₀(x,y) = exp(-((x-0.5)^2 + (y-0.5)^2) / 0.02)

# The boundary conditions are defined by creating `Boundary` objects, which will then be fed to the PDE structure
BoundaryLeft = Boundary(Dirichlet,t->0.0,Left,1)
BoundaryRight = Boundary(Neumann,t->0.0,Right,1)
BoundaryUpDown = PeriodicBoundary(2)

# Create a few more things we'll need for the PDE and the solver

order = 2
method = :cgie

Kx = Ky = ones(Float64,nx,ny);

## The parallel term

function P∥(u,uₒ,)
end


# NOTE: currently only conjugate gradient implicit Euler (`:cgie`) works as a solver
#
# Now we can create a PDE object to pass to the solver, in this case a `VariableCoefficientPDE1D` ([link](@ref VariableCoefficientPDE1D)),

P = VariableCoefficientPDE2D(u₀,Kx,Ky,order,BoundaryLeft,BoundaryRight,BoundaryUpDown)

# Lastly before solving we define our time step and simulation time,

Δt = 0.01grid.Δx;
t_f = 100Δt;

# Finally we call the solver (currently not working)
# 
# `soln = solve(P,grid,Δt,t_f,method);`

#