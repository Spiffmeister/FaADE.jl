# # 2D example
# 
# The 2D code works similarly to the 1D version, with a few different function calls.

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
# We first need to create a domain to solve the PDE using [`Grid2D`](@ref SBP_operators.Helpers.Grid2D)
#

𝒟x = [0.0,1.0]
𝒟y = [0.0,1.0]
nx = ny = 41
grid = Grid2D(𝒟x,𝒟y,nx,ny)

# The initial condition is a simple function
u₀(x,y) = exp(-((x-0.5)^2 + (y-0.5)^2) / 0.02)

# The boundary conditions are defined by creating [`Boundary`](@ref SBP_operators.Helpers.Boundary) objects, which will then be fed to the PDE structure
BoundaryLeft = Boundary(Dirichlet,t->0.0,Left,1)
BoundaryRight = Boundary(Neumann,t->0.0,Right,1)
BoundaryUpDown = PeriodicBoundary(2)

# Create a few more things we'll need for the PDE and the solver

order = 2
method = :cgie

Kx = Ky = ones(Float64,nx,ny);

# NOTE: currently only conjugate gradient implicit Euler (`:cgie`) works as a solver
#
# Now we can create a PDE object to pass to the solver, in this case a [`VariableCoefficientPDE2D`](@ref SBP_operators.Helpers.VariableCoefficientPDE2D),

P = VariableCoefficientPDE2D(u₀,Kx,Ky,order,BoundaryLeft,BoundaryRight,BoundaryUpDown)

# Lastly before solving we define our time step and simulation time,

Δt = 0.01grid.Δx;
t_f = 100Δt;

# Finally we call the solver (currently not working with `Documenter.jl`)
# 
# `soln = solve(P,grid,Δt,t_f,method);`

#
# The solver ourputs a [`solution`](@ref SBP_operators.solvers.solution) data structure, with everything packaged in that we would need to reconstruct
# the problem from the final state if we wanted to restart.
# 
# No visualisation routines are written at the moment but we imported the `Plots.jl` package earlier so we'll use that

# `using Plots`
# `plot(soln.grid.grid,soln.u[2])`

