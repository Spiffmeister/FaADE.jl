# # 2D example
# 
# The 2D code works similarly to the 1D version, with a few different function calls.

using FaADE

#
# This example gives a 1D solution to a 2D problem.
#
# We'll solve the 2D heat equation.
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
# We first need to create a domain to solve the PDE using [`Grid2D`](@ref FaADE.Helpers.Grid2D)
#

𝒟x = [0.0,1.0]
𝒟y = [0.0,1.0]
nx = ny = 41
grid = Grid2D(𝒟x,𝒟y,nx,ny)

# The initial condition
u₀(x,y) = exp(-((x-0.5)^2 + (y-0.5)^2) / 0.02)

# The boundary conditions are defined by creating [`Boundary`](@ref FaADE.Helpers.Boundary) objects, which will then be fed to the PDE structure
BoundaryLeft = Boundary(Dirichlet,(y,t)->0.0,Left)
BoundaryRight = Boundary(Neumann,(y,t)->0.0,Right)
BoundaryUpDown = PeriodicBoundary(2)

# The `2` input to the periodic boundary ensures it is along the y-axis.
#
# Set the FD order to 2 and use the conjugate gradient implicit euler (`:cgie`) solver,

order = 2
method = :cgie

# Forward Euler and RK4 are also available.
#
# Set the diffusion in $x$ and $y$ directions to 1

Kx(x,y) = 1.0
Ky(x,y) = 1.0

# Now we can create a PDE object to pass to the solver, in this case a [`VariableCoefficientPDE2D`](@ref FaADE.Helpers.VariableCoefficientPDE2D),

P = VariableCoefficientPDE2D(u₀,Kx,Ky,order,BoundaryLeft,BoundaryRight,BoundaryUpDown)

# Lastly before solving we define our time step and simulation time,

Δt = 0.01grid.Δx;
t_f = 100Δt;

# Finally we call the solver (currently not working with `Documenter.jl`)
# 
soln = solve(P,grid,Δt,t_f,method)

#
# The solver ourputs a [`solution`](@ref FaADE.solvers.solution) data structure, with everything packaged in that we would need to reconstruct
# the problem from the final state if we wanted to restart.
# 
# No visualisation routines are written at the moment but we imported the `Plots.jl` package earlier so we'll use that


