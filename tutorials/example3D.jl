# # Magnetic field example
#
# For this example we will solve the head equation with a mangetic field aligned with the grid 
#
# ```math
#   \mathbf{B} = (0,0,1)
# ```
# In this case we expect the parallel operator to do nothing since $\mathbf{P}_f=\mathbf{P}_b=I$
#

using SBP_operators

#
# For this we'll solve the field aligned equation is
# ```math
#   \frac{\partial u}{\partial t} = \kappa_\perp \nabla_\perp^2 u + \mathcal{P}_\parallel u
# ```
# with Dirichlet boundaries in ``x``
# ```math
#       u(0,y,t) = 0, \qquad \partial_x u(1,y,t) = 0
# ```
# periodic in ``y``, and initial condition
# ```math
#       u(x,0) = \exp\left(\frac{-(x-0.5)^2}{0.02}\right)
# ```
#
# 
# We first need to create a domain to solve the PDE using [Grid2D](@ref SBP_operators.Helpers.Grid2D)
#

𝒟x = [0.0,1.0]
𝒟y = [0.0,1.0]
nx = ny = 21
grid = Grid2D(𝒟x,𝒟y,nx,ny)

# The initial condition is
u₀(x,y) = exp(-((x-0.5)^2 + (y-0.5)^2) / 0.02)

# The boundary conditions are defined by creating [Boundary](@ref SBP_operators.Helpers.Boundary) objects, which will then be fed to the PDE structure
BoundaryLeft = Boundary(Dirichlet,t->0.0,Left)
BoundaryRight = Boundary(Neumann,t->0.0,Right)
BoundaryUpDown = PeriodicBoundary(2)

# The `2` input to the periodic boundary ensures it is along the y-axis.
#
# Set the FD order to 2 and use the conjugate gradient implicit euler (`:cgie`) solver,

order = 2
method = :cgie

# Forward Euler and RK4 are also available.
#
# Set the diffusion in $x$ and $y$ directions to 1

Kx(y) = 1.0
Ky(x) = 1.0

# NOTE: currently only conjugate gradient implicit Euler (`:cgie`) works as a solver
#
# Now we can create a PDE object to pass to the solver, in this case a [VariableCoefficientPDE2D](@ref SBP_operators.Helpers.VariableCoefficientPDE2D),

P = VariableCoefficientPDE2D(u₀,Kx,Ky,order,BoundaryLeft,BoundaryRight,BoundaryUpDown)

# The parallel penalty function can be generated by providing the code the ODE for magnetic field lines,

function Bfield(X,x,p,t)
    X[2] = 0.0
    X[1] = 0.0
end

# Assuming a $2\pi$ periodicity then we can construct a parallel grid object with [construct_grid](@ref)

PGrid = construct_grid(Bfield,Dom,[-2π,2π])

# Lastly before solving we define our time step and simulation time,

Δt = 0.01grid.Δx;
t_f = 100Δt;

# Finally we call the solver (currently not working with `Documenter.jl`)
# 
# `soln = solve(P,grid,Δt,t_f,method,penalty_func=Ppar);`
