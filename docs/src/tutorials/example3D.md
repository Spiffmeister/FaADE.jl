```@meta
EditURL = "<unknown>/tutorials/example3D.jl"
```

# 3D Example

For this example we will solve the head equation with a mangetic field aligned with the grid

```math
  \mathbf{B} = (0,0,1)
```
In this case we expect the parallel operator to do nothing since ``\mathbf{P}_f=\mathbf{P}_b=I``

````@example example3D
using SPADE
````

For this we'll solve the field aligned equation is
```math
  \frac{\partial u}{\partial t} = \kappa_\perp \nabla_\perp^2 u + \mathcal{P}_\parallel u
```
with Dirichlet boundaries in ``x``
```math
      u(0,y,t) = 0, \qquad \partial_x u(1,y,t) = 0
```
periodic in ``y``, and initial condition
```math
      u(x,0) = \exp\left(\frac{-(x-0.5)^2}{0.02}\right)
```


We first need to create a domain to solve the PDE using [`Grid2D`](@ref SPADE.Helpers.Grid2D)

````@example example3D
𝒟x = [0.0,1.0]
𝒟y = [0.0,1.0]
nx = ny = 21
grid = Grid2D(𝒟x,𝒟y,nx,ny)
````

The initial condition is

````@example example3D
u₀(x,y) = exp(-((x-0.5)^2 + (y-0.5)^2) / 0.02)
````

The boundary conditions are defined by creating [`Boundary`](@ref SPADE.Helpers.Boundary) objects, which will then be fed to the PDE structure

````@example example3D
BoundaryLeft = Boundary(Dirichlet,(y,t)->0.0,Left)
BoundaryRight = Boundary(Neumann,(y,t)->0.0,Right)
BoundaryUpDown = PeriodicBoundary(2)
````

The `2` input to the periodic boundary ensures it is along the y-axis.

Set the FD order to 2 and use the conjugate gradient implicit euler (`:cgie`) solver,

````@example example3D
order = 2
method = :cgie
````

Forward Euler and RK4 are also available.

Set the diffusion in $x$ and $y$ directions to 1

````@example example3D
Kx(x,y) = 1.0
Ky(x,y) = 1.0
````

NOTE: currently only conjugate gradient implicit Euler (`:cgie`) works as a solver

Now we can create a PDE object to pass to the solver, in this case a [`VariableCoefficientPDE2D`](@ref SPADE.Helpers.VariableCoefficientPDE2D),

````@example example3D
P = VariableCoefficientPDE2D(u₀,Kx,Ky,order,BoundaryLeft,BoundaryRight,BoundaryUpDown)
````

The parallel penalty function can be generated by providing the code the ODE for magnetic field lines,

````@example example3D
function Bfield(X,x,p,t)
    X[2] = 0.0
    X[1] = 0.0
end
````

Assuming a ``2\pi`` periodicity then we can construct a parallel grid object with [`construct_grid`](@ref SPADE.Parallel.construct_grid)

````@example example3D
PGrid = construct_grid(Bfield,grid,[-2π,2π])
````

It is not strictly necessary (since you can simply provide the `PGrid` to the `solve` call), but we can construct the penalty function manually with

````@example example3D
Pfn = generate_parallel_penalty(PGrid,grid,order)
````

Lastly before solving we define our time step and simulation time,

````@example example3D
Δt = 0.01grid.Δx;
t_f = 100Δt;
nothing #hide
````

Finally we call the solver

````@example example3D
soln = solve(P,grid,Δt,t_f,method,penalty_func=Pfn)
````

To provide the parallel grid instead use `soln = solve(P,grid,Δt,t_f,method,Pgrid=PGrid)`

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*
