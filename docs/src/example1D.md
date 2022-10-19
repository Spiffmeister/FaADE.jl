```@meta
EditURL = "<unknown>/../tutorials/example1D.jl"
```

# 1D example

````@example example1D
using SBP_operators
````

As an example we'll solve the 1D heat equation
```math
      \frac{\partial u}{\partial t} = K\frac{\partial}{\partial x}\frac{\partial u}{\partial x}
```
with boundary conditions
```math
      u(0,t) = 0, \qquad \partial_x u(1,t) = 0
```
and initial condition
```math
      u(x,0) = \exp\left(\frac{-(x-0.5)^2}{0.02}\right)
```


We first need to create a domain to solve the PDE,
`Grid1D` ([link](@ref Grid1D))

````@example example1D
𝒟 = [0.0,1.0]
n = 41
grid = Grid1D(𝒟,n)
````

The initial condition is a simple function

````@example example1D
u₀(x) = exp.(-(x.-0.5).^2 ./ 0.02)
````

The boundary conditions are defined by creating `Boundary` objects, which will then be fed to the PDE structure

````@example example1D
BoundaryLeft = Boundary(Dirichlet,t->0.0,Left,1)
BoundaryRight = Boundary(Neumann,t->0.0,Right,1)
````

Create a few more things we'll need for the PDE and the solver

````@example example1D
order = 2
method = :cgie

K = ones(Float64,n);
nothing #hide
````

NOTE: currently only conjugate gradient implicit Euler (`:cgie`) works as a solver

Now we can create a PDE object to pass to the solver, in this case a `VariableCoefficientPDE1D` ([link](@ref VariableCoefficientPDE1D)),

````@example example1D
P = VariableCoefficientPDE1D(u₀,K,order,BoundaryLeft,BoundaryRight)
````

Lastly before solving we define our time step and simulation time,

````@example example1D
Δt = 0.01grid.Δx;
t_f = 100Δt;
nothing #hide
````

Finally we call the solver (currently not working)

`soln = solve(P,grid,Δt,t_f,method);`

The solver ourputs a [`solution`](@ref solution) data structure, with everything packaged in that we would need to reconstruct
the problem from the final state if we wanted to restart.

No visualisation routines are written at the moment but we imported the `Plots.jl` package earlier so we'll use that

using Plots
plot(soln.grid.grid,soln.u[2])

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

