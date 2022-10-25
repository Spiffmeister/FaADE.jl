```@meta
EditURL = "<unknown>/../tutorials/exampleBField.jl"
```

# Magnetic field example

For this example we will solve the head equation with a mangetic field aligned with the grid

```math
  \mathbf{B} = (0,0,1)
```

```math
  \frac{\partial u}{\partial t} = (\kappa_\perp \nabla_\perp^2 + \kappa_\parallel \nabla_\parallel^2)u
```

````@example exampleBField
using SBP_operators
````

For this we'll solve the 2D heat equation
```math
      \frac{\partial u}{\partial t} = K\frac\Delta u
```
with boundary conditions
```math
      u(0,t) = 0, \qquad \partial_x u(1,t) = 0
```
and initial condition
```math
      u(x,0) = \exp\left(\frac{-(x-0.5)^2}{0.02}\right)
```


We first need to create a domain to solve the PDE using [`Grid2D`](@ref SBP_operators.Helpers.Grid2D)

````@example exampleBField
𝒟x = [0.0,1.0]
𝒟y = [0.0,1.0]
nx = ny = 41
grid = Grid2D(𝒟x,𝒟y,nx,ny)
````

The initial condition is a simple function

````@example exampleBField
u₀(x,y) = exp(-((x-0.5)^2 + (y-0.5)^2) / 0.02)
````

The boundary conditions are defined by creating [`Boundary`](@ref SBP_operators.Helpers.Boundary) objects, which will then be fed to the PDE structure

````@example exampleBField
BoundaryLeft = Boundary(Dirichlet,t->0.0,Left,1)
BoundaryRight = Boundary(Neumann,t->0.0,Right,1)
BoundaryUpDown = PeriodicBoundary(2)
````

Create a few more things we'll need for the PDE and the solver

````@example exampleBField
order = 2
method = :cgie

Kx = Ky = ones(Float64,nx,ny);

# The parallel term

H_x = SBP_operators.build_H(ny,order)
H_x = 1.0 ./H_x.^2

H_y = SBP_operators.build_H(nx,order)
H_y = 1.0 ./H_y.^2

function Ppar(u,uₒ,Δt)
    for j = 1:ny
        for i = 1:nx
            u[i,j] = 1.0/(1.0 - κ_para * τ_para/2.0 * Δt * (H_y[i] + H_x[j])) * (uₒ[i,j] - Δt*τ_para/4.0 *(H_y[i] + H_x[j])*([i,j] + [i,j]))
        end
    end
end
````

NOTE: currently only conjugate gradient implicit Euler (`:cgie`) works as a solver

Now we can create a PDE object to pass to the solver, in this case a [`VariableCoefficientPDE2D`](@ref SBP_operators.Helpers.VariableCoefficientPDE2D),

````@example exampleBField
P = VariableCoefficientPDE2D(u₀,Kx,Ky,order,BoundaryLeft,BoundaryRight,BoundaryUpDown)
````

Lastly before solving we define our time step and simulation time,

````@example exampleBField
Δt = 0.01grid.Δx;
t_f = 100Δt;
nothing #hide
````

Finally we call the solver (currently not working with `Documenter.jl`)

`soln = solve(P,grid,Δt,t_f,method,penalty_func=Ppar);`

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

