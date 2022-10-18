# SBP_Operators.jl Documentation

```@docs
SBP_operators
```

```@contents
Pages = ["SATs.md","Derivatives.md","Helpers.md","solvers.md"]
Depth=2
```

# User interface and examples

```
push!("./path/to/SBP_operators")
using SBP_operators
```

## 1D

```
D = [0.0,1.0]
n = 41
grid = Grid1D(D,n)
```

```
BoundaryLeft = Boundary(Dirichlet,t -> 0.0,Left)
BoundaryRight = Boundary(Neumann,t -> 0.0,Right)
```

```
u0 = exp.(-(x.-0.5).^2 ./ 0.02)
```


```
Δt = 0.01*grid.Δx^2
t_f = 100Δt
```


```
P = VariableCoefficientPDE1D(u0,K,order,BoundaryLeft,BoundaryRight)
```

```
soln = solve(P,grid,Δt,t_f,:cgie)
```



```
using Plots
plot(soln.grid.grid,soln.u[2],xlims=(0.0,1.0),ylims=(0.0,1.0))
```


## 2D


