# Simultaneous Approximation Terms


Users can create a generic SBP operator with the function,
```@docs
FaADE.SATs.SAT
```
Note that this requires users to first create a boundary condition with [boundary conditions](@ref User_interaction)

## Boundary operators

The following are for boundary conditions
```@docs
FaADE.SATs.SAT_Dirichlet
FaADE.SATs.SAT_Neumann
```

If periodic boundaries are applied
```@docs
FaADE.SATs.SAT_Periodic
```

If a domain is split along some boundary then the following matches the values at the interface.




## Boundary derivatives and penalties

The `D_x^T B` and `B D_x` operators are given by,
```@docs
FaADE.SATs.BoundaryDerivativeTranspose
FaADE.SATs.BoundaryDerivative
```

Penalties are constructed with the following,
```@docs
FaADE.SATs.SATpenalties
FaADE.SATs.hval
```

