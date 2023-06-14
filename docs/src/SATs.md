# Simultaneous Approximation Terms

```@docs
SBP_operators.SATs
```

Users can create a generic SBP operator with the function,
```@docs
SBP_operators.SATs.SAT
```
Note that this requires users to first create a boundary condition with [`BoundaryConditions`](@ref)

## Boundary operators

The following are for boundary conditions
```@docs
SBP_operators.SATs.SAT_Dirichlet
SBP_operators.SATs.SAT_Neumann
```

If periodic boundaries are applied
```@docs
SBP_operators.SATs.SAT_Periodic
```

If a domain is split along some boundary then the following matches the values at the interface.




## Boundary derivatives and penalties

The `D_x^T B` and `B D_x` operators are given by,
```@docs
SBP_operators.SATs.BoundaryDerivativeTranspose
SBP_operators.SATs.BoundaryDerivative
```

Penalties are constructed with the following,
```@docs
SBP_operators.SATs.SATpenalties
SBP_operators.SATs.hval
```

