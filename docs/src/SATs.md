# Simultaneous Approximation Terms

```@docs
SBP_operators.SATs
```


Users can use the generic SBP operator
```@docs
SBP_operators.SATs.SAT
```

## Boundary operators

The following are for boundary conditions on either side of the domain.
```@docs
SBP_operators.SATs.SAT_Dirichlet
SBP_operators.SATs.SAT_Neumann
```

If periodic boundaries are applied to both sides of the domain.
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