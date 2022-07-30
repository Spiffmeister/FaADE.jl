# Simultaneous Approximation Terms

Users can use the generic SBP operator
```@docs
SBP_operators.SAT
```

## Boundary operators

The following are for boundary conditions on either side of the domain.
```@docs
SBP_operators.SAT_Dirichlet
SBP_operators.SAT_Neumann
SBP_operators.SAT_Robin
```

If periodic boundaries are applied to both sides of the domain.
```@docs
SBP_operators.SAT_Periodic
```

If a domain is split along some boundary then the following matches the values at the interface.

```@docs
SBP_operators.Split_domain
```
