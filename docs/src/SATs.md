# Simultaneous Approximation Terms

```@docs
FaADE.SATs
```

All `SATs` belong to the `abstract type`

```@docs
FaADE.SATs.SimultanousApproximationTerm
```


## Boundary operators

The following are for boundary conditions
```@docs
FaADE.SATs.SAT_Dirichlet
FaADE.SATs.SAT_Neumann
```

Robin (mixed) boundary conditions are to be implemented still.

If periodic boundaries are applied
```@docs
FaADE.SATs.SAT_Periodic
```

Interface conditions are not currently implemented.



---

## Boundary derivatives and penalties

The following are all internally accessed when setting up the SATs

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


