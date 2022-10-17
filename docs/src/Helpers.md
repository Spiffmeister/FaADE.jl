# Helpers

```@docs
SBP_operators.Helpers
```

## User interaction

```@docs
SBP_operators.Helpers.Boundary
SBP_operators.Helpers.PeriodicBoundary
```

```@docs
SBP_operators.Helpers.VariableCoefficientPDE1D
SBP_operators.Helpers.VariableCoefficientPDE2D
```

```@docs
SBP_operators.Helpers.Grid1D
SBP_operators.Helpers.Grid2D
```


## Types

```@docs
SBP_operators.Helpers.NodeType
SBP_operators.Helpers.BoundaryConditionType
```




## Internal structures

These data blocks are used internally and are not interacted with by the user,
```@docs
SBP_operators.Helpers.DataBlock
SBP_operators.Helpers.ConjGradBlock
SBP_operators.Helpers.BoundaryData1D
SBP_operators.Helpers.BoundaryData2D
```

Moving data between blocks uses the following functions,
```@docs
SBP_operators.Helpers.copyUtoSAT!
SBP_operators.Helpers.copySATtoU!
SBP_operators.Helpers.addSATtoU!
```


## Nifty internal functions
```@docs
SBP_operators.Helpers.check_order
SBP_operators.Helpers.check_boundary
SBP_operators.Helpers.halforder
SBP_operators.Helpers.BoundaryNodeInput
SBP_operators.Helpers.BoundaryNodeOutput
SBP_operators.Helpers.SATNodeOutput
SBP_operators.Helpers.SelectLoopDirection
```

```@docs
SBP_operators.Helpers.GetAxis
SBP_operators.Helpers.GetDim
SBP_operators.Helpers.GetMinΔ
```

