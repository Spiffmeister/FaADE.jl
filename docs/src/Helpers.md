# Helpers

The `Helpers` sub-package contains both user interface function and functions useful for many of the internals such as data structure objects. 

## User interaction (@id User_interaction)

```@docs
SPADE.Helpers.Boundary
SPADE.Helpers.PeriodicBoundary
```

```@docs
SPADE.Helpers.VariableCoefficientPDE1D
SPADE.Helpers.VariableCoefficientPDE2D
```

```@docs
SPADE.Helpers.Grid1D
SPADE.Helpers.Grid2D
```


## Types

```@docs
SPADE.Helpers.NodeType
SPADE.Helpers.BoundaryConditionType
```




## Internal structures

These data blocks are used internally and are not interacted with by the user,
```@docs
SPADE.Helpers.DataBlock
SPADE.Helpers.BoundaryData1D
SPADE.Helpers.BoundaryData2D
```

Moving data between blocks uses the following functions,
```@docs
SPADE.Helpers.copyUtoSAT!
SPADE.Helpers.copySATtoU!
SPADE.Helpers.addSATtoU!
```


## Nifty internal functions
```@docs
SPADE.Helpers.check_order
SPADE.Helpers.check_boundary
SPADE.Helpers.halforder
SPADE.Helpers.BoundaryNodeInput
SPADE.Helpers.BoundaryNodeOutput
SPADE.Helpers.SATNodeOutput
SPADE.Helpers.SelectLoopDirection
```

```@docs
SPADE.Helpers.GetAxis
SPADE.Helpers.GetDim
SPADE.Helpers.GetMinΔ
```

