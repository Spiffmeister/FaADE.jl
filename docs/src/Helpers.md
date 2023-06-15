# Helpers

The `Helpers` sub-package contains both user interface function and functions useful for many of the internals such as data structure objects. 

## User interaction (@id User_interaction)

```@docs
FaADE.Helpers.Boundary
FaADE.Helpers.PeriodicBoundary
```

```@docs
FaADE.Helpers.VariableCoefficientPDE1D
FaADE.Helpers.VariableCoefficientPDE2D
```

```@docs
FaADE.Helpers.Grid1D
FaADE.Helpers.Grid2D
```


## Types

```@docs
FaADE.Helpers.NodeType
FaADE.Helpers.BoundaryConditionType
```




## Internal structures

These data blocks are used internally and are not interacted with by the user,
```@docs
FaADE.Helpers.DataBlock
FaADE.Helpers.BoundaryData1D
FaADE.Helpers.BoundaryData2D
```

Moving data between blocks uses the following functions,
```@docs
FaADE.Helpers.copyUtoSAT!
FaADE.Helpers.copySATtoU!
FaADE.Helpers.addSATtoU!
```


## Nifty internal functions
```@docs
FaADE.Helpers.check_order
FaADE.Helpers.check_boundary
FaADE.Helpers.halforder
FaADE.Helpers.BoundaryNodeInput
FaADE.Helpers.BoundaryNodeOutput
FaADE.Helpers.SATNodeOutput
FaADE.Helpers.SelectLoopDirection
```

```@docs
FaADE.Helpers.GetAxis
FaADE.Helpers.GetDim
FaADE.Helpers.GetMinΔ
```

