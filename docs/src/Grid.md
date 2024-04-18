# Grid.md

```@docs
FaADE.Grid
```

One thing to note if editing code is that the `Grid.jl` module has its own internal copy of first derivative operators for constructing curvilinear coordinates.

## 1D Problems

```@docs
Grid1D
```



## 2D

```@docs
Grid2D
```


## Multiblock grids

**Solvers for these methods are not implemented yet (SATs required)**

The following is used for constructing multiblock problems. 


```@docs
GridMultiBlock
```




## Curvilinear grid generation

The `meshgrid` function can be used to generate domains in curvilinear coordinates.

```@docs
FaADE.Grid.meshgrid
```

### Toroidal surfaces

For generating meshgrid functions in an annulus one can use the following commands to generate two torii.

```@docs
FaADE.Grid.Torus
```



