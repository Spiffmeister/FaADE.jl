# Derivative operators

Below are the first and second derivative summation by parts operators. All stencils are written in 1D and then extended to two dimensions by looping over the required axis.

```@docs
FaADE.Derivatives
```




Derivative operators with structs belong to the `abstract type`.

```@docs
FaADE.Derivatives.DerivativeOperatorType
```




# First derivative operators

The below operators are for first derivatives

```@autodocs
Modules = [FaADE.Derivatives]
Pages = ["Op_Dx.jl"]
```


# Second Derivative operators

```@autodocs
Modules = [FaADE.Derivatives]
Pages = ["Op_Dxx.jl"]
```


# Internals

## First derivative internals

```@autodocs
Modules = [FaADE.Derivatives]
Pages = ["DerivativeFirst.jl"]
```


## Second Derivative internals

```@autodocs
Modules = [FaADE.Derivatives]
Pages = ["DerivativeSecond.jl"]
```

