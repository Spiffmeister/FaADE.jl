# Summation by Parts operators

Below are the first and second derivative summation by parts operators. All stencils are written in 1D and then extended to two dimensions by looping over the required axis.


## First Derivative operators

The below operators are for first derivatives

```@docs
FaADE.Derivatives.D₁
```



### In place
```@docs
FaADE.Derivatives.D₁!
FaADE.Derivatives.D₁ᵀ!
```



## Second Derivative operators


```@docs
FaADE.Derivatives.D₂
```


---

### In place

The following functions are also available as iterators

```@docs
FaADE.Derivatives.D₂!
```




## Internal Functions

The following are all used internally for the various stencils.


```@autodocs
Modules = [FaADE.Derivatives]
Pages = ["DerivativeFirst.jl","DerivativeSecond.jl"]
```

