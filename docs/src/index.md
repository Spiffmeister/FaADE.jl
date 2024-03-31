# FaADE.jl

_A Summation by Parts code for solving the field aligned Anisotropic Diffusion Equation._ 


## Features

`FaADE.jl` is a code for solving the field aligned anisotropic diffusion equation

$$\frac{\partial u}{\partial t} = \nabla\cdot(\mathbf{K} \nabla ) u$$

where

$$\mathbf{K} = k_\perp\mathbf{I} + (k_\parallel - k_\perp)\frac{\mathbf{B}\mathbf{B}^T}{\|\mathbf{B}\|^2}.$$


- Uses the Summation by Parts formulation with Simultaneous Approximation Terms (SBP-SAT) for boundary conditions [[1](https://doi.org/10.1007/s10915-011-9525-z)].
- Parallel penalty operator used to apply diffusion along vector field lines.
- Currently arbitrary parallel mappings can be provided in Cartesian geometry or an ODE for mapping grid points along field lines.
- Provides solutions in 1D for:
    - diffusion problems
    - with a parallel map in the second dimension
- and solutions in 2D for:
    - diffusion problems
    - with a parallel map in the third dimension



## Examples

The best place to start is in the Examples section in the navigation bar to the left.

## Modules


```@contents
Pages = ["SATs.md","Derivatives.md","Helpers.md","solvers.md"]
Depth=2
```

## Similar software

- [SummationByParts.jl](https://github.com/ranocha/SummationByPartsOperators.jl): A Julia implementation of a wide range of SBP operators
- [pyoculus](https://github.com/zhisong/pyoculus): A magnetic field diagnostic package in python based on an earlier FORTRAN implementation [oculus](https://github.com/SRHudson/Oculus)




## References

The mathematical background for this package can be found in:
- D. Muir, K. Duru, M. Hole, and S. Hudson, ‘Provably stable numerical method for the anisotropic diffusion equation in toroidally confined magnetic fields’. arXiv, Jun. 01, 2023. doi: [10.48550/arXiv.2306.00423](http://arxiv.org/abs/2306.00423)




