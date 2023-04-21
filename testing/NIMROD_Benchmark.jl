using LinearAlgebra

cd("..")
using Interpolations
push!(LOAD_PATH,"./plas_diff")
push!(LOAD_PATH,"./SBP_operators")
using SBP_operators
using plas_diff




Ψ(x,y) = cos(π*x)*cos(π*y)

B(x,y,z) = [π*cos(π*x)*sin(π*y),
    π*sin(π*x)*cos(π*y),
    0]

# Source term
F(x,y) = -2π*cos(π*x)*cos(π*y)

# Domain
𝒟x = [-0.5,0.5]
𝒟y = [-0.5,0.5]
nx = ny = 32

Dom = Grid2D(𝒟x,𝒟y,nx,ny)


# Homogeneous boundary conditions
BoundaryLeft    = Boundary(Dirichlet,(y,t) -> 0.0, Left, 1)
BoundaryRight   = Boundary(Dirichlet,(y,t) -> 0.0, Right, 1)
BoundaryUp      = Boundary(Dirichlet,(x,t) -> 0.0, Up, 1)
BoundaryDown    = Boundary(Dirichlet,(x,t) -> 0.0, Down, 1)


# Initial condition
u₀(x,y) = 0.0

# Perpendicular diffusion coefficient
k(x,y) = 1.0

# Build PDE problem
P = VariableCoefficientPDE2D(u₀,k,k,2,BoundaryLeft,BoundaryRight,BoundaryUp,BoundaryDown)

# Time domain
Δt = 1.0Dom.Δx^2
t_f = Inf

# Solve
soln = solve(P,Dom,Δt,t_f,:cgie,adaptive=true)


# Exact solution
T(x,y,t) = (1.0 - exp(-2.0*k(x,y)*π^2*t) )/( k(x,y) )*Ψ(x,y)

T_exact = zeros(Dom.nx,Dom.ny)
for j = 1:ny
    for i = 1:nx
        T_exact[i,j] = T(Dom.gridx[i],Dom.gridy[j],soln.t)
    end
end


