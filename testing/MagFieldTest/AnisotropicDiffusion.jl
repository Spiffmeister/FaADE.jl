#=
Testing


∂ₜu = ∂ₓ(k ∂ₓu) = ∂ₓk ∂ₓu + k ∂ₓ∂ₓu
=#
using Distributed

addprocs(4)
@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere using SBP_operators

@everywhere using DifferentialEquations





# Initalise the things
n = 100

domain = [0,2π]

x = collect(range(domain[1],domain[2],length=n))

u = InitialCondition(x)

struct sim
    u   :: Matrix{Float64}  # Solution
    x   :: Vector{Float64}  # x grid
    y   :: Vector{Float64}  # y grid
    xy  :: Matrix{Float64}  # xyz grid
    Δx  :: Float64          # Spatial scale
    Δy  :: Float64          # Spatial scale
    Δz  :: Float64          # Spatial scale
    Δt  :: Float64          # Time step size
end







function magnetic_field(x,ϵₘₙ,m,n)
    # Magnetic field, assumes B(ψ,θ,ζ) = ∇×A(ψ,θ,ζ)
    B = zeros(Float64,3)
    B[1] = x[1]*(x[1]-1) * sum(ϵₘₙ .* sin.(m*x[2] - n*x[3]) .* m)
    B[2] = x[1] - 2*(x[1]-1) * sum(ϵₘₙ .* cos.(m*x[2] - n*x[3]))
    B[3] = 1.0
    return B
end

function field_line_hamiltonian(χ,x::Array{Float64},ϵₘₙ::Array{Float64},m::Array{Float64},n::Array{Float64},t::Float64)
    # The field line Hamiltonian for the above magnetic field
    χ[1] = x[2]
    χ[2] = -sum(ϵₘₙ .* sin.(m*x[1] - n*t) .* m)
    return χ
end





function timesolve()
end



#===
SET SIMULATION PARAMS
===#


# Domain
𝒟x = [0.5,0.68]
𝒟y = [0.0,2π]
# 𝒟z = [0.0,2π]

# Grid spacing
xₙ = 20
yₙ = 20
zₙ = 1

Δx = (𝒟x[2] - 𝒟x[1])/xₙ
Δy = (𝒟y[2] - 𝒟y[2])/yₙ


x = collect(range(𝒟x[1],𝒟x[2],step=Δx))
y = collect(range(𝒟y[1],𝒟y[2],step=Δy))



# Diffusion coefficient
κ = 1.0

# Perturbation parameters
k = 2.1e-3
ϵₘₙ = [k/2.0, k/3.0]
m = [2.0, 3.0]
n = [1.0, 2.0]


H(χ,x,t) = field_line_hamiltonian(χ,x,ϵₘₙ,m,n,t)
B(x) = magnetic_field(x,ϵₘₙ,m,n)

construct_grid(x,y,xₙ,yₙ,zₙ,Δx,Δy,H,B)





