#=
Testing


∂ₜu = ∂ₓ(k ∂ₓu) = ∂ₓk ∂ₓu + k ∂ₓ∂ₓu
=#


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





function construct_grid(x,y,mnk,grid_fn;symmetry=true)
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



# Diffusion coefficient
κ = 1.0

# Perturbation parameters
k = 2.1e-3
ϵₘₙ = [k/2.0, k/3.0]
m = [2.0, 3.0]
n = [1.0, 2.0]


