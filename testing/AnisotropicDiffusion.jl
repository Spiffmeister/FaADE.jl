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





function construct_grid(x,y,nx,ny,grid_fn;start=0.0,stop=2π)
    # Constructs a ψ,θ plane by tracing a given plane in a given direction from a given start point

    # X = zeros(Float64,ny,ny)
    # Y = zeros(Float64,nx,ny)
    plane = zeros(Float64,2,nx*ny)

    x₀ = [[y₀,x₀] for y₀ in y for x₀ in x]

    function prob_fn(prob,i,repeat)
        remake(prob,u₀=x₀[i])
    end

    P = ODEProblem(grid_fn,x₀[1],(start,stop))
    EP = EnsembleProblem(P,prob_func=prob_fn)

    sim = solve(EP,Tsit5(),EnsembleDistributed(),trajectories=nx*ny,batch_size=floor(Int64,nx*ny/nworkers()),save_on=false,save_end=true)

    for i = 1:length(sim.u)
        plane[:,i] = mod.(sim.u[i][2:-1:1,2],2π)
    end

    return plane
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


