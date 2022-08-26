



struct grid
    x   :: Array{Float64}
    Δx  :: Float64
    n   :: Int64

    function grid(𝒟,n)
        Δx = (𝒟[2] - 𝒟[1])/(n-1)
        x = collect(range(𝒟[1],𝒟[2],step=Δx))
        new(x,Δx,n)
    end
end



mutable struct solution
    u       :: Vector{Vector{Float64}}
    x       :: Vector{Float64}
    Δt      :: Union{Float64,Vector{Float64}}
    t       :: Vector{Float64}

    function solution(u₀,x,Δx,t,Δt;preallocate=false)
        if preallocate
            N = ceil(Int64,t/Δt)
            n = length(x)
            u = [zeros(Float64,n) for _ in 1:N]

            u[1] = u₀

            new(u,x,Δt,collect(range(0.0,t,length=N)))
        else
            u = u₀
            new([u],x,[Δt],[t])
        end

    end
end

mutable struct solution_2d
    u   :: Vector{Matrix{Float64}}
    x   :: Vector{Float64}
    y   :: Vector{Float64}
    Δt  :: Vector{Float64}
    t   :: Vector{Float64}
    function solution_2d(u₀,x,y,t,Δt)
        new([u₀],x,y,[Δt],[t])
    end
end


"""
    Struct for storing checkpoints for 2D simulations
"""
mutable struct checkpoint_2d
    # Solution info
    soln        :: solution_2d
    Δx          :: Float64
    Δy          :: Float64
    t_f         :: Float64
    # Diffusion coefficient matricies
    kx          :: Matrix{Float64}
    ky          :: Matrix{Float64}
    # Boundary functions
    gx          :: Function
    gy          :: Function
    # Parallel penalty function if provided
    parallel    :: Bool
    penalty_fn  :: Union{Function,Nothing}
    # Simulation parameters
    order_x     :: Int64
    order_y     :: Int64
    method      :: Symbol
    maxIT       :: Int64
    samplefactor:: Float64
    tol         :: Float64
    rtol        :: Float64
    adaptive    :: Bool
end

