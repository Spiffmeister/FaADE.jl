
struct Boundary_Periodic <: SimultanousApproximationTerm
    BDₓᵀ        :: Vector{Real}
    B₀Dₓ        :: Vector{Real}
    BₙDₓ        :: Vector{Real}
    type        :: BoundaryConditionType
    axis        :: Int
    Δx          :: Real
    penalties   :: NamedTuple

    function Boundary_Periodic(Δx::Real,axis::Int,order::Int)

        BDₓᵀ = BoundaryDerivativeTranspose(order,Δx)
        B₀Dₓ = BoundaryDerivative(Left,order,Δx)
        BₙDₓ = BoundaryDerivative(Right,order,Δx)

        α₀, τ₁, τ₀ = SATpenalties(Periodic,Δx,order)
        penalties = (α₀=α₀, τ₁=τ₁, τ₀=τ₀)

        new(BDₓᵀ,B₀Dₓ,BₙDₓ,Periodic,axis,Δx,penalties)
    end
end




"""
    generate_Periodic
"""
function generate_Periodic(SATP::Boundary_Periodic,solver)

    loopdirection = select_SAT_direction(SATP.axis)

    p = SATP.penalties

    Term(cache,u,c) = 
        SAT_Periodic!(cache,u,c,p.α₀,p.τ₀,p.τ₁,SATP.B₀Dₓ,SATP.BₙDₓ,SATP.BDₓᵀ,SATP.order,loopdirection)
    
    return Term
end







"""
    SAT_Periodic(u::Vector{Float64},Δx::Float64;c=[1.0,1.0],order::Int64=2)

Simulatenous approximation term for Periodic boundary conditions
"""
function SAT_Periodic(u::AbstractVector,Δx::Float64;
        c=[1.0,1.0],order::Int64=2)

    # Get h value
    h = hval(order)

    # Penalties
    α₀ = 0.5/(h * Δx) # Derivative penatly
    τ₁ = -0.5/(h * Δx) # Symmeteriser penalty
    τ₀ = -max(c[1]/2(h*Δx),c[end]/2(h*Δx))/(h*Δx) # Dirichlet penalty

    SAT = zeros(Float64,2order)

    # Dirichlet terms
    SAT[1]  += τ₀ * (u[1] - u[end])
    SAT[end]+= τ₀ * (u[end] - u[1])

    # Symmeteriser
    L₁u = [(u[1] - u[end])]
    SAT[1:order]        = τ₁ * c[1] * BDₓᵀ(L₁u,Left,Δx,order)
    SAT[order+1:2order] = -τ₁ * c[end] * BDₓᵀ(L₁u,Right,Δx,order)

    # Neumann term
    BDₓu₁ = BDₓ(u,Left,Δx,order)
    BDₓuₙ = BDₓ(u,Right,Δx,order)

    SAT[1]  += α₀ * (c[1]*BDₓu₁ - c[end]*BDₓuₙ)
    SAT[end]+= α₀ * (c[1]*BDₓu₁ - c[end]*BDₓuₙ)

    return SAT[1:order], SAT[order+1:end]

end


function SAT_Periodic!(SAT::AbstractArray,u::AbstractArray,Δ::Float64;
        c=[1.0,1.0], order::Int=2)
    map((A,U,K) -> SAT_Periodic!(A,U,Δ,c=K,order=order), eachrow(SAT),eachrow(u),eachrow(c))
    SAT
end


function SAT_Periodic!(uₓₓ::AbstractVector,u::AbstractVector,Δx::Float64;
        c=[1.0,1.0],order::Int64=2)

    # Get h value
    h = hval(order)

    # Penalties
    α₀ = 0.5/(h * Δx) # Derivative penatly
    τ₁ = -0.5/(h * Δx) # Symmeteriser penalty
    τ₀ = -max(c[1]/2(h*Δx),c[end]/2(h*Δx))/(h*Δx) # Dirichlet penalty

    uₓₓ[1]  += τ₀ * (u[1] - u[end])
    uₓₓ[end]+= τ₀ * (u[end] - u[1])

    # Symmeteriser
    L₁u = (u[1] - u[end])

    uₓₓ[1:order]        .+= τ₁ * c[1] * BDₓᵀ(L₁u,Left,Δx,order)
    uₓₓ[end-order+1:end].+= τ₁ * c[end] * BDₓᵀ(L₁u,Right,Δx,order)

    # Neumann term
    BDₓu₁ = BDₓ(u,Left,Δx,order)
    BDₓuₙ = BDₓ(u,Right,Δx,order)

    uₓₓ[1]  += α₀ * (c[1]*BDₓu₁ - c[end]*BDₓuₙ)
    uₓₓ[end]+= α₀ * (c[1]*BDₓu₁ - c[end]*BDₓuₙ)
    # return uₓₓ

end



"""
    SAT_Periodic!


"""
function SAT_Periodic!(cache::AbstractArray,u::AbstractArray,c::AbstractArray,
        α₀::Real,τ₀::Function,τ₁::Real,
        E₀Dₓ::Vector{Real},EₙDₓ::Vector{Real},BDₓᵀ::Vector{Real},
        order::Int,loopaxis::Function)

    for (S,U,K) in zip(loopaxis(cache),loopaxis(u),loopaxis(c))
        # Dirichlet terms
        S[1]  += τ₀(K) * (U[1] - U[end])
        S[end]+= τ₀(K) * (U[end] - U[1])
        # Symmeteriser
        L₁u = (U[1] - U[end])
        S[1:order]        .+= τ₁ * K[1] * BDₓᵀ*L₁u
        S[end-order+1:end].+= τ₁ * K[end] * BDₓᵀ*L₁u
        # Neumann terms
        S[1]  += α₀ * K[1]*dot(E₀Dₓ,U[1:order]) - K[end]*dot(EₙDₓ,U[end-order+1:end])
        S[end]+= α₀ * K[1]*dot(E₀Dₓ,U[1:order]) - K[end]*dot(EₙDₓ,U[end-order+1:end])
    end
    cache
end
