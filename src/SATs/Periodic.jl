
struct SAT_Periodic{T} <: SimultanousApproximationTerm
    BDₓᵀ        :: Vector{T}
    E₀Dₓ        :: Vector{T}
    EₙDₓ        :: Vector{T}
    type        :: BoundaryConditionType
    axis        :: Int
    Δx          :: T
    α₀          :: T
    τ₁          :: T
    τ₀          :: T

    function Boundary_Periodic(Δx::T,axis::Int,order::Int) where T

        BDₓᵀ = BoundaryDerivativeTranspose(order,Δx)
        E₀Dₓ = BoundaryDerivative(Left,order,Δx)
        EₙDₓ = BoundaryDerivative(Right,order,Δx)

        α₀, τ₁, τ₀ = SATpenalties(Periodic,Δx,order)
        penalties = (α₀=α₀, τ₁=τ₁, τ₀=τ₀)

        new{T}(BDₓᵀ,E₀Dₓ,EₙDₓ,Periodic,axis,Δx,penalties)
    end
end




"""
    generate_Periodic
"""
function generate_Periodic(SATP::SAT_Periodic,solver)

    loopdirection = SelectLoopDirection(SATP.axis)

    let α₀ = SATP.α₀, τ₀ = SATP.τ₀, τ₁ = SATP.τ₁,
        E₀Dₓ = SATP.E₀Dₓ, EₙDₓ = SATP.EₙDₓ, BDₓᵀ = SATP.BDₓᵀ, 
        order=SATP.order

        Term(cache,u,c) = 
            SAT_Periodic!(cache,u,c,α₀,τ₀,τ₁,E₀Dₓ,EₙDₓ,BDₓᵀ,order,loopdirection)
    end
    
    return Term
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
