
"""
    SAT_Periodic

Storage of all objects needed for a Periodic SAT ``u(x_0) = u(x_N)`` and ``\\left.\\partial_x u\\right|_{x_0} = \\left.\\partial_x u\\right|_{x_N}``.
"""
struct SAT_Periodic{
        TT<:Real,
        VT<:Vector{TT},
        F1<:Function,F2<:Function} <: SimultanousApproximationTerm
    
    type        :: BoundaryConditionType
    axis        :: Int
    order       :: Int
    DₓᵀE₀       :: VT
    DₓᵀEₙ       :: VT
    E₀Dₓ        :: VT
    EₙDₓ        :: VT
    Δx          :: TT
    α₀          :: TT
    τ₁          :: TT
    τ₀          :: F1
    loopaxis    :: F2

    function SAT_Periodic(Δx::TT,axis::Int,order::Int) where TT

        DₓᵀE₀ = BoundaryDerivativeTranspose(Left,order,Δx)
        DₓᵀEₙ = BoundaryDerivativeTranspose(Right,order,Δx)
        E₀Dₓ = BoundaryDerivative(Left,Δx,order)
        EₙDₓ = BoundaryDerivative(Right,Δx,order)

        α₀, τ₁, τ₀ = SATpenalties(Periodic,Δx,order)

        loopaxis = SelectLoopDirection(axis)


        new{TT,Vector{TT},typeof(τ₀),typeof(loopaxis)}(
            Periodic,axis,order,DₓᵀE₀,DₓᵀEₙ,E₀Dₓ,EₙDₓ,Δx,α₀,τ₁,τ₀)
    end
end




"""
    generate_Periodic
"""
function generate_Periodic(SATP::SAT_Periodic,solver)

    loopdirection = SelectLoopDirection(SATP.axis)

    let α₀ = SATP.α₀, τ₀ = SATP.τ₀, τ₁ = SATP.τ₁,
        E₀Dₓ = SATP.E₀Dₓ, EₙDₓ = SATP.EₙDₓ, DₓᵀE₀ = SATP.DₓᵀE₀, DₓᵀEₙ = SATP.DₓᵀEₙ,
        order=SATP.order

        Term!(cache::Array,u::Array,c::Array) = 
            SAT_Periodic!(cache,u,c,α₀,τ₀,τ₁,E₀Dₓ,EₙDₓ,DₓᵀE₀,DₓᵀEₙ,order,loopdirection)
        return Term!
    end
end






"""
    SAT_Periodic!
"""
function SAT_Periodic!(cache::AT,u::AT,c::AT,
        α₀::T,τ₀::Function,τ₁::T,
        E₀Dₓ::Vector{T},EₙDₓ::Vector{T},DₓᵀE₀::Vector{T},DₓᵀEₙ::Vector{T},
        order::Int,loopaxis::F) where {T,F,AT}

    for (S,U,K) in zip(loopaxis(cache),loopaxis(u),loopaxis(c))
        # Dirichlet terms
        S[1]  += τ₀(K) * (U[1] - U[end])
        S[end]+= τ₀(K) * (U[end] - U[1])
        # Symmeteriser
        L₁u = (U[1] - U[end])
        for i = 1:order
            S[i]            += τ₁*K[1]*DₓᵀE₀[i]*L₁u
            S[end-order+i]  += τ₁*K[end]*DₓᵀEₙ[i]*L₁u
            #Neumann terms
            # S[1]  += α₀ * (K[1]*E₀Dₓ[i]*U[i] - K[end]*EₙDₓ[i]*U[end-order+i])
            # S[end]+= α₀ * (K[1]*E₀Dₓ[i]*U[i] - K[end]*EₙDₓ[i]*U[end-order+i])
        end
        S[1]  += α₀ * (K[1]*dot(E₀Dₓ,U[1:order]) - K[end]*dot(EₙDₓ,U[end-order+1:end]))
        S[end]+= α₀ * (K[1]*dot(E₀Dₓ,U[1:order]) - K[end]*dot(EₙDₓ,U[end-order+1:end]))
        # S[1:order]        .+= τ₁ * K[1] * BDₓᵀ*L₁u
        # S[end-order+1:end].+= -τ₁ * K[end] * BDₓᵀ*L₁u
        # # Neumann terms
        # S[1]  += α₀ * (K[1]*dot(E₀Dₓ,U[1:order]) - K[end]*dot(EₙDₓ,U[end-order+1:end]))
        # S[end]+= α₀ * (K[1]*dot(E₀Dₓ,U[1:order]) - K[end]*dot(EₙDₓ,U[end-order+1:end]))
    end
    cache
end



function (SP::SAT_Periodic{T})(cache::AT,u::AT,c::AT) where {T,AT}
    
    # loopaxis = SelectLoopDirection(SP.axis)
    # stable_eachslice(x) = eachslice(x,dims=SP.axis) :: 
    
    for (S,U,K) in zip(SP.loopaxis(cache),SP.loopaxis(u),SP.loopaxis(c))
        # Dirichlet terms
        S[1]  += SP.τ₀(K) * (U[1] - U[end])
        S[end]+= SP.τ₀(K) * (U[end] - U[1])
        # Symmeteriser
        L₁u = (U[1] - U[end]) :: T
        for i = 1:SP.order
            S[i]            += SP.τ₁*K[1]*SP.DₓᵀE₀[i]*L₁u
            S[end-SP.order+i]  += SP.τ₁*K[end]*SP.DₓᵀEₙ[i]*L₁u
        end
        S[1]  += SP.α₀ * (K[1]*dot(SP.E₀Dₓ,U[1:SP.order]) - K[end]*dot(SP.EₙDₓ,U[end-SP.order+1:end]))
        S[end]+= SP.α₀ * (K[1]*dot(SP.E₀Dₓ,U[1:SP.order]) - K[end]*dot(SP.EₙDₓ,U[end-SP.order+1:end]))
    end
    cache
end
