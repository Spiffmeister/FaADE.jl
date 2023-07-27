
"""
    SAT_Periodic

Storage of all objects needed for a Periodic SAT ``u(x_0) = u(x_N)`` and ``\\left.\\partial_x u\\right|_{x_0} = \\left.\\partial_x u\\right|_{x_N}``.
"""
struct SAT_Periodic{
        TT<:Real,
        VT<:Vector{TT},
        F1<:Function,F2<:Function} <: SimultanousApproximationTerm{:Periodic}
    
    type        :: BoundaryConditionType
    axis        :: Int
    order       :: Int
    D₁ᵀE₀       :: VT
    D₁ᵀEₙ       :: VT
    E₀D₁        :: VT
    EₙD₁        :: VT
    Δx          :: TT
    α₀          :: TT
    τ₁          :: TT
    τ₀          :: F1
    loopaxis    :: F2

    function SAT_Periodic(Δx::TT,axis::Int,order::Int) where TT

        D₁ᵀE₀ = BoundaryDerivativeTranspose(Left,order,Δx)
        D₁ᵀEₙ = BoundaryDerivativeTranspose(Right,order,Δx)
        E₀D₁ = BoundaryDerivative(Left,Δx,order)
        EₙD₁ = BoundaryDerivative(Right,Δx,order)

        α₀, τ₁, τ₀ = SATpenalties(Periodic,Δx,order)

        loopaxis = SelectLoopDirection(axis)


        new{TT,Vector{TT},typeof(τ₀),typeof(loopaxis)}(
            Periodic,axis,order,D₁ᵀE₀,D₁ᵀEₙ,E₀D₁,EₙD₁,Δx,α₀,τ₁,τ₀)
    end
end




"""
    generate_Periodic
"""
function generate_Periodic(SATP::SAT_Periodic,solver)

    loopdirection = SelectLoopDirection(SATP.axis)

    let α₀ = SATP.α₀, τ₀ = SATP.τ₀, τ₁ = SATP.τ₁,
        E₀D₁ = SATP.E₀D₁, EₙD₁ = SATP.EₙD₁, D₁ᵀE₀ = SATP.D₁ᵀE₀, D₁ᵀEₙ = SATP.D₁ᵀEₙ,
        order=SATP.order

        Term!(cache::Array,u::Array,c::Array) = 
            SAT_Periodic!(cache,u,c,α₀,τ₀,τ₁,E₀D₁,EₙD₁,D₁ᵀE₀,D₁ᵀEₙ,order,loopdirection)
        return Term!
    end
end






"""
    SAT_Periodic!
"""
function SAT_Periodic!(cache::AT,u::AT,c::AT,
        α₀::T,τ₀::Function,τ₁::T,
        E₀D₁::Vector{T},EₙD₁::Vector{T},D₁ᵀE₀::Vector{T},D₁ᵀEₙ::Vector{T},
        order::Int,loopaxis::F) where {T,F,AT}

    for (S,U,K) in zip(loopaxis(cache),loopaxis(u),loopaxis(c))
        # Dirichlet terms
        S[1]  += τ₀(K) * (U[1] - U[end])
        S[end]+= τ₀(K) * (U[end] - U[1])
        # Symmeteriser
        L₁u = (U[1] - U[end])
        for i = 1:order
            S[i]            += τ₁*K[1]*D₁ᵀE₀[i]*L₁u
            S[end-order+i]  += τ₁*K[end]*D₁ᵀEₙ[i]*L₁u
            #Neumann terms
            S[1]  += α₀ * (K[1]*E₀D₁[i]*U[i] - K[end]*EₙD₁[i]*U[end-order+i])
            S[end]+= α₀ * (K[1]*E₀D₁[i]*U[i] - K[end]*EₙD₁[i]*U[end-order+i])
        end
        # S[1]  += α₀ * (K[1]*dot(E₀D₁,U[1:order]) - K[end]*dot(EₙD₁,U[end-order+1:end]))
        # S[end]+= α₀ * (K[1]*dot(E₀D₁,U[1:order]) - K[end]*dot(EₙD₁,U[end-order+1:end]))
        # S[1:order]        .+= τ₁ * K[1] * BD₁ᵀ*L₁u
        # S[end-order+1:end].+= -τ₁ * K[end] * BD₁ᵀ*L₁u
        # # Neumann terms
        # S[1]  += α₀ * (K[1]*dot(E₀D₁,U[1:order]) - K[end]*dot(EₙD₁,U[end-order+1:end]))
        # S[end]+= α₀ * (K[1]*dot(E₀D₁,U[1:order]) - K[end]*dot(EₙD₁,U[end-order+1:end]))
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
            S[i]            += SP.τ₁*K[1]*SP.D₁ᵀE₀[i]*L₁u
            S[end-SP.order+i]  += SP.τ₁*K[end]*SP.D₁ᵀEₙ[i]*L₁u
        end
        S[1]  += SP.α₀ * (K[1]*dot(SP.E₀D₁,U[1:SP.order]) - K[end]*dot(SP.EₙD₁,U[end-SP.order+1:end]))
        S[end]+= SP.α₀ * (K[1]*dot(SP.E₀D₁,U[1:SP.order]) - K[end]*dot(SP.EₙD₁,U[end-SP.order+1:end]))
    end
    cache
end
