
"""
    SAT_Periodic

Storage of all objects needed for a Periodic SAT ``u(x_0) = u(x_N)`` and ``\\left.\\partial_x u\\right|_{x_0} = \\left.\\partial_x u\\right|_{x_N}``.
"""
struct SAT_Periodic{
        TN<:NodeType,
        TT<:Real,
        VT<:Vector{TT},
        F1<:Function,F2<:Function} <: SimultanousApproximationTerm{:Periodic}
    type        :: BoundaryConditionType
    side        :: TN
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

    Δy          :: TT
    coordinates :: Symbol

    
end
function SAT_Periodic(Δx::TT,axis::Int,order::Int) where TT

    D₁ᵀE₀ = BoundaryDerivativeTranspose(Left,order,Δx)
    D₁ᵀEₙ = BoundaryDerivativeTranspose(Right,order,Δx)
    E₀D₁ = BoundaryDerivative(Left,Δx,order)
    EₙD₁ = BoundaryDerivative(Right,Δx,order)

    α₀, τ₁, τ₀ = SATpenalties(Periodic,Δx,order)

    loopaxis = SelectLoopDirection(axis)

    SAT_Periodic{typeof(Internal),TT,Vector{TT},typeof(τ₀),typeof(loopaxis)}(
        Periodic,Internal,axis,order,D₁ᵀE₀,D₁ᵀEₙ,E₀D₁,EₙD₁,Δx,α₀,τ₁,loopaxis,0.0,:Cartesian)
end
function SAT_Periodic(Δx::TT,axis::Int,order::Int,side::NodeType) where TT
    D₁ᵀE₀ = BoundaryDerivativeTranspose(Left,order,Δx)
    D₁ᵀEₙ = BoundaryDerivativeTranspose(Right,order,Δx)
    E₀D₁ = BoundaryDerivative(Left,Δx,order)
    EₙD₁ = BoundaryDerivative(Right,Δx,order)

    α₀, τ₁, τ₀ = SATpenalties(Periodic,Δx,order)


    # E₀ = _BoundaryOperator(TT,Left)
    # Eₙ = _BoundaryOperator(TT,Right)

    # HinvL = _InverseMassMatrix(order,Δx,side)
    # HinvR = _InverseMassMatrix(order,Δx,side)

    # D₁ᵀ = _DerivativeTranspose(order,Δx,side)
    # D₁ = _BoundaryDerivative(order,Δx,side)


    # H⁻¹L₀ = HinvL[1].*[1,-1]

    # H⁻¹ = HinvL
    # H⁻¹ = HinvL

    # D₁ᵀE₀ = _BoundaryDerivative(order,Δx,Left)
    # D₁ᵀEₙ = _BoundaryDerivative(order,Δx,Right)


    # α₀ = 1/2
    # τ₁ = 1/2
    # τ₀ = c -> -(1+1/c) * c/(2HinvL[1])



    loopaxis = SelectLoopDirection(axis)

    SAT_Periodic{typeof(side),TT,Vector{TT},typeof(τ₀),typeof(loopaxis)}(
        Periodic,side,axis,order,D₁ᵀE₀,D₁ᵀEₙ,E₀D₁,EₙD₁,Δx,α₀,τ₁,τ₀,loopaxis,0.0,:Cartesian)
end
function SAT_Periodic(Δx::TT,axis::Int,order::Int,side::NodeType,Δy::TT,coord::Symbol) where TT
    D₁ᵀE₀ = BoundaryDerivativeTranspose(Left,order,Δx)
    D₁ᵀEₙ = BoundaryDerivativeTranspose(Right,order,Δx)
    E₀D₁ = BoundaryDerivative(Left,Δx,order)
    EₙD₁ = BoundaryDerivative(Right,Δx,order)

    α₀, τ₁, τ₀ = SATpenalties(Periodic,Δx,order)

    loopaxis = SelectLoopDirection(axis)

    SAT_Periodic{typeof(side),TT,Vector{TT},typeof(τ₀),typeof(loopaxis)}(
        Periodic,side,axis,order,D₁ᵀE₀,D₁ᵀEₙ,E₀D₁,EₙD₁,Δx,α₀,τ₁,τ₀,loopaxis,Δy,coord)
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



"""
    SAT_Periodic!(cache::AT,u::AT,c::AT,SP::SAT_Periodic{NodeType{:Internal}}) where {AT}
Periodic SAT for a single block
"""
function SAT_Periodic!(cache::AT,u::AT,c::AT,SP::SAT_Periodic{NodeType{:Internal}}) where {AT}
    for (S,U,K) in zip(SP.loopaxis(cache),SP.loopaxis(u),SP.loopaxis(c))
        # Dirichlet terms
        S[1]  += SP.τ₀(K) * (U[1] - U[end])
        S[end]+= SP.τ₀(K) * (U[end] - U[1])
        # Symmeteriser
        L₁u = (U[1] - U[end])
        for i = 1:order
            S[i]            += SP.τ₁*SP.K[1]*SP.D₁ᵀE₀[i]*L₁u
            S[end-order+i]  += SP.τ₁*SP.K[end]*SP.D₁ᵀEₙ[i]*L₁u
            #Neumann terms
            S[1]  += SP.α₀ * (K[1]*SP.E₀D₁[i]*U[i] - K[end]*SP.EₙD₁[i]*U[end-order+i])
            S[end]+= SP.α₀ * (K[1]*SP.E₀D₁[i]*U[i] - K[end]*SP.EₙD₁[i]*U[end-order+i])
        end
    end
    cache
end


"""
    SAT_Periodic!(cache::AT,u::AT,c::AT,SP::SAT_Periodic{NodeType{:Left}},::SATMode{:SolutionMode}) where {AT}
Periodic SAT for left boundary of a split domain
"""
function SAT_Periodic!(dest::VT,u::VT,c::VT,SP::SAT_Periodic{TN}) where {VT<:AbstractVector,TN<:NodeType{:Left}}
    # for (S,U,K) in zip(SP.loopaxis(dest),SP.loopaxis(u),SP.loopaxis(c))
        # Dirichlet terms
        dest[1]  += SP.τ₀(c) * (u[1] - u[end])
        # Symmeteriser
        L₁u = (u[1] - u[end])
        for i = 1:SP.order
            dest[i]    += SP.τ₁*c[1]*SP.D₁ᵀE₀[i]*L₁u
            #Neumann terms
            dest[1]    += SP.α₀ * (c[1]*SP.E₀D₁[i]*u[i] - c[end]*SP.EₙD₁[i]*u[end-SP.order+i])
        end
        # println(dest[1])
    # end
    dest
end
"""
    SAT_Periodic!(cache::AT,u::AT,c::AT,SP::SAT_Periodic{NodeType{:Right}},::SATMode{:SolutionMode}) where {AT}
Periodic SAT for right boundary of a split domain
"""
function SAT_Periodic!(dest::VT,u::VT,c::VT,SP::SAT_Periodic{TN}) where {VT<:AbstractVector,TN<:NodeType{:Right}}
    # for (S,U,K) in zip(SP.loopaxis(dest),SP.loopaxis(u),SP.loopaxis(c))
        # Dirichlet terms
        dest[end]+= SP.τ₀(c) * (u[end] - u[1])
        # Symmeteriser
        L₁u = (u[1] - u[end])
        for i = 1:SP.order
            dest[end-SP.order+i]  += SP.τ₁*c[end]*SP.D₁ᵀEₙ[i]*L₁u
            #Neumann terms
            dest[end]   += SP.α₀ * (c[1]*SP.E₀D₁[i]*u[i] - c[end]*SP.EₙD₁[i]*u[end-SP.order+i])
        end
    # end
    dest
end
function SAT_Periodic!(dest::AT,u::AT,c::AT,SP::SAT_Periodic{TN,TT}) where {TT,AT<:AbstractMatrix{TT},TN<:NodeType}
    for (DEST,U,C) in zip(SP.loopaxis(dest),SP.loopaxis(u),SP.loopaxis(c))
        SAT_Periodic!(DEST,U,C,SP)
    end

    if SP.coordinates == :Curvilinear
        n = size(dest,SP.axis)
        m = size(dest,mod1(SP.axis+1,2))

        # @show n, m
        # @show SP.side, SP.axis
        # @show size(dest), size(u)

        if SP.side == Left
            DEST = view(dest,   1,1:m)
            SRC = view(u[1,1:m]-u[m,1:m], :,:)
        elseif SP.side == Right
            DEST = view(dest,   n,1:m)
            SRC = view(u,       n,1:m)
        elseif SP.side == Up
            DEST = view(dest,   1:m,1)
            SRC = view(u[:,1]-u[:,n], :,1)
        else
            DEST = view(dest,   1:m,n)
            SRC = view(u[:,1]-u[:,n],:,1)
        end
    # @show size(DEST),size(SRC)
        FirstDerivativeTranspose!(DEST,SRC,m,SP.Δy,SP.order,TT(1))
    end

    dest
end


"""
    SAT_Periodic!(cache::AT,u::AT,c::AT,SP::SAT_Periodic{NodeType{:Left}},::SATMode{:DataMode}) where {AT}
Periodic SAT for the left hand buffer
"""
function SAT_Periodic!(dest::AT,u::AT,c::AT,SP::SAT_Periodic{NodeType{:Left}},::SATMode{:DataMode}) where {AT}
    for (S,U,K) in zip(SP.loopaxis(dest),SP.loopaxis(u),SP.loopaxis(c))
        S[end] = U[end]
        S[end-1] = TT(0)
        for i = 1:SP.order
            S[end-1] += K[end]*SP.EₙD₁[i]*U[end-SP.order+i]
        end
    end
    dest
end
"""
    SAT_Periodic!(cache::AT,u::AT,c::AT,SP::SAT_Periodic{NodeType{:Right}},::SATMode{:DataMode}) where {AT}
Periodic SAT for the right hand buffer
"""
function SAT_Periodic!(dest::AT,u::AT,c::AT,SP::SAT_Periodic{NodeType{:Right}},::SATMode{:DataMode}) where {AT}
    for (S,U,K) in zip(SP.loopaxis(dest),SP.loopaxis(u),SP.loopaxis(c))
        S[1] = U[1]
        S[2] = TT(0)
        for i = 1:SP.order
            S[2] += K[1]*SP.E₀D₁[i]*U[i]
        end
    end
    dest
end









