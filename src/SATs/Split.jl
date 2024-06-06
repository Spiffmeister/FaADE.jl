



struct SAT_Interface{
        TN<:NodeType,
        COORD,
        TT<:Real,
        TV<:Vector{TT},
        F1<:Function,
        F2<:Function} <: SimultanousApproximationTerm{:Interface}

    side    :: TN
    axis    :: Int
    order   :: Int
    D₁ᵀE₀   :: TV
    D₁ᵀEₙ   :: TV
    D₁E₀    :: TV
    D₁Eₙ    :: TV
    τ₀      :: F1
    τ₁      :: TT
    τ₂      :: TT
    loopaxis :: F2

    Δy      :: TT
    coordinates :: Symbol

    function SAT_Interface(Δx₁::TT,Δx₂,side::TN,axis::Int,order::Int;Δy=TT(0),coordinates=:Cartesian) where {TT,TN}
        # Δxₗ = Δx₁
        # Δxᵣ = Δx₂
        if TN <: NodeType{:Left} #If the left boundary it should be this way around
            Δxₗ = Δx₁
            Δxᵣ = Δx₂
            
            # τ₀, τ₁, τ₂ = SATpenalties(Interface,Δxₗ,Δxᵣ,order)
        else
            Δxₗ = Δx₂
            Δxᵣ = Δx₁
            
            # τ₀, τ₁, τ₂ = SATpenalties(Interface,Δxᵣ,Δxₗ,order)
        end

        # τ₁ = -TT(1//2) / h
        # τ₂ = TT(1//2)
        # τ₀(c) = -TT(1//2) * c * 1/(h * min(Δx₁,Δx₂)^2)

        # τ₂ penalties
        D₁ᵀE₀ = BoundaryDerivativeTranspose(Left,order,Δxₗ^2) # H⁻¹D₁ᵀE₀
        D₁ᵀEₙ = BoundaryDerivativeTranspose(Right,order,Δxᵣ^2) # H⁻¹D₁ᵀEₙ
        # τ₁ penalties
        E₀D₁ = BoundaryDerivative(Left,Δxₗ,order)
        EₙD₁ = BoundaryDerivative(Right,Δxᵣ,order)

        τ₀, τ₁, τ₂ = SATpenalties(Interface,Δx₁,Δx₂,order)

        loopaxis = SelectLoopDirection(axis)

        new{TN,coordinates,TT,Vector{TT},typeof(τ₀),typeof(loopaxis)}(side,axis,order,
            D₁ᵀE₀,D₁ᵀEₙ,E₀D₁,EₙD₁,τ₀,τ₁,τ₂,loopaxis,Δy,coordinates)
    end

end




function (SI::SAT_Interface{NodeType{:Left,DIM},TT})(cache::AT,c::AT,u::AT,buffer::AT,::SATMode{:SolutionMode}) where {TT,DIM,AT}
    # rhs = u⁻ - u⁺
    # println("Left ",u," ",buffer)
    for (S⁺,U⁺,K⁺,U⁻) in zip(SI.loopaxis(cache),SI.loopaxis(u),SI.loopaxis(c),SI.loopaxis(buffer))
        S⁺[1] += SI.τ₀(K⁺[1]) * (U⁻[end] - U⁺[1])
        for i = 1:SI.order
            # τ₁ K D₁ᵀL₀ u + αL₀KD₁u
            S⁺[i] += SI.τ₁*K⁺[1]*SI.D₁ᵀE₀[i]*(U⁻[end] - U⁺[1])
            S⁺[1] += SI.α₀ * K⁺[1] * (SI.D₁E₀[i]*U⁻[end-SI.order+i] - SI.D₁Eₙ[i]*U⁺[i])
        end
        # S⁺[1] += -τ₀(K⁺[1]) * (U⁻[1] - S⁺[1]) # L₀u = u⁻ - u⁺
    end
end
function (SI::SAT_Interface{NodeType{:Right,DIM},TT})(cache::AT,c::AT,u::AT,buffer::AT,::SATMode{:SolutionMode}) where {TT,DIM,AT}
    # println("Right ",u," ",buffer)
    for (S⁻,U⁻,K⁻,U⁺) in zip(SI.loopaxis(cache),SI.loopaxis(u),SI.loopaxis(c),SI.loopaxis(buffer))
        S⁻[end] += SI.τ₀(K⁻[end]) * (U⁻[end] - U⁺[1])
        for i = 1:SI.order
            S⁻[end-SI.order+i] += SI.τ₁ * K⁻[end] * SI.D₁ᵀE₀[i] * (U⁻[end] - U⁺[1])
            S⁻[end] += SI.α₀ * K⁻[end] * (SI.D₁E₀[i]*U⁻[end-SI.order+i] - SI.D₁Eₙ[i]*U⁺[i])
        end
        # S⁻[end] += SI.τ₀ * U⁻[end]
    end
end

#######

"""
    SAT_Interface!
"""
function SAT_Interface! end
"""
    SAT_Interface!
Left handed SAT for interface conditions. Correspond to block 2 in the setup

|---|---|
| 1 | 2 |
|---|---|

Superscript + is the current block - is the joining block
"""
function SAT_Interface!(dest::AT,u::AT,c::AT,buffer::AT,SI::SAT_Interface{TN},::SATMode{:SolutionMode}) where {AT,TN<:Union{NodeType{:Left},NodeType{:Down}}}
    for (S⁺,U⁺,K⁺,U⁻) in zip(SI.loopaxis(dest),SI.loopaxis(u),SI.loopaxis(c),SI.loopaxis(buffer))
        S⁺[1] += SI.τ₀(K⁺[1]) * (U⁺[1] - U⁻[1])
        for i = 1:SI.order
            # τ₁ K D₁ᵀL₀ u + αL₀KD₁u
            S⁺[i] += SI.τ₂ * K⁺[1] * SI.D₁ᵀE₀[i]*(U⁻[1] - U⁺[1])
            S⁺[1] += SI.τ₁ * K⁺[1] * -SI.D₁E₀[i]*U⁺[i]
            # S⁺[1] += SI.τ₁ * K⁺[1] * (SI.D₁Eₙ[i]*U⁻[end-SI.order+i] - SI.D₁E₀[i]*U⁺[i])
        end
        S⁺[1] += SI.τ₁ * U⁻[2]
        # @show K⁺[1] * SI.D₁E₀[1]*U⁺[1] + K⁺[1] * SI.D₁E₀[2]*U⁺[2]
    end
    dest
end
"""
    SAT_Interface!
Right handed SAT for interface conditions. Correspond to block 1 in the setup

|---|---|
| 1 | 2 |
|---|---|

Superscript - is the current block + is the joining block
"""
function SAT_Interface!(dest::AT,u::AT,c::AT,buffer::AT,SI::SAT_Interface{TN},::SATMode{:SolutionMode}) where {AT,TN<:Union{NodeType{:Right},NodeType{:Up}}}
    for (S⁻,U⁻,K⁻,U⁺) in zip(SI.loopaxis(dest),SI.loopaxis(u),SI.loopaxis(c),SI.loopaxis(buffer))
        S⁻[end] += SI.τ₀(K⁻[end]) * (U⁻[end] - U⁺[1])
        # tmp = 0.0
        # tmp2 = 0.0
        for i = 1:SI.order
            S⁻[end-SI.order+i]  += SI.τ₂ * K⁻[end] * SI.D₁ᵀEₙ[i] * (U⁻[end] - U⁺[1])
            S⁻[end]             += SI.τ₁ * K⁻[end] * SI.D₁Eₙ[i]*U⁻[end-SI.order+i]
            # S⁻[end]             += SI.τ₁ * K⁻[end] * (SI.D₁Eₙ[i]*U⁻[end-SI.order+i] - SI.D₁E₀[i]*U⁺[i])
            # tmp = tmp + K⁻[end-SI.order+i] * SI.D₁Eₙ[i]*U⁻[end-SI.order+i]
            # tmp = tmp + K⁻[end] * SI.D₁Eₙ[i]*U⁻[end-SI.order+i]
        end
        S⁻[end] += -SI.τ₁ * U⁺[2]
        # @show K⁻[end-1] * SI.D₁Eₙ[1]*U⁻[end-1] + K⁻[end] * SI.D₁Eₙ[2]*U⁻[end]
        # @show "VALUE", tmp
    end
    dest
end


function SAT_Interface!(dest::AT,u::AT,cx::AT,cxy::AT,buffer::AT,SI::SAT_Interface{TN,:Curvilinear,TT},::SATMode{:SolutionMode}) where {AT,TN,TT}
    # @show TN, "SAT", cx[1], cxy[1]
    SAT_Interface!(dest,u,cx,buffer,SI,SolutionMode)
    n = size(dest,SI.axis)
    m = size(dest,mod1(SI.axis+1,2))

    if SI.side == Left
        DEST= view(dest,1, 1:m)
        SRC = view(u,  1,  1:m)
        C   = view(TT(-1)*SI.τ₁*cxy, 1,  1:m)

        # sgn = TT(1) * SI.τ₁
    elseif SI.side == Right
        DEST= view(dest,n,  1:m)
        SRC = view(u,   n,  1:m)
        C   = view(TT(1)*SI.τ₁*cxy, n,  1:m)

        # sgn = TT(-1) * SI.τ₁
    elseif SI.side == Up
        DEST = view(dest,   1:m, n)
        SRC = view(u,       1:m, n)
        C = view(cxy[:,1],  1:m, n) # needs to be fixed, should pull from both regions
    else
        DEST = view(dest,   1:m, n)
        SRC = view(u[:,1],  1:m, n)
        C = view(cxy[:,1],  1:m, n) # needs to be fixed, should pull from both regions
    end
    # @show TN
    # @show SRC
    # @show DEST
    D₁!(DEST,C,SRC,m,SI.Δy,Val(SI.order),TT(1))
    # FirstDerivativeTranspose!(DEST,buff,C,m,SI.Δy,SI.order,TT(1))
    # @show DEST
    dest
end




"""
    SAT_Interface_cache!
Computes the required values for sending to the buffer
"""
function SAT_Interface_cache! end
"""
    SAT_Interface_cache!
Computes the required values from the RIGHT handed block for sending to the buffer for LEFT handed interface conditions
"""
function SAT_Interface_cache!(dest::AT,u::AT,c::AT,SI::SAT_Interface{TN,COORD,TT}) where {TT,AT,COORD,TN<:Union{NodeType{:Left},NodeType{:Down}}}
    for (S⁻,U⁻,K⁻) in zip(SI.loopaxis(dest),SI.loopaxis(u),SI.loopaxis(c))
        S⁻[1] = U⁻[end]
        S⁻[2] = TT(0)
        for i = 1:SI.order
            S⁻[2] += K⁻[end] * SI.D₁Eₙ[i] * U⁻[end-SI.order+i]
        end
    end
    dest
end
"""
    SAT_Interface_cache!
Computes the required values from the LEFT handed block for sending to the buffer for RIGHT handed interface conditions
"""
function SAT_Interface_cache!(dest::AT,u::AT,c::AT,SI::SAT_Interface{TN,COORD,TT}) where {TT,AT,COORD,TN<:Union{NodeType{:Right},NodeType{:Up}}}
    for (S⁺,U⁺,K⁺) in zip(SI.loopaxis(dest),SI.loopaxis(u),SI.loopaxis(c))
        S⁺[1] = U⁺[1]
        S⁺[2] = TT(0)
        # tmp = TT(0)
        for i = 1:SI.order
            S⁺[2] += 1/(J[1] * sqrt(qx[1]^2 + qy[1]^2)) * K⁺[1] * SI.D₁E₀[i] * U⁺[i]
            # tmp = tmp + K⁺[1] * SI.D₁E₀[i] * U⁺[i]
        end
        # @show "CACHE", tmp
        # @show K⁺[1] * SI.D₁E₀[1] * U⁺[1] + K⁺[2] * SI.D₁E₀[2] * U⁺[2]
        # K_q[1,1] * D₁E₀[1] * U[1,1] + K_q[1,2] * D₁E₀[2] * U[1,2]
    end
    dest
end
"""
    SAT_Interface_cache!
Computes the required values from the LEFT handed block for sending to the buffer for RIGHT handed interface conditions
"""
function SAT_Interface_cache!(dest::AT,u::AT,c::AT,cxy::AT,SI::SAT_Interface{TN,:Curvilinear,TT}) where {TT,AT,TN}
    # @show TN, "CACHE", c[1], cxy[1]
    SAT_Interface_cache!(dest,u,c,SI)

    m = size(dest,mod1(SI.axis+1,2))

    if SI.side == Left
        @views DEST = dest[2,:]
        @views SRC = u[end,:]
        @views C = cxy[end,:]
    elseif SI.side == Right
        @views DEST= dest[2,:]
        @views SRC = u[1,:] 
        @views C   = cxy[1,:]
    elseif SI.side == Up
    else
    end

    # @show TN, "cache"
    # @show DEST
    # @show SRC
    # @show c[1,:]
    D₁!(DEST,C,SRC,m,SI.Δy,Val(SI.order),TT(1))
    # @show DEST
    # @show ""

    # FirstDerivativeTranspose!(DEST,SRC,C,m,SI.Δy,SI.order,TT(1))

    dest
end








