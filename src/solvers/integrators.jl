
"""
    (RK::ExplicitBlock{T,N,4})
RK4 integrator
"""
function (RK::ExplicitBlock{T,N,AT,4})(RHS::Function,DBlock::DataBlock,t::Float64) where {T,N,AT}
    
    DBlock.uₙ₊₁ .= DBlock.u

    RHS(RK.k1,DBlock.u,         DBlock.K,t)
    RHS(RK.k2,DBlock.u+0.5RK.k1,DBlock.K,t+0.5RK.Δt)
    RHS(RK.k3,DBlock.u+0.5RK.k2,DBlock.K,t+0.5RK.Δt)
    RHS(RK.k4,DBlock.u+RK.k3,   DBlock.K,t+RK.Δt)
    DBlock.uₙ₊₁ .+= DBlock.u + RK.Δt/6.0 * (RK.k1 + 2RK.k2 + 2RK.k3 + RK.k4)

end
"""
    (ForwardEuler::ExplicitBlock{T,N,1})
Forward Euler integrator
"""
function (ForwardEuler::ExplicitBlock{T,N,AT,1})(RHS::Function,DBlock::DataBlock,t::Float64) where {T,N,AT}
    
    DBlock.uₙ₊₁ .= DBlock.u

    RHS(RK.k1,DBlock.u,DBlock.K,t)
    DBlock.uₙ₊₁ .+= DBlock.u + RK.Δt*RK.k1

end



#= IMPLICIT METHODS =#
"""
    implicit_euler
TODO: REMOVE
"""
function implicit_euler end
# Vector form
function implicit_euler(uₙ::Vector,uₒ::Vector,RHS::Function,n::Int,Δx::Float64,Δt::Float64,k::Vector,t::Float64,x::Vector,boundary;maxIT::Int=100,α::Float64=1.5)
    error("Do not use")
    uⱼ = uₒ
    for j = 1:maxIT
        uⱼ = uⱼ - α * Δt * (uⱼ - uₒ - Δt*RHS(uₙ,uⱼ,n,x,Δx,t,Δt,k,boundary))
    end
    uⱼ = uₒ + Δt*RHS(uₙ,uⱼ,n,x,Δx,t,Δt,k,boundary)
    return uⱼ
end
# Matrix form
function implicit_euler(uₙ::Matrix,uₒ::Matrix,RHS::Function,nx::Int,ny::Int,Δx::Float64,Δy::Float64,Δt::Float64,kx::Matrix,ky::Matrix,t::Float64,x::Vector,y::Vector,boundary_x,boundary_y;maxIT::Int=100,α::Float64=1.5)
    error("Do not use")
    uⱼ = uₒ
    for j = 1:maxIT
        uⱼ = uⱼ - α * Δt * (uⱼ - uₒ - Δt*RHS(uₙ,uⱼ,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,boundary_x,boundary_y))
    end
    uⱼ = uₒ + Δt*RHS(uₙ,uⱼ,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,boundary_x,boundary_y)
    return uⱼ
end



"""
    conj_grad!
In-place conjugate gradient method.

See also [`build_H`](@ref), [`A!`](@ref), [`innerH`](@ref)
"""
function conj_grad!(RHS,DBlock::DataBlock{TT,DIM,AT},CGB::ConjGradBlock{TT,DIM,AT},Δt::TT;
        atol=1.e-5,rtol=1.e-10,maxIT::Int=10,warnings=false) where {TT,DIM,AT}
    
    local rnorm::TT
    local unorm::TT
    local dₖAdₖ::TT
    local βₖ::TT
    local αₖ::TT
    
    
    # x₀ = uₙ #Initial guess
    # CGB.b .= DBlock.uₙ₊₁ #uₙ₊₁ is our initial guess and RHS
    # rₖ = (uₙ₊₁ - Δt*uₓₓⁿ⁺¹) - (uₙ + F)
    A!(RHS,CGB.rₖ,DBlock.u,Δt,DBlock.K)
    CGB.rₖ .= CGB.rₖ .- CGB.b
    # DBlock.uₙ₊₁ .= DBlock.u

    CGB.dₖ .= -CGB.rₖ # dₖ = -rₖ
    
    i = 0
    rnorm = sqrt(CGB.innerprod(CGB.rₖ,CGB.rₖ))
    unorm = max(sqrt(CGB.innerprod(DBlock.uₙ₊₁,DBlock.uₙ₊₁)),1e-14)
    while (rnorm > rtol*unorm) & (i < maxIT)
        A!(RHS,CGB.Adₖ,CGB.dₖ,Δt,DBlock.K) # Adₖ = dₖ - Δt*D(dₖ)
        dₖAdₖ = CGB.innerprod(CGB.dₖ,CGB.Adₖ)
        αₖ = -CGB.innerprod(CGB.rₖ,CGB.dₖ)/dₖAdₖ
        @. DBlock.uₙ₊₁ = DBlock.uₙ₊₁ + αₖ*CGB.dₖ #xₖ₊₁ = xₖ + αₖ*dₖ
        
        # rₖ = (xₖ₊₁ - Δt*Dxₖ₊₁ - b)
        A!(RHS,CGB.rₖ,DBlock.uₙ₊₁,Δt,DBlock.K)
        @. CGB.rₖ = CGB.rₖ - CGB.b
        
        A!(RHS,CGB.Drₖ,CGB.rₖ,Δt,DBlock.K) # Drₖ = rₖ - Δt*D(rₖ)
        
        βₖ = CGB.innerprod(CGB.rₖ,CGB.Drₖ)/dₖAdₖ
        @. CGB.dₖ = βₖ*CGB.dₖ - CGB.rₖ

        rnorm = sqrt(CGB.innerprod(CGB.rₖ,CGB.rₖ))
        i += 1
    end
    if (rnorm>rtol*unorm) & warnings
        CGB.scalar.converged = false
        warnstr = string("CG did not converge with Δt=",Δt,", rel error=",rnorm/unorm,", rel tolerance=",rtol,".")
        @warn warnstr
    end
end
function conj_grad!(DBlock::newLocalDataBlock{TT,DIM,AT},Δt::TT;
    atol=1.e-5,rtol=1.e-10,maxIT::Int=10,warnings=false) where {TT,DIM,AT}

    local rnorm::TT
    local unorm::TT
    local dₖAdₖ::TT
    local βₖ::TT
    local αₖ::TT

    # x₀ = uₙ #Initial guess
    # CGB.b .= DBlock.uₙ₊₁ #uₙ₊₁ is our initial guess and RHS
    # rₖ = (uₙ₊₁ - Δt*uₓₓⁿ⁺¹) - (uₙ + F)
    A!(DBlock.rₖ,DBlock.u,DBlock,Δt)
    @. DBlock.rₖ = DBlock.rₖ - DBlock.b

    @. DBlock.dₖ = -DBlock.rₖ # dₖ = -rₖ

    i = 0
    rnorm = sqrt(DBlock.innerprod(DBlock.rₖ,DBlock.rₖ))
    unorm = max(sqrt(DBlock.innerprod(DBlock.uₙ₊₁,DBlock.uₙ₊₁)),1e-14)
    while (rnorm > rtol*unorm) & (i < maxIT)
        A!(DBlock.cache,DBlock.dₖ,DBlock,Δt) # Adₖ = dₖ - Δt*D(dₖ)
        dₖAdₖ = CGB.innerprod(DBlock.dₖ,DBlock.cache)
        αₖ = -DBlock.innerprod(DBlock.rₖ,DBlock.dₖ)/dₖAdₖ
        @. DBlock.uₙ₊₁ = DBlock.uₙ₊₁ + αₖ*DBlock.dₖ #xₖ₊₁ = xₖ + αₖ*dₖ

        # rₖ = (xₖ₊₁ - Δt*Dxₖ₊₁ - b)
        A!(DBlock.rₖ,DBlock.uₙ₊₁,DBlock,Δt)
        @. DBlock.rₖ = DBlock.rₖ - DBlock.b

        A!(DBlock.cache,DBlock.rₖ,DBlock,Δt) # Drₖ = rₖ - Δt*D(rₖ)

        βₖ = DBlock.innerprod(DBlock.rₖ,DBlock.cache)/dₖAdₖ
        @. DBlock.dₖ = βₖ*DBlock.dₖ - DBlock.rₖ

        rnorm = sqrt(DBlock.innerprod(DBlock.rₖ,DBlock.rₖ))
        i += 1
    end
    if (rnorm>rtol*unorm) & warnings
        CGB.scalar.converged = false
        warnstr = string("CG did not converge with Δt=",Δt,", rel error=",rnorm/unorm,", rel tolerance=",rtol,".")
        @warn warnstr
    end
end
function conj_grad!(DBlock::DataMultiBlock{TT,DIM},Δt::TT;
    atol=1.e-5,rtol=1.e-10,maxIT::Int=10,warnings=false) where {TT,DIM}

    local rnorm::TT
    local unorm::TT
    local dₖAdₖ::TT
    local βₖ::TT
    local αₖ::TT

    # x₀ = uₙ #Initial guess
    # CGB.b .= DBlock.uₙ₊₁ #uₙ₊₁ is our initial guess and RHS
    # rₖ = (uₙ₊₁ - Δt*uₓₓⁿ⁺¹) - (uₙ + F)
    # A!(CGB.rₖ,DBlock.u,DBlock)
    A!(:u,DBlock)
    setValue(:rₖ,:cache,DBlock)
    muladd!(:rₖ,:b,DBlock,β=TT(-1)) #rₖ = rₖ - b

    setValue(:dₖ,:rₖ,DBlock,TT(-1)) #dₖ = -rₖ

    i = 0
    rnorm = sqrt(innerprod(:rₖ,:rₖ,DBlock)) #√(rₖ,rₖ)ₕ
    unorm = max(sqrt(innerprod(:rₖ,:rₖ,DBlock)),1e-14) #√(rₖ,rₖ)ₕ
    # println("un+1 ",DBlock[1].uₙ₊₁)
    while (rnorm > rtol*unorm) & (i < maxIT)
        # Comm dₖ boundaries
        A!(:dₖ,DBlock) # Adₖ = dₖ - Δt*D(dₖ)
        dₖAdₖ = innerprod(:dₖ,:cache,DBlock) #dₖAdₖ = (dₖ,Ddₖ)
        αₖ = -innerprod(:rₖ,:dₖ,DBlock)/dₖAdₖ #αₖ = -√(rₖ,dₖ)/dₖAdₖ
        muladd!(:uₙ₊₁,:dₖ,DBlock,β=αₖ) #uₙ₊₁ = uₙ + αₖdₖ
        # rₖ = (xₖ₊₁ - Δt*Dxₖ₊₁ - b)
        # Comm uₙ₊₁ boundaries
        A!(:uₙ₊₁,DBlock)
        # @. CGB.rₖ = CGB.cache - CGB.b
        setValue(:rₖ,:cache,DBlock) # rₖ = cache
        muladd!(:rₖ,:b,DBlock,β=TT(-1))  # rₖ = rₖ - b

        # Comm rₖ boundaries
        A!(:rₖ,DBlock) # Drₖ = rₖ - Δt*D(rₖ)
        βₖ = innerprod(:rₖ,:cache,DBlock)/dₖAdₖ

        # @. CGB.dₖ = βₖ*CGB.dₖ - CGB.rₖ
        muladd!(:dₖ,:rₖ,DBlock,α=βₖ,β=TT(-1))

        # rnorm = sqrt(CGB.innerprod(CGB.rₖ,CGB.rₖ))
        rnorm = sqrt(innerprod(:rₖ,:rₖ,DBlock))
        i += 1
    end
    if (rnorm>rtol*unorm) & warnings
        DBlock.SC.converged = false
        warnstr = string("CG did not converge with Δt=",Δt,", rel error=",rnorm/unorm,", rel tolerance=",rtol,".")
        @warn warnstr
    end
end

"""
    A!
Mutable ``u - ΔtD(u)``
"""
function A!(PDE!::Function,tmp::AT,uⱼ::AT,Δt::TT,k::KT) where {TT,AT,KT}
    PDE!(tmp,uⱼ,k)
    @. tmp = uⱼ - Δt*tmp
end

@inline function A!(write::Symbol,read::Symbol,D::newLocalDataBlock{TT},Δt::TT) where {AT,TT}
    # Compute the derivatives
    W = getproperty(D,:cache)
    R = getproperty(D,write)
    D.Derivative(W,R,D.K)
    # Communicate the boundaries
    # Apply needed SATs
    applySATs(W,R,D,SolutionMode)
    # u = r - Δt Du
    @. W = R - Δt*W
end

function A!(source::Symbol,DB::DataMultiBlock{TT,DIM}) where {TT,DIM}
    fillBuffers(source,DB)
    for i in eachblock(DB)
        # cache = u - Δt Du
        C = getproperty(DB[i],:cache) :: AbstractArray{TT,DIM}
        U = getproperty(DB[i],source) :: AbstractArray{TT,DIM}
        # Compute the derivatives
        DB[i].Derivative(C,U,DB[i].K)
        applySATs(C,U,DB[i],SolutionMode) #Apply the SATs
        muladd!(C,U,-DB[i].SC.Δt,TT(1))
        # @. DB[i].cache = U - Δt*DB[i].cache

        # A!(DB[i].cache,U,DB[i],Δt)
    end
end


