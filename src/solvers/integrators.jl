
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



"""
    conj_grad!
In-place conjugate gradient method.

See also [`build_H`](@ref), [`A!`](@ref), [`innerH`](@ref)
"""
function conj_grad!(RHS,DBlock::DataBlock{TT,DIM,AT},CGB::ConjGradBlock{TT,DIM,AT},Δt::TT;
        atol=1.e-5,rtol=1.e-14,maxIT::Int=10,warnings=false) where {TT,DIM,AT}
    # LEGACY
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
function conj_grad!(DBlock::DataMultiBlock{TT,DIM};
    atol=1.e-5,rtol=1.e-14,maxIT::Int=10,warnings=false) where {TT,DIM}

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

    # mul!(DBlock[1].cache,DBlock[1].u,DBlock[1].K,DBlock[1].Derivative)
    # A!(:u,DBlock)
    # DBlock[1].cache .= DBlock.SC.Δt*DBlock[1].cache
    # addSource!(DBlock[1].source,DBlock[1].cache,DBlock[1].grid,DBlock[1].SC.t,DBlock[1].SC.Δt,0.0)
    # rnorm = sqrt(innerprod(:cache,:cache,DBlock)) #√(rₖ,rₖ)ₕ

    i = 0
    rnorm = sqrt(innerprod(:rₖ,:rₖ,DBlock)) #√(rₖ,rₖ)ₕ
    unorm = max(sqrt(innerprod(:uₙ₊₁,:uₙ₊₁,DBlock)),1e-14) #√(uₙ₊₁,uₙ₊₁)ₕ
    # println("norms ",rnorm," ",unorm)
    while (rnorm > rtol*unorm) & (i < maxIT)
        # Comm dₖ boundaries
        A!(:dₖ,DBlock) # Adₖ = dₖ - Δt*D(dₖ)
        dₖAdₖ = innerprod(:dₖ,:cache,DBlock) #dₖAdₖ = (dₖ,Ddₖ)
        αₖ = -innerprod(:rₖ,:dₖ,DBlock)/dₖAdₖ #αₖ = -(rₖ,dₖ)/dₖAdₖ
        muladd!(:uₙ₊₁,:dₖ,DBlock,β=αₖ) #uₙ₊₁ = uₙ + αₖdₖ
        # rₖ = (xₖ₊₁ - Δt*Dxₖ₊₁ - b)
        A!(:uₙ₊₁,DBlock)
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
    # println("unp1 ",DBlock[1].uₙ₊₁[100])
    if (rnorm>rtol*unorm) & warnings
        DBlock.SC.converged = false
        warnstr = string("CG did not converge with Δt=",DBlock.SC.Δt,", rel error=",rnorm/unorm,", rel tolerance=",rtol,".")
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
    tmp
end
# b = (1+(1-θ)ΔtD)u^n + (1-θ)ΔtF^n + θΔtF^{n+1} + θSAT_data^{n+1} + (1-θ)SAT_data^{n}
# rk = (1-θΔtD)u^n
#b = (1+ΔtD)u^n + ΔtF^n + SAT_data^{n}
#rk = u^n
#b - rk = ΔtD(u^n) + ΔtF^n
function A!(read::Symbol,D::newLocalDataBlock{TT,DIM,AT,KT,DCT,GT,BT,DT,PT}) where {TT,DIM,AT,KT,DCT,GT,BT,DT,PT}
    # Compute the derivatives
    W = getarray(D,:cache)
    R = getarray(D,read)
    # D.Derivative(W,R,D.K)x
    mul!(W,R,D.K,D.Derivative)
    # Communicate the boundaries
    # Apply needed SATs
    applySATs(W,R,D,SolutionMode)
    # r = u - Δt Du
    @. W = R - D.SC.θ*D.SC.Δt*W #(I - θΔtD⟂)u
end
function Awonky!(read::Symbol,D::newLocalDataBlock{TT,DIM,AT,KT,DCT,GT,BT,DT,PT}) where {TT,DIM,AT,KT,DCT,GT,BT,DT,PT}
    # Compute the derivatives
    W = getarray(D,:cache)
    R = getarray(D,read)
    
    mul!(W,R,D.K,D.Derivative)      # D₂u
    applySATs(W,R,D,SolutionMode)   # D⟂ = D₂u + SATu

    @. W = R + (1-D.SC.θ)*D.SC.Δt*W # u + (1-θ)ΔtD⟂u
    addSource!(D.source,W,D.grid,D.SC.t,D.SC.Δt,D.SC.θ) # + source
    applySATs(W,D,DataMode) # + RHS SATs
    @. R = W
end

function A!(source::Symbol,DB::DataMultiBlock{TT}) where {TT}
    fillBuffers(source,DB)
    for i in eachblock(DB)
        A!(source,DB[i])
    end
end


"""
    implicit_euler
Implicit Euler with conjugate gradient method
"""
function implicit_euler(DBlock::DataMultiBlock,t::TT,Δt::TT) where TT
    # advance time
    for I in eachblock(DBlock)
        # DBlock[I].SC.t += DBlock[I].SC.Δt
        DBlock[I].SC.t = t
        DBlock[I].SC.Δt = Δt
    end
    setBoundaryConditions!(DBlock)
    # update diff coefficients
    setDiffusionCoefficient!(DBlock)
    # Set and update boundary data
    applySATs(:b,DBlock,DataMode)
    # add source terms
    addSource!(:b,DBlock)
    # Compute uₙ₊₁
    conj_grad!(DBlock)
end


"""
    theta_method
Use the Crank-Nicholson method to solve the PDE
"""
function theta_method(DBlock::DataMultiBlock,t::TT,Δt::TT) where TT
    # b = u_n already set; SATs applied after so compute 
    # b = (I - Δt/2 D_\perp)uₙ
    for I in eachblock(DBlock)
        DBlock[I].SC.t = t
        DBlock[I].SC.Δt = Δt
    end

    setDiffusionCoefficient!(DBlock)
    setBoundaryConditions!(DBlock)
    
    # b += Δt (SAT_{data}^{n+1} + SAT_{data}^{n}) + Δt/2 (S^{n+1} + S^{n})
    for I in eachblock(DBlock)
        Awonky!(:b,DBlock[I])
    end
    
    
    # Compute uₙ₊₁ = (I - D)^{-1} b
    conj_grad!(DBlock)

    for I in eachblock(DBlock)
        DBlock[I].SC.t += DBlock.SC.Δt
        # DBlock[I].SC.Δt = DBlock.SC.Δt*DBlock.SC.θ
    end
end






function RHS!(cache,source,D)
    
    if (θ - 1) != TT(0)
        mul!(cache,source,D.K,D.Derivative)
        @. cache = u + (TT(1)-DBlock.SC.θ)*DBlock[I].SC.Δt*cache    # I - Δt/2 D₂u
    end

    applySATs(cache,source,D,SolutionMode)

    cache = source
end

