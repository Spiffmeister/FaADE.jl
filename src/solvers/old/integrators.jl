#=
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
=#


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
    atol=1.e-5,rtol=1.e-12,maxIT::Int=1000,warnings=false) where {TT,DIM}

    local rnorm::TT
    local unorm::TT
    local dₖAdₖ::TT
    local βₖ::TT
    local αₖ::TT

    # testA!(:u,DBlock[1]) # cache ← Au
    A!(:u,DBlock) # cache ← Au = u - θ*Δt*Du
    setValue(:rₖ,:cache,DBlock) # r ← cache
    muladd!(:rₖ,:b,DBlock,α=TT(-1)) # r = b - cache
    setValue(:dₖ,:rₖ,DBlock) # d = r
    
    i = 0
    rnorm = sqrt(innerprod(:rₖ,:rₖ,DBlock)) #√(rₖ,rₖ)ₕ
    unorm = max(sqrt(innerprod(:u,:u,DBlock)),1e-14) #√(uₙ₊₁,uₙ₊₁)ₕ
    while (rnorm ≥ rtol*unorm) & (i < maxIT)
        # testA!(:dₖ,DBlock[1]) # cache ← Adₖ
        A!(:dₖ,DBlock[1]) # cache ← Adₖ = dₖ - θ*Δt*D(dₖ)
        dₖAdₖ = innerprod(:dₖ,:cache,DBlock) #αₖ = (dₖ,Ddₖ)
        αₖ = rnorm^2 / dₖAdₖ

        # @. DBlock[1].uₙ₊₁ = DBlock[1].uₙ₊₁ + αₖ*DBlock[1].dₖ
        muladd!(:uₙ₊₁,:dₖ,DBlock,β=αₖ)
        
        # testA!(:uₙ₊₁,DBlock[1]) # cache ← Auₙ₊₁ = uₙ₊₁ - θ*Δt*Duₙ₊₁
        # A!(:uₙ₊₁,DBlock) # cache ← Auₙ₊₁ = uₙ₊₁ - θ*Δt*Duₙ₊₁
        # DBlock[1].rₖ .= DBlock[1].b - DBlock[1].cache #rₖ = b - cache

        # @. DBlock[1].rₖ = DBlock[1].rₖ - αₖ*DBlock[1].cache #rₖ = rₖ - αₖAdₖ
        muladd!(:rₖ,:cache,DBlock,β=-αₖ) #rₖ = rₖ - αₖAdₖ

        βₖ = innerprod(:rₖ,:rₖ,DBlock) #βₖ*(rₖ₋₁,rₖ₋₁) = (rₖ,rₖ)

        # @. DBlock[1].dₖ = DBlock[1].rₖ + βₖ/rnorm^2 * DBlock[1].dₖ #dₖ = rₖ + βₖ dₖ
        muladd!(:dₖ,:rₖ,DBlock,α=βₖ/rnorm^2) #dₖ = rₖ + βₖ/rnorm^2 * dₖ
        
        rnorm = sqrt(βₖ)
        # println(rnorm)
        i += 1
    end
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

function A!(read::Symbol,D::newLocalDataBlock{TT,DIM,AT,KT,DCT,GT,BT,DT,PT}) where {TT,DIM,AT,KT,DCT,GT,BT,DT,PT}
    # Compute the derivatives
    W = getarray(D,:cache)
    R = getarray(D,read)

    mul!(W,R,D.K,D.Derivative) # W ← D(u)
    applySATs(W,R,D,SolutionMode) # W += SATu
    
    @. W = R - D.SC.θ*D.SC.Δt*W #(I - θΔtD⟂)u

    W
end
function A!(source::Symbol,DB::DataMultiBlock{TT}) where {TT}
    fillBuffers(source,DB)
    for i in eachblock(DB)
        A!(source,DB[i])
    end
end


function CGRHS!(D::newLocalDataBlock{TT,DIM,AT,KT,DCT,GT,BT,DT,PT}) where {TT,DIM,AT,KT,DCT,GT,BT,DT,PT}
    # Compute the derivatives
    cache = getarray(D,:cache)
    b = getarray(D,:b)
    u = getarray(D,:u)
    mul!(cache,u,D.K,D.Derivative)      # D₂u
    applySATs(cache,u,D,SolutionMode)   # D⟂ = D₂u + SATu
    @. b = u + (1-D.SC.θ)*D.SC.Δt*cache # u + (1-θ)ΔtD⟂u
    # setBoundaryConditions!(D)
    applySATs(b,D,DataMode) # + RHS SATs
    
    addSource!(D.source,b,D.grid,D.SC.t,D.SC.Δt,D.SC.θ) # + source
end
function CGRHS!(DB::DataMultiBlock{TT}) where {TT}
    fillBuffers(:u,DB)
    for i in eachblock(DB)
        CGRHS!(DB[i])
    end
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
    
    # b = Δt (SAT_{data}^{n+1} + SAT_{data}^{n}) + Δt/2 (S^{n+1} + S^{n})
    CGRHS!(DBlock)
    # Compute uₙ₊₁ = (I - D)^{-1} b
    conj_grad!(DBlock)
end




function testA!(read::Symbol,D::newLocalDataBlock)
    #uses whatever is stored in K[1] as the A operator
    R = getarray(D,read)
    D.cache = D.K*R
end

