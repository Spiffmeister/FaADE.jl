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
In-place conjugate gradient method. Designed for multiblock problems

See also [`build_H`](@ref), [`A!`](@ref), [`innerH`](@ref)
"""
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
    # unorm = max(sqrt(innerprod(:u,:u,DBlock)),1e-14) #√(uₙ₊₁,uₙ₊₁)ₕ
    bnorm = max(sqrt(innerprod(:b,:b,DBlock)),1e-14) #√(b,b)ₕ
    while (rnorm ≥ rtol*bnorm) & (i < maxIT)
        # testA!(:dₖ,DBlock[1]) # cache ← Adₖ
        A!(:dₖ,DBlock) # cache ← Adₖ = dₖ - θ*Δt*D(dₖ)
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
    if (rnorm>rtol*bnorm) & warnings
        DBlock.SC.converged = false
        warnstr = string("CG did not converge at t=",DBlock.SC.t," with Δt=",DBlock.SC.Δt," i=",i," rel error=",rnorm/bnorm,", rel tolerance=",rtol,".")
        @warn warnstr
    end 
    # warnstr = string("CG did not converge at t=",DBlock.SC.t," with Δt=",DBlock.SC.Δt," i=",i," rel error=",rnorm/bnorm,", rel tolerance=",rtol,".")
    # @warn warnstr
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

function A!(read::Symbol,D::newLocalDataBlock{TT,DIM,COORD,AT,KT,DCT,GT,BT,DT,ST,PT}) where {TT,DIM,COORD,AT,KT,DCT,GT,BT,DT,ST,PT}
    # Compute the derivatives
    W = getarray(D,:cache)
    R = getarray(D,read)

    mul!(W,R,D.K,D.Derivative,TT(0)) # W ← D(u)
    applySATs(W,R,D,SolutionMode) # W += SATu
    
    if COORD == :Variable
        @. W = W/D.grid.J
    end
    
    @. W = R - D.SC.θ*D.SC.Δt*W #(I - θΔtD⟂)u

    W
end
function A!(source::Symbol,DB::DataMultiBlock{TT}) where {TT}
    fillBuffers(source,DB)
    for i in eachblock(DB)
        A!(source,DB[i])
    end
end


function CGRHS!(D::newLocalDataBlock{TT,DIM,COORD,AT,KT,DCT,GT,BT,DT,ST,PT}) where {TT,DIM,COORD,AT,KT,DCT,GT,BT,DT,ST,PT}
    # Compute the derivatives
    cache = getarray(D,:cache)
    b = getarray(D,:b)
    u = getarray(D,:u)
    mul!(cache,u,D.K,D.Derivative,TT(0))      # D₂u
    applySATs(cache,u,D,SolutionMode)   # D⟂ = D₂u + SATu
    if COORD == :Constant
        @. b = u + (1-D.SC.θ)*D.SC.Δt*cache # u + (1-θ)ΔtD⟂u
    else
        @. b = u*D.grid.J + (1-D.SC.θ)*D.SC.Δt*cache # u + (1-θ)ΔtD⟂u
    end
    # setBoundaryConditions!(D)
    applySATs(b,D,DataMode) # + RHS SATs

    if COORD == :Variable
        @. b = b/D.grid.J
    end
    
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
    # map(setDiffusionCoefficient!,DBlock.Block)
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

