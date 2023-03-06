#= EXPLICIT METHODS =#

"""
    forward_euler
"""
function forward_euler end
function forward_euler(DBlock::DataBlock,RHS::Function)
    # Simple forward euler method
    RHS(DBlock.uₙ₊₁,DBlock.u,DBlock.K)
    DBlock.uₙ₊₁ .= DBlock.u .+ DBlock.Δt*DBlock.uₙ₊₁
end
function forward_euler(uₒ::AbstractMatrix,RHS::Function,Δt::Float64)
    # Simple forward euler method
    uₙ = uₒ + Δt*RHS(uₙ,uₒ,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,gx,gy)
    return uₙ
end


"""
    RK4
"""
mutable struct ExplicitBlock{T,N} <: DataBlockType{T,N}
    k1  :: AbstractArray{T,N}
    k2  :: AbstractArray{T,N}
    k3  :: AbstractArray{T,N}
    k4  :: AbstractArray{T,N}
    Δt  :: T

    function ExplicitBlock{T}(grid::GridType,Δt::T,integrator::Symbol=:RK4) where T
        if typeof(grid) <: Grid1D
            n = grid.n
        elseif typeof(grid) <: Grid2D
            n = (grid.nx,grid.ny)
        end

        dims = length(n)

        k1 = zeros(T,n)
        if integrator == :RK4
            k2 = zeros(T,n)
            k3 = zeros(T,n)
            k4 = zeros(T,n)
        elseif integrator == :euler
            k2 = k3 = k4 = zeros(T,0)
        end

        new{T,dims}(k1,k2,k3,k4,Δt)
    end
end
# function RK4 end
# function RK4(uₙ::Vector,uₒ::Vector,RHS::Function,n::Int,Δx::Float64,Δt::Float64,k::Vector,t::Float64,x::Vector,boundary)
function (RK::ExplicitBlock)(RHS::Function,DBlock::DataBlock,t::Float64)
    
    DBlock.uₙ₊₁ .= DBlock.u

    RHS(RK.k1,DBlock.u,         DBlock.K,t)
    RHS(RK.k2,DBlock.u+0.5RK.k1,DBlock.K,t+0.5RK.Δt)
    RHS(RK.k3,DBlock.u+0.5RK.k2,DBlock.K,t+0.5RK.Δt)
    RHS(RK.k4,DBlock.u+RK.k3,   DBlock.K,t+RK.Δt)
    DBlock.uₙ₊₁ .+= DBlock.u + RK.Δt/6.0 * (RK.k1 + 2RK.k2 + 2RK.k3 + RK.k4)

end
function (ForwardEuler::ExplicitBlock)(RHS::Function,DBlock::DataBlock,t::Float64)
    
    DBlock.uₙ₊₁ .= DBlock.u

    RHS(RK.k1,DBlock.u,DBlock.K,t)
    DBlock.uₙ₊₁ .+= DBlock.u + RK.Δt*RK.k1

end



#= IMPLICIT METHODS =#
"""
    implicit_euler
"""
function implicit_euler end
# Vector form
function implicit_euler(uₙ::Vector,uₒ::Vector,RHS::Function,n::Int,Δx::Float64,Δt::Float64,k::Vector,t::Float64,x::Vector,boundary;maxIT::Int=100,α::Float64=1.5)
    uⱼ = uₒ
    for j = 1:maxIT
        uⱼ = uⱼ - α * Δt * (uⱼ - uₒ - Δt*RHS(uₙ,uⱼ,n,x,Δx,t,Δt,k,boundary))
    end
    uⱼ = uₒ + Δt*RHS(uₙ,uⱼ,n,x,Δx,t,Δt,k,boundary)
    return uⱼ
end
# Matrix form
function implicit_euler(uₙ::Matrix,uₒ::Matrix,RHS::Function,nx::Int,ny::Int,Δx::Float64,Δy::Float64,Δt::Float64,kx::Matrix,ky::Matrix,t::Float64,x::Vector,y::Vector,boundary_x,boundary_y;maxIT::Int=100,α::Float64=1.5)
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
function conj_grad!(RHS::Function,DBlock::DataBlock,CGB::ConjGradBlock,Δt::Float64;
        atol::Float64=1.e-5,rtol::Float64=1.e-10,maxIT::Int=10,warnings=true)
    
    # x₀ = uₙ #Initial guess
    # CGB.b .= DBlock.uₙ₊₁ #uₙ₊₁ is our initial guess and RHS
    # rₖ = (uₙ₊₁ - Δt*uₓₓⁿ⁺¹) - (uₙ + F)
    A!(RHS,CGB.rₖ,DBlock.u,Δt,DBlock.K)
    CGB.rₖ .= CGB.rₖ .- CGB.b
    # DBlock.uₙ₊₁ .= DBlock.u

    CGB.dₖ .= -CGB.rₖ # dₖ = -rₖ

    
    i = 0
    rnorm = sqrt(CGB.innerprod(CGB.rₖ,CGB.rₖ))
    unorm = sqrt(CGB.innerprod(DBlock.uₙ₊₁,DBlock.uₙ₊₁))
    while (rnorm > rtol*unorm) & (i < maxIT)
        A!(RHS,CGB.Adₖ,CGB.dₖ,Δt,DBlock.K) # Adₖ = dₖ - Δt*D(dₖ)
        dₖAdₖ = CGB.innerprod(CGB.dₖ,CGB.Adₖ)
        αₖ = -CGB.innerprod(CGB.rₖ,CGB.dₖ)/dₖAdₖ
        DBlock.uₙ₊₁ .= DBlock.uₙ₊₁ .+ αₖ*CGB.dₖ #xₖ₊₁ = xₖ + αₖ*dₖ

        # rₖ = (xₖ₊₁ - Δt*Dxₖ₊₁ - b)
        A!(RHS,CGB.rₖ,DBlock.uₙ₊₁,Δt,DBlock.K)
        CGB.rₖ .= CGB.rₖ .- CGB.b

        A!(RHS,CGB.Drₖ,CGB.rₖ,Δt,DBlock.K) # Drₖ = rₖ - Δt*D(rₖ)

        βₖ = CGB.innerprod(CGB.rₖ,CGB.Drₖ)/dₖAdₖ
        CGB.dₖ .= -CGB.rₖ .+ βₖ*CGB.dₖ
        rnorm = sqrt(CGB.innerprod(CGB.rₖ,CGB.rₖ))
        i += 1
    end
    if (rnorm>rtol*unorm) & warnings
        CGB.converged = false
        warnstr = string("CG did not converge with Δt=",Δt,", rel error=",rnorm/unorm,", rel tolerance=",rtol,".")
        @warn warnstr
    end
end


#= SUPPORT =#


"""
    A!
Mutable ``u - ΔtD(u)``
"""
function A!(PDE!::Function,tmp::AbstractArray,uⱼ::AbstractArray,Δt::Float64,k::AbstractArray)
    PDE!(tmp,uⱼ,k)
    tmp .= uⱼ .- Δt*tmp
end

