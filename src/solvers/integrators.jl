#= EXPLICIT METHODS =#

"""
    forward_euler(uₙ::Vector,uₒ::Vector,RHS::Function,n::Int,Δx::Float64,Δt::Float64,k::Vector,t::Float64,x::Vector,g)
or
    function forward_euler(DBlock::DataBlock,RHS::Function)
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
    RK4(uₙ::Vector,uₒ::Vector,RHS::Function,n::Int,Δx::Float64,Δt::Float64,k::Vector,t::Float64,x::Vector,boundary)
"""
function RK4 end
function RK4(uₙ::Vector,uₒ::Vector,RHS::Function,n::Int,Δx::Float64,Δt::Float64,k::Vector,t::Float64,x::Vector,boundary)
    k1 = RHS(uₙ,uₒ        ,n,x,Δx,t,Δt,k,       boundary)
    k2 = RHS(uₙ,uₒ+Δt/2*k1,n,x,Δx,t+0.5Δt,Δt,k, boundary)
    k3 = RHS(uₙ,uₒ+Δt/2*k2,n,x,Δx,t+0.5Δt,Δt,k, boundary)
    k4 = RHS(uₙ,uₒ+Δt*k3  ,n,x,Δx,t+Δt,Δt,k,    boundary)
    uₙ = uₒ + Δt/6 * (k1 + 2k2 + 2k3 + k4)
    return uₙ
end



#= IMPLICIT METHODS =#
"""
    implicit_euler(uₙ::Vector,uₒ::Vector,RHS::Function,n::Int,Δx::Float64,Δt::Float64,k::Vector,t::Float64,x::Vector,boundary;maxIT::Int=100,α::Float64=1.5)
or
    function implicit_euler(uₙ::Matrix,uₒ::Matrix,RHS::Function,nx::Int,ny::Int,Δx::Float64,Δy::Float64,Δt::Float64,kx::Matrix,ky::Matrix,t::Float64,x::Vector,y::Vector,boundary_x,boundary_y;maxIT::Int=100,α::Float64=1.5)
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
function conj_grad!(DBlock::DataBlock,CGB::ConjGradBlock,RHS::Function,Δt::Float64,order::Int;
        atol::Float64=1.e-5,rtol::Float64=1.e-10,maxIT::Int=10,warnings=true)
    
    # x₀ = uₙ #Initial guess
    # CGB.b .= DBlock.uₙ₊₁ #uₙ₊₁ is our initial guess and RHS
    # rₖ = (uₙ₊₁ - Δt*uₓₓⁿ⁺¹) - (uₙ + F)
    A!(CGB.rₖ,DBlock.u,RHS,Δt,DBlock.K)
    CGB.rₖ .= CGB.rₖ .- CGB.b
    # DBlock.uₙ₊₁ .= DBlock.u

    CGB.dₖ .= -CGB.rₖ # dₖ = -rₖ
    
    i = 0
    rnorm = sqrt(innerH(CGB.rₖ,CGB.rₖ,order,CGB.Δ))
    unorm = sqrt(innerH(DBlock.uₙ₊₁,DBlock.uₙ₊₁,order,CGB.Δ))
    while (rnorm > rtol*unorm) & (i < maxIT)
        A!(CGB.Adₖ,CGB.dₖ,RHS,Δt,DBlock.K) # Adₖ = dₖ - Δt*D(dₖ)
        dₖAdₖ = innerH(CGB.dₖ,CGB.Adₖ, order,CGB.Δ)
        αₖ = -innerH(CGB.rₖ,CGB.dₖ, order,CGB.Δ)/dₖAdₖ

        DBlock.uₙ₊₁ .= DBlock.uₙ₊₁ .+ αₖ*CGB.dₖ #xₖ₊₁ = xₖ + αₖ*dₖ

        # rₖ = (xₖ₊₁ - Δt*Dxₖ₊₁ - b)
        A!(CGB.rₖ,DBlock.uₙ₊₁,RHS,Δt,DBlock.K)
        CGB.rₖ .= CGB.rₖ .- CGB.b

        A!(CGB.Drₖ,CGB.rₖ,RHS,Δt,DBlock.K) # Drₖ = rₖ - Δt*D(rₖ)

        βₖ = innerH(CGB.rₖ,CGB.Drₖ, order,CGB.Δ)/dₖAdₖ
        CGB.dₖ .= -CGB.rₖ .+ βₖ*CGB.dₖ
        rnorm = sqrt(innerH(CGB.rₖ,CGB.rₖ, order,CGB.Δ))
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
function A!(tmp::AbstractArray,uⱼ::AbstractArray,PDE!::Function,Δt::Float64,k::AbstractArray)
    PDE!(tmp,uⱼ,k)
    tmp .= uⱼ .- Δt*tmp
end


"""
    innerH
Computes the ``H``-inner product between two vectors or matrices ``u`` and ``v``.
"""
function innerH end
# H inner product for 1D problems
@views function innerH(u::AbstractVector,v::AbstractVector,order::Int,Δ::Float64)
    tmp = 0.0
    if order == 2
        tmp += 0.5*Δ * u[1]*v[1]
        tmp += 0.5*Δ * u[end]*v[end]
        tmp += Δ*dot(u[2:end-1],v[2:end-1])
        return tmp
    elseif order == 4
        tmp += sum([17.0/48.0, 59.0/48.0, 43.0/48.0, 49.0/48.0]*Δ .* u[1:4] .* v[1:4])
        tmp += sum([49.0/48.0, 43.0/48.0, 59.0/48.0, 17.0/48.0]*Δ .* u[end-3:end] .* v[end-3:end])
        tmp += Δ*dot(u[5:end-4],v[5:end-4])
        return tmp
    elseif order == 6
        tmp += sum([13649.0/43200.0, 12013.0/8640.0, 2711.0/4320.0, 5359.0/4320.0, 7877.0/8640.0, 43801.0/43200.0])
        tmp += sum([43801.0/43200.0, 7877.0/8640.0, 5359.0/4320.0, 2711.0/4320.0, 12013.0/8640.0, 13649.0/43200.0])
        tmp += Δ*dot(u[7:end-6],v[7:end-6])
    end
end
# H inner product for 2D problems
@views function innerH(u::AbstractMatrix,v::AbstractMatrix,order::Int,Δ::Float64)
    tmp = 0.0
    if order == 2
        tmp += 0.25*Δ * u[1,1]*v[1,1]
        tmp += 0.25*Δ * u[end,end]*v[end,end]
        tmp += 0.25*Δ * u[1,end]*v[1,end]
        tmp += 0.25*Δ * u[end,1]*v[end,1]

        tmp += 0.5*Δ * dot(u[2:end-1,1],v[2:end-1,1])
        # tmp += 0.5*Δ * dot(@views(u[2:end-1,1]),@views(v[2:end-1,1]))
        tmp += 0.5*Δ * dot(u[2:end-1,end],v[2:end-1,end])
        tmp += 0.5*Δ * dot(u[1,2:end-1],v[1,2:end-1])
        tmp += 0.5*Δ * dot(u[end,2:end-1],v[end,2:end-1])

        # Internal sum
        # for (uᵢ,vⱼ) in zip(u[2:end-1,2:end-1],v[2:end-1,2:end-1]) #TODO: Write this in nifty.jl
        #     tmp += Δ * uᵢ*vⱼ
        # end
        tmp += Δ * dot(u[2:end-1,2:end-1],v[2:end-1,2:end-1])  #see LinearAlgebra.dot
    elseif order == 4
    elseif order == 6
    end
    return tmp
end

"""
    build_H(n::Int64,order::Int64)
"""
function build_H(n::Int64,order::Int64)
    H = ones(n)
    if order == 2
        H[1] = H[end] = 0.5
    elseif order == 4
        H[1] = H[end] = 17.0/48.0
        H[2] = H[end-1] = 59.0/48.0
        H[3] = H[end-2] = 43.0/48.0
        H[4] = H[end-3] = 49.0/48.0
    elseif order == 6
        H[1] = H[end] = 13649.0/43200.0
        H[2] = H[end-1] = 12013.0/8640.0
        H[3] = H[end-2] = 2711.0/4320.0
        H[4] = H[end-3] = 5359.0/4320.0
        H[5] = H[end-4] = 7877.0/8640.0
        H[6] = H[end-5] = 43801.0/43200.0
    else
        error("Order must be 2,4 or 6")
    end
    return H
end

function innerH(u::AbstractMatrix,Hx::AbstractVector,Hy::AbstractVector,v::AbstractMatrix)
    #= DEPRECATED =#
    nx,ny = size(u)
    tmp = 0.0
    for j = 1:ny
        for i = 1:nx
            tmp += u[i,j]*Hx[i]*Hy[j]*v[i,j]
        end
    end
    return tmp
end

