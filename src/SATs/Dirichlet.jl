


struct Boundary_Dirichlet{T} <: SimultanousApproximationTerm
    BDₓᵀ    :: Vector{T}
    RHS     :: Function
    type    :: BoundaryCondition
    side    :: NodeType
    axis    :: Int
    order   :: Int
    Δx      :: T
    α       :: T
    τ       :: T
    expression  :: String

    function Boundary_Dirichlet(RHS,order,Δx,side,axis)

        BD = BoundaryDerivativeTranspose(order,Δx)



        fullsat = "τH⁻¹ E H⁻¹E(u-f) + α H⁻¹ (K H Dₓᵀ) H⁻¹ E (u-f)"


        new(BD,Dirichlet,side,axis,order,Δx,α,τ,fullsat)

    end

end


"""
    SAT_Dirichlet

1. SAT_Dirichlet(::NodeType{:Left},u::AbstractVector,Δ::Float64;
        c=1.0,order::Int64=2,forcing=false)
2. SAT_Dirichlet(::NodeType{:Left,:Right},u::Vector{Float64},Δx::Float64,RHS;
        c=1.0,order::Int64=2,separate_forcing::Bool=false)
 

Simulatenous approximation term for Dirichlet boundary conditions 
    ``u(xᵢ)=g``
where ``i\\in \\{0,N\\}``.

`NodeType` is either `Left` or `Right`.

See [`NodeType`](@ref)
"""
function SAT_Dirichlet end
# Dirichlet for implicit integrators
function SAT_Dirichlet(node,u::AbstractArray,Δ::Float64; c=1.0,order::Int64=2,forcing=false)
    SAT = zeros(Float64,order)
    SAT = SAT_Dirichlet!(SAT,node,u,Δ,c=c,order=order,forcing=forcing)
end
function SAT_Dirichlet(node,u::AbstractVector,Δ::Float64,RHS; c=1.0,order::Int64=2,forcing=false)
    SAT = zeros(Float64,order)
end
# Dirichlet for explicit integrators




function (SAT::Boundary_Dirichlet)()




"""
    SAT_Dirichlet!


"""
function SAT_Dirichlet!(SAT::AbstractVector,::NodeType{:Left},u::AbstractVector,Δ::Float64,RHS;
        c=1.0,order::Int64=2,forcing::Bool=false)
    α,τ = SATpenalties(Dirichlet,Δ,order)

    SAT[1:order] += BDₓᵀ(u,Left,Δ,order) - BDₓᵀ(RHS,Left,Δ,order)

    SAT[1:order] .*= α * c[1]
    SAT[1]  += τ*(u[1] - RHS[1])
end
function SAT_Dirichlet(SAT::AbstractVector,::NodeType{:Right},u::AbstractVector,Δ::Float64,RHS;
        c=1.0,order::Int64=2,forcing::Bool=false)
    α,τ = SATpenalties(Dirichlet,Δ,order)

    SAT[end-order+1:end] .+= BDₓᵀ(u,Right,Δ,order) + BDₓᵀ(-RHS,Right,Δ,order)

    SAT[end-order+1:end] .*= α * c[end]
    SAT[end]  += τ*(u[end] - RHS[end])
end




#=== Explicit methods ===#
function SAT_Dirichlet_internal!(SAT::AbstractArray,::NodeType{:Left},u::AbstractArray,c::AbstractArray,RHS::Float64,Δ::Float64,α::Float64,τ::Float64,BD::AbstractVector,order::Int,loopaxis::Function)
    for (S,C,U) in zip(loopaxis(SAT),loopaxis(c),loopaxis(u))
        S[1:order] .+= α*C[1]*(BD*U[1] .- BD*RHS)
        S[1] += τ*(U[1] - RHS)
    end
end
function SAT_Dirichlet_internal!(SAT::AbstractArray,::NodeType{:Right},u::AbstractArray,c::AbstractArray,Δ::Float64,RHS::Float64,α::Float64,τ::Float64,BD::AbstractVector,order::Int,loopaxis::Function)
    for (S,C,U) in zip(loopaxis(SAT),loopaxis(c),loopaxis(u))
        S[end-order+1:end] .+= α*C[end]*(BD*U[end] .- BD*RHS)
        S[end] += τ*(U[end] - RHS)
    end
end

#=== Implicit methods ===#
function SAT_Dirichlet_internal!(SAT::AbstractArray,::NodeType{:Left},u::AbstractArray,c::AbstractArray,Δ::Float64,α::Float64,τ::Float64,BD::AbstractVector,order::Int,loopaxis::Function)
    for (S,C,U) in zip(loopaxis(SAT),loopaxis(c),loopaxis(u))
        S[1:order] .+= α*C[1]*BD*U[1]
        S[1] += τ*U[1]
    end
end
function SAT_Dirichlet_internal!(SAT::AbstractArray,::NodeType{:Right},u::AbstractArray,c::AbstractArray,Δ::Float64,α::Float64,τ::Float64,BD::AbstractVector,order::Int,loopaxis::Function)
    for (S,C,U) in zip(loopaxis(SAT),loopaxis(c),loopaxis(u))
        S[end-order+1:end] .+= α*C[end]*BD*U[end]
        S[end] += τ*U[end]
    end
end
function SAT_Dirichlet_internal_forcing!(SAT::AbstractArray,::NodeType{:Left},u::AbstractArray,c::AbstractArray,Δ::Float64,α::Float64,τ::Float64,BD::AbstractVector,order::Int,loopaxis::Function)
    for (S,C) in zip(loopaxis(SAT),loopaxis(c))
        S[1:order] .+= α*C[1]*BD*u[1]
        S[1] += τ*u[1]
    end
end
function SAT_Dirichlet_internal_forcing!(SAT::AbstractArray,::NodeType{:Right},u::AbstractArray,c::AbstractArray,Δ::Float64,α::Float64,τ::Float64,BD::AbstractVector,order::Int,loopaxis::Function)
    for (S,C) in zip(loopaxis(SAT),loopaxis(c))
        S[end-order+1:end] .+= α*C[end]*BD*u[end]
        S[end] += τ*u[end]
    end
end





function SAT_Dirichlet!(uₓₓ::AbstractArray,::NodeType{:Right},u::AbstractArray,Δ::Float64,α,τ,BD;
        c=1.0,order::Int64=2,forcing=false)

    # α,τ = SATpenalties(Dirichlet,Δ,order)

    if !forcing
        uₓₓ[end-order+1:end] .+= α*c[end]*[-u[1], u[1]]/Δ
        uₓₓ[end] += τ*u[end]
    else
        uₓₓ[end-order+1:end] .-= α*c[end]*[-u[1], u[1]]/Δ
        uₓₓ[end] -= τ*u[end]
    end
end














