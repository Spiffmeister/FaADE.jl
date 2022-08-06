
"""
function Dₓ(u::AbstractVector{Float64},n::Int64,Δx::Float64;order::Int64=2)
or
Dₓ(u::Matrix{Float64},n::Vector{Int64},Δx::Float64;order::Int64=2,dims::Union{Int64,Vector{Int64}}=1)

- If `typeof(u) <: AbstractVector{Float64}` then uses a 1D first derivative SBP operator.
- Returns a `Vector{Float64}` of length `n`.
- If `typeof(u) <: AbstractMatrix{Float64}` then runs over a matrix.
- If `dims==1` then takes derivative along rows (`u[:,i]`)
- If `dims==2` then takes derivative along columns (`u[i,:]`)
- Returns a `Matrix{Float64}` of size `nx × ny`

Internally uses [`FirstDerivative`](@ref)
"""
function Dₓ end
function Dₓ(u::AbstractVector{Float64},n::Int64,Δx::Float64;order::Int64=2)
    #= First derivative for 1D systems =#
    # Computes the first derivative for 1D systems, the default order is 2
    uₓ = zeros(Float64,n)
    Dₓ!(uₓ,u,n,Δx,order=order)

    return uₓ
end
function Dₓ(u::Matrix{Float64},nx::Int64,ny::Int64,Δx::Float64;order::Int64=2,dims::Union{Int64,Vector{Int64}}=1)
#= First derivative for multidimensional systems =#

    uₓ = zeros(Float64,nx,ny)
    if dims == 1
        # Derivative in the 1st dimension
        for i = 1:ny
            uₓ[:,i] = Dₓ!(uₓ[:,i],u[:,i],nx,Δ,order=order)
        end
    elseif dims == 2
        # Derivative in the 2nd dimension
        for i = 1:nx
            uₓ[i,:] = Dₓ!(uₓ[i,:],u[i,:],ny,Δ,order=order)
        end
    else
        error("dim must be 1 or 2.")
    end

    return u
end

