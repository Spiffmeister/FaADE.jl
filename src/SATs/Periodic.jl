
"""
    SAT_Periodic(u::Vector{Float64},Δx::Float64;c=[1.0,1.0],order::Int64=2)

Simulatenous approximation term for Periodic boundary conditions
"""
function SAT_Periodic(u::AbstractVector,Δx::Float64;
        c=[1.0,1.0],order::Int64=2)

    # Get h value
    h = hval(order)

    # Penalties
    α₀ = 0.5/(h * Δx) # Derivative penatly
    τ₁ = -0.5/(h * Δx) # Symmeteriser penalty
    τ₀ = -max(c[1]/2(h*Δx),c[end]/2(h*Δx))/(h*Δx) # Dirichlet penalty

    SAT = zeros(Float64,2order)

    # Dirichlet terms
    SAT[1]  += τ₀ * (u[1] - u[end])
    SAT[end]+= τ₀ * (u[end] - u[1])

    # Symmeteriser
    L₁u = [(u[1] - u[end])]
    SAT[1:order]        = τ₁ * c[1] * BDₓᵀ(L₁u,Left,Δx,order)
    SAT[order+1:2order] = -τ₁ * c[end] * BDₓᵀ(L₁u,Right,Δx,order)

    # Neumann term
    BDₓu₁ = BDₓ(u,Left,Δx,order)
    BDₓuₙ = BDₓ(u,Right,Δx,order)

    SAT[1]  += α₀ * (c[1]*BDₓu₁ - c[end]*BDₓuₙ)
    SAT[end]+= α₀ * (c[1]*BDₓu₁ - c[end]*BDₓuₙ)

    return SAT[1:order], SAT[order+1:end]

end


function SAT_Periodic!(uₓₓ::AbstractVector,u::AbstractVector,Δx::Float64;
        c=[1.0,1.0],order::Int64=2)

    # Get h value
    h = hval(order)

    # Penalties
    α₀ = 0.5/(h * Δx) # Derivative penatly
    τ₁ = -0.5/(h * Δx) # Symmeteriser penalty
    τ₀ = -max(c[1]/2(h*Δx),c[end]/2(h*Δx))/(h*Δx) # Dirichlet penalty

    uₓₓ[1]  += τ₀ * (u[1] - u[end])
    uₓₓ[end]+= τ₀ * (u[end] - u[1])

    # Symmeteriser
    L₁u = (u[1] - u[end])

    uₓₓ[1:order]        .+= τ₁ * c[1] * BDₓᵀ(L₁u,Left,Δx,order)
    uₓₓ[end-order+1:end].+= τ₁ * c[end] * BDₓᵀ(L₁u,Right,Δx,order)

    # Neumann term
    BDₓu₁ = BDₓ(u,Left,Δx,order)
    BDₓuₙ = BDₓ(u,Right,Δx,order)

    uₓₓ[1]  += α₀ * (c[1]*BDₓu₁ - c[end]*BDₓuₙ)
    uₓₓ[end]+= α₀ * (c[1]*BDₓu₁ - c[end]*BDₓuₙ)
    return uₓₓ

end



# function SAT_Periodic!(uₓₓ::AbstractMatrix,u::AbstractVector,Δx::Float64,n;
#         c=[1.0,1.0],order::Int64=2)
    
#     for i = 1:n
#     end
    
# end