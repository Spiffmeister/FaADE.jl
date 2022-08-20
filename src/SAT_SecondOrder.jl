#=
    Second derivative boundary operators
=#

"""
    SAT(type::BoundaryCondition,::NodeType{:Left},u::AbstractVector{Float64},Δx::Float64,g;
        order=2::Int,c::Union{Float64,AbstractVector{Float64}}=1.0,αβ::Vector{Float64}=[1.0,1.0],separate_forcing::Bool=false)
Returns the left SAT (or SAT+Forcing) term as a vector or SAT and Forcing term as separate vectors
Inputs: 
- type: Dirichlet, Neumann, Robin (see ? BoundaryCondition)
- ::NodeType: Left or Right (see ? NodeType)
- u: Solution vector
- Δx: step size
- g: boundary value
Optional inputs:
- order: 2 (default), 4, 6
- coefficient: c in ``\\frac{\\partial}{\\partial x}\\left(c\\frac{\\partial u}{\\partial x}\\right)``
- `αβ`: `[α,β]` Robin boundary coefficient ``\\alpha \\partial_x u + \\beta u = f``
- separate_forcing: true, false (default) -- controls the return type.
    false; returns a single vector of the SAT
    true; returns SAT and F(orcing) vectors individually
For periodic conditions call Periodic(u,Δx,c;order), for a split domain call Split_domain(u⁻,u⁺,Δx⁻,Δx⁺,c⁻,c⁺;order=2,order⁻=2,order⁺=2)
"""
function SAT1 end
function SAT1(type::BoundaryCondition,::NodeType{:Left},u::AbstractVector{Float64},Δx::Float64,g;
        order=2::Int,c::Union{Float64,AbstractVector{Float64}}=1.0,αβ::Vector{Float64}=[1.0,1.0],separate_forcing::Bool=false)
    if type == Dirichlet # Dirichlet boundaries
        SAT,F = SAT_Dirichlet1(Left,u,Δx,g,c=c,order=order,separate_forcing=true)
    elseif type == Neumann # Neumann boundary conditions
        SAT,F = SAT_Neumann1(Left,u,Δx,g,c=c,order=order,separate_forcing=true)
    elseif type == Robin # Robin boundaries
        SAT,F = SAT_Robin1(Left,u,Δx,g,order=order,a=αβ[1],b=αβ[2],separate_forcing=true)
    else
        error("Must specify either Dirichlet, Neumann or Robin. For periodic boundaries use the SAT_Periodic() function")
    end
    # If the forcing term is part of the SAT (i.e. not using GC method)
    if !separate_forcing
        SAT += F
        return SAT
    else
        return SAT, F
    end
end
function SAT1(type::BoundaryCondition,::NodeType{:Right},u::AbstractVector{Float64},Δx::Float64,g;
        order=2::Int,c::Union{Float64,AbstractVector{Float64}}=1.0,αβ::Vector{Float64}=[1.0,1.0],separate_forcing::Bool=false)
    if type == Dirichlet
        SAT, F = SAT_Dirichlet1(Right,u,Δx,g,c=c,order=order,separate_forcing=true)
    elseif type == Neumann # Neumann boundary conditions
        SAT, F = SAT_Neumann1(Right,u,Δx,g,c=c,order=order,separate_forcing=true)
    elseif type == Robin
        SAT, F = SAT_Robin1(Right,u,Δx,g,order=order,a=αβ[1],b=αβ[2],separate_forcing=true)
    else
        error("Must specify either Dirichlet, Neumann or Robin. For Periodic use the SAT_Periodic() function")
    end
    # If the forcing term is part of the SAT (i.e. not using GC method)
    if !separate_forcing
        SAT += F
        return SAT
    else
        return SAT, F
    end
end


function SAT(type::BoundaryCondition,node::NodeType,u::AbstractVector{Float64},Δ::Float64,RHS;
        order=2::Int,c::Union{Float64,AbstractVector{Float64}}=1.0,αβ::Vector{Float64}=[1.0,1.0],forcing=false)
    if type == Dirichlet
        SAT = SAT_Dirichlet(node,u,Δ,RHS,c=c,order=order)
    elseif type == Neumann
        SAT = SAT_Neumann(node,u,Δ,RHS,c=c,order=order)
    elseif type == Robin
        SAT = SAT_Robin(node,u,Δ,RHS,a=αβ[1],b=αβ[2],order=order)
    end
    return SAT
end
function SAT(type::BoundaryCondition,node::NodeType,u::AbstractVector{Float64},Δ::Float64;
        order=2::Int,c::Union{Float64,AbstractVector{Float64}}=1.0,αβ::Vector{Float64}=[1.0,1.0],forcing=false)
    if type == Dirichlet
        SAT = SAT_Dirichlet(node,u,Δ,c=c,order=order,forcing=forcing)
    elseif type == Neumann
        SAT = SAT_Neumann(node,u,Δ,c=c,order=order,forcing=forcing)
    elseif type == Robin
        SAT = SAT_Robin(node,u,Δ,a=αβ[1],b=αβ[2],order=order,forcing=forcing)
    end
    return SAT
end






### MATRIX FORM #TODO: REMOVE THIS/CHANGE TO MATRIX SAT()
function SAT_left(type::BoundaryCondition,u::AbstractMatrix{Float64},Δx::Float64,Δy::Float64,nx::Int64,ny::Int64,g,dim;
    order::Int64=2,c::Union{Float64,AbstractVector{Float64}}=1.0,αβ::Vector{Float64}=[1.0,1.0],separate_forcing::Bool=false)
    if dim == 1
        SAT = zeros(Float64,ny,order)
        F = zeros(Float64,ny,order)
        for i = 1:ny
            SAT[i,:],F[i,:] = SAT(type,u[i,:],Δx,g,order=order,c=c,αβ=αβ,separate_forcing=false)
        end
    elseif dim == 2
        SAT = zeros(Float64,nx,order)
        F = zeros(Float64,nx,order)
        for i = 1:nx
            SAT[i,:],F[i,:] = SAT(type,u[i,:],Δy,g,order=order,c=c,αβ=αβ,separate_forcing=false)
        end
    end
    if separate_forcing
        return SAT, F
    else
        return SAT + F
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
function SAT_Dirichlet(::NodeType{:Left},u::AbstractVector,Δ::Float64;
        c=1.0,order::Int64=2,forcing=false)

    α,τ = SATpenalties(Dirichlet,Δ,order)
    SAT = zeros(Float64,order)

    SAT = α*c[1]*BDₓᵀ(u,Left,Δ,order)
    SAT[1] += τ*u[1]
    if !forcing
        return SAT
    else
        return -SAT
    end
end
function SAT_Dirichlet(::NodeType{:Right},u::AbstractVector,Δ::Float64;
        c=1.0,order::Int64=2,forcing=false)

    α,τ = SATpenalties(Dirichlet,Δ,order)
    SAT = zeros(Float64,order)

    SAT = α*c[end]*BDₓᵀ(u,Right,Δ,order)
    SAT[end] += τ*u[end]
    if !forcing
        return SAT
    else
        return -SAT
    end
end
# Dirichlet for explicit integrators
function SAT_Dirichlet(::NodeType{:Left},u::AbstractVector,Δ::Float64,RHS;
        c=1.0,order::Int64=2,forcing::Bool=false)
    α,τ = SATpenalties(Dirichlet,Δ,order)
    SAT = zeros(Float64,order)

    SAT += BDₓᵀ(u,Left,Δ,order) - BDₓᵀ(RHS,Left,Δ,order)
    
    SAT .*= α * c[1]
    SAT[1]  += τ*(u[1] - RHS[1])
    
    return SAT
end
function SAT_Dirichlet(::NodeType{:Right},u::AbstractVector,Δ::Float64,RHS;
        c=1.0,order::Int64=2,forcing::Bool=false)
    α,τ = SATpenalties(Dirichlet,Δ,order)
    SAT = zeros(Float64,order)

    SAT += BDₓᵀ(u,Right,Δ,order) + BDₓᵀ(-RHS,Right,Δ,order)

    SAT .*= α * c[end]
    SAT[end]  += τ*(u[end] - RHS[end])
   
    return SAT
end


#=
====================== Neumann boundary conditions ======================
=#
"""
    SAT_Neumann

1. SAT_Neumann(::NodeType{:Left},u::Vector{Float64},Δx::Float64;
        c=1.0,order::Int64=2,forcing::Bool=false)
2. SAT_Neumann(::NodeType{:Left},u::AbstractVector,Δx::Float64,RHS;
        c=1.0,order::Int64=2,forcing::Bool=false)

Simulatenous approximation term for Neumann boundary conditions
    ``\\left.\\frac{\\partial u}{\\partial x}\\right|_{x_i} = g(t)``
where ``i\\in\\{0,N\\}``
`NodeType` is either `Left` or `Right`
"""
function SAT_Neumann end
function SAT_Neumann(::NodeType{:Left},u::Vector{Float64},Δx::Float64;
        c=1.0,order::Int64=2,forcing::Bool=false)
    # Penalties
    τ = SATpenalties(Neumann,Δx,order)
    SAT = zeros(Float64,order)

    if !forcing
        SAT[1] = τ*c[1]*BDₓ(u,Left,Δx,order)
        return SAT
    else
        SAT[1] += -τ*u[1]
        return SAT
    end
end
function SAT_Neumann(::NodeType{:Right},u::Vector{Float64},Δx::Float64;
        c=1.0,order::Int64=2,forcing::Bool=false)
    # Penalties
    τ = SATpenalties(Neumann,Δx,order)
    SAT = zeros(Float64,order)

    if !forcing
        SAT[end] = -τ*c[end]*BDₓ(u,Right,Δx,order)
        return SAT
    else
        SAT[end] -= -τ*u[end] 
        return SAT
    end
end
function SAT_Neumann(::NodeType{:Left},u::AbstractVector,Δx::Float64,RHS;
        c=1.0,order::Int64=2,forcing::Bool=false)
    τ = SATpenalties(Neumann,Δx,order)
    SAT = zeros(Float64,order)
    
    SAT[1] = τ*(c[1]*BDₓ(u,Left,Δx,order) - RHS[1])
    return SAT
end
function SAT_Neumann(::NodeType{:Right},u::AbstractVector,Δx::Float64,RHS;
        c=1.0,order::Int64=2,forcing::Bool=false)
    τ = SATpenalties(Neumann,Δx,order)
    SAT = zeros(Float64,order)

    SAT[end] = -τ*(c[end]*BDₓ(u,Right,Δx,order) - RHS[end])
    return SAT
end

"""
    SAT_Robin

1. SAT_Robin(::NodeType{:Left},u::Vector{Float64},Δx::Float64;
        order=2,a=1.0,b=1.0,forcing=false)
2. SAT_Robin(::NodeType{:Left},u::AbstractVector,Δx::Float64,RHS;
        order=2,a=1.0,b=1.0,forcing=false)


Simulatenous approximation term for Robin boundary conditions
`NodeType` is either `Left` or `Right`
"""
function SAT_Robin end
function SAT_Robin(::NodeType{:Left},u::Vector{Float64},Δx::Float64;
        order=2,a=1.0,b=1.0,forcing=false)
    # Get penalties
    τ = SATpenalties(Robin,a,Δx,order)
    SAT = zeros(Float64,order)

    if !forcing
        SAT[1] = b*u[1] + a*BDₓ(u,Left,Δx,order)
        return SAT
    else
        SAT[1] = -τ*u[1]
        return SAT
    end
end
function SAT_Robin(::NodeType{:Right},u::Vector{Float64},Δx::Float64;
        order=2,a=1.0,b=1.0,forcing=false)
    # Get penalties
    τ = SATpenalties(Robin,a,Δx,order)
    SAT = zeros(Float64,order)

    if !forcing
        SAT[end] = b*u[end] - a*BDₓ(u,Right,Δx,order)
        return SAT
    else
        SAT[end] = τ*u[end]
        return SAT
    end
end
function SAT_Robin(::NodeType{:Left},u::AbstractVector,Δx::Float64,RHS;
        order=2,a=1.0,b=1.0,forcing=false)
    # Get penalties
    τ = SATpenalties(Robin,a,Δx,order)
    SAT = zeros(Float64,order)

    SAT[1] = τ * (b*u[1] + a*BDₓ(u,Left,Δx,order) - RHS[1])
    return SAT
end
function SAT_Robin(::NodeType{:Right},u::AbstractVector,Δx::Float64,RHS;
        order=2,a=1.0,b=1.0,forcing=false)
    # Get penalties
    τ = SATpenalties(Robin,a,Δx,order)
    SAT = zeros(Float64,order)

    SAT[end] = τ * (b*u[end] - a*BDₓ(u,Right,Δx,order) + RHS[end])
    return SAT
end


"""
    SAT_Periodic(u::Vector{Float64},Δx::Float64,c::Vector{Float64};order::Int64=2)

Simulatenous approximation term for Periodic boundary conditions
"""
function SAT_Periodic(u::Vector{Float64},Δx::Float64,c::Vector{Float64};order::Int64=2)

    # Get h value
    h = hval(order)

    # Penalties
    α₀ = 0.5/(h * Δx) # Derivative penatly
    τ₁ = -0.5/(h * Δx) # Symmeteriser penalty
    τ₀ = -max(c[1]/2(h*Δx),c[end]/2(h*Δx))/(h*Δx) # Dirichlet penalty

    SAT = zeros(Float64,2order)

    # Dirichlet terms
    SAT[1] += τ₀ * (u[1] - u[end])
    SAT[end] += τ₀ * (u[end] - u[1])

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


"""
    Split_domain(u⁻::Vector{Float64},u⁺::Vector{Float64},Δx⁻::Float64,Δx⁺::Float64,c⁻,c⁺;
        order::Int64=2,order⁻::Int64=2,order⁺::Int64=2,separate_forcing::Bool=false)
Simulatenous approximation term for a split domain
    TODO: Fix for when `order⁺ != order⁻`
"""
function Split_domain(u⁻::Vector{Float64},u⁺::Vector{Float64},Δx⁻::Float64,Δx⁺::Float64,c⁻,c⁺;order::Int64=2,order⁻::Int64=2,order⁺::Int64=2,separate_forcing::Bool=false)

    h⁻ = hval(order⁻)
    h⁺ = hval(order⁺)
    
    # Left side
    α₀ = -0.5
    τ₁ = 0.5
    τ₀ = max(c⁻[end]/(2*h⁻*Δx⁻), c⁺[1]/(2*h⁺*Δx⁺))

    SAT₀ = zeros(Float64,order⁻+order⁺)
    SAT₁ = zeros(Float64,order⁻+order⁺)
    SAT₂ = zeros(Float64,order⁻+order⁺)
    
    F = zeros(Float64,order⁻+order⁺)

    # Function condition
    L₀u = zeros(Float64,order⁻+order⁺)
    L₀u[order⁻] += u⁻[end]
    L₀u[order⁻+1] += u⁻[end]

    SAT₀[1:order⁻] = τ₀/(h⁻ * Δx⁻) * L₀u[1:order⁻]
    SAT₀[order⁻+1:end] = τ₀/(h⁺ * Δx⁺) * L₀u[order⁻+1:end]

    F[1:order⁻] = -τ₀/(h⁻ * Δx⁻) * u⁺[1]
    F[order⁻+1:end] = τ₀/(h⁺ * Δx⁺) * u⁻[end]

    
    # Symmeteriser
    DₓᵀL₀u⁻ = boundary_Dₓᵀ(L₀u[1:order⁻],Δx⁻,order⁻)
    DₓᵀL₀u⁺ = boundary_Dₓᵀ(L₀u[order⁻+1:end],Δx⁺,order⁺)

    SAT₁[1:order⁻] = τ₁/(h⁻ * Δx⁻) * c⁻[end] * DₓᵀL₀u⁻[order⁻+1:end]
    SAT₁[order⁻+1:end] = -τ₁/(h⁺ * Δx⁺) * c⁺[1] * DₓᵀL₀u⁺[1:order⁺]

    # Derivative condition
    Dₓu⁻ = boundary_Dₓ(u⁻,Δx⁻,order⁻)
    Dₓu⁺ = boundary_Dₓ(u⁺,Δx⁺,order⁺)

    SAT₂[1:order⁻] = α₀/(h⁻ * Δx⁻) * c⁻[end] * Dₓu⁻[order⁻+1:end]
    SAT₂[order⁻+1:end] = α₀/(h⁺ * Δx⁺) * c⁺[1] * Dₓu⁺[1:order⁺]

    SAT = SAT₀ + SAT₁ + SAT₂
    
    SAT⁻ = SAT[1:order⁻]
    SAT⁺ = SAT[order⁻+1:end]


    if !separate_forcing
        SAT⁻ += F[1:order⁻]
        SAT⁺ += F[order⁻+1:end]
        return SAT⁻, SAT⁺
    else
        return SAT⁻, SAT⁺, F[1:order⁻], F[order⁻+1:end]
    end

end




"""
    BDₓᵀ

1. BDₓᵀ(u::AbstractVector,::NodeType{:Left},Δx::Float64,order::Int64=2)
2. BDₓᵀ(u::AbstractVector,::NodeType{:Right},Δx::Float64,order::Int64=2)

Transpose of first order derivative operator at the boundary needed to compute certain SATs
"""
function BDₓᵀ end
function BDₓᵀ(u::AbstractVector,::NodeType{:Left},Δx::Float64,order::Int64=2)
    if order == 2
        BOp = [-1.0, 1.0]
    elseif order == 4
        BOp = [-24.0/17.0, 59.0/34.0, -4.0/17.0, -3.0/34.0]
    elseif order == 6
        BOp = [-1.582533518939116, 2.033378678700676, -0.141512858744873, -0.450398306578272, 0.104488069284042, 0.036577936277544]
    else
        error("Order must be 2, 4, or 6")
    end
    BOp *= u[1]/Δx
    return BOp
end
function BDₓᵀ(u::AbstractVector,::NodeType{:Right},Δx::Float64,order::Int64=2)
    if order == 2
        BOp = [-1.0, 1.0]
    elseif order == 4
        BOp = [-24.0/17.0, 59.0/34.0, -4.0/17.0, -3.0/34.0]
    elseif order == 6
        BOp = [-1.582533518939116, 2.033378678700676, -0.141512858744873, -0.450398306578272, 0.104488069284042, 0.036577936277544]
    else
        error("Order must be 2, 4, or 6")
    end
    BOp *= u[end]/Δx
    return BOp
end






"""
    boundary_Dₓ(u::AbstractVector{Float64},Δx::Float64,order::Int64=2)
Transpose of first order derivative operator at the boundary needed to compute certain SATs
"""
function boundary_Dₓ(u::Vector{Float64},Δx::Float64,order::Int64=2)
    # Implementation of the 

    uₓ = zeros(Float64,2*order)

    if order == 2
        #Boundary for the second order case

        uₓ[1] = (u[2] - u[1])/Δx

        uₓ[end] = (u[end] - u[end-1])/Δx
        
        
    elseif order == 4
        #Boundary for the fourth order case

        B0p = [-24.0/17.0, 59.0/34.0, -4.0/17.0, -3.0/34.0]
        uₓ[1] = sum(B0p.*u[1:4])/Δx
        
        Bnp = [3.0/34.0, 4.0/17.0, -59.0/34.0, 24.0/17.0]
        uₓ[end] = sum(Bnp.*u[end-3:end])/Δx
        
        
    elseif order == 6
        #Boundary for the sixth order case

        B0p = [-1.582533518939116, 2.033378678700676, -0.141512858744873, -0.450398306578272, 0.104488069284042, 0.036577936277544]
        uₓ[1] = sum(B0p.*u[1:6])/Δx
        
        Bnp = [-0.036577936277544, -0.104488069284042, 0.450398306578272, 0.141512858744873, -2.033378678700676, 1.582533518939116]
        uₓ[end] = sum(Bnp.*u[end-5:end])/Δx

    end

    return uₓ
end



"""
    boundary_Dₓ(u::AbstractVector{Float64},Δx::Float64,order::Int64=2)

1. BDₓ(u::AbstractVector,::NodeType{:Left},Δx::Float64,order::Int64=2)
2. BDₓ(u::AbstractVector,::NodeType{:Right},Δx::Float64,order::Int64=2)

Transpose of first order derivative operator at the boundary needed to compute certain SATs
"""
function BDₓ end
function BDₓ(u::AbstractVector,::NodeType{:Left},Δx::Float64,order::Int64=2)
    if order == 2
        #Boundary for the second order case
        return (u[2] - u[1])/Δx
    elseif order == 4
        #Boundary for the fourth order case
        B0p = [-24.0/17.0, 59.0/34.0, -4.0/17.0, -3.0/34.0]
        return sum(B0p.*u[1:4])/Δx
    elseif order == 6
        #Boundary for the sixth order case
        B0p = [-1.582533518939116, 2.033378678700676, -0.141512858744873, -0.450398306578272, 0.104488069284042, 0.036577936277544]
        return sum(B0p.*u[1:6])/Δx
    end
    return uₓ
end
function BDₓ(u::AbstractVector,::NodeType{:Right},Δx::Float64,order::Int64=2)
    if order == 2
        #Boundary for the second order case
        return (u[end] - u[end-1])/Δx
    elseif order == 4
        #Boundary for the fourth order case
        Bnp = [3.0/34.0, 4.0/17.0, -59.0/34.0, 24.0/17.0]
        return sum(Bnp.*u[end-3:end])/Δx
    elseif order == 6
        #Boundary for the sixth order case
        Bnp = [-0.036577936277544, -0.104488069284042, 0.450398306578272, 0.141512858744873, -2.033378678700676, 1.582533518939116]
        return sum(Bnp.*u[end-5:end])/Δx
    end
end





"""
    SATpenalties

1. `SATpenalties(::BoundaryCondition{:Dirichlet,:Neumann},Δx::Float64,order::Int64)`
2. `SATpenalties(::BoundaryCondition{:Robin},a,Δx::Float64,order::Int64)`

Determines the penatly parameters for the given boundary conditions.
"""
function SATpenalties end
function SATpenalties(::BoundaryCondition{:Dirichlet},Δx::Float64,order::Int64)
    # For reading in penalty parameters for Dirichlet SATs
    h = hval(order)

    α = 1.0 * (h * Δx)^-1 # α*H^{-1}

    τ = 1.0
    τ = -(1.0 + τ) * (h * Δx)^-2 # τ*H^{-1}H^{-1}
    return α, τ
end
function SATpenalties(::BoundaryCondition{:Neumann},Δx::Float64,order::Int64)
    # For reading in penalty parameters for Neumann SATs
    h = hval(order)

    τ = 1.0/(h * Δx) # τ*H^{-1}
    return τ
end
function SATpenalties(::BoundaryCondition{:Robin},a,Δx::Float64,order::Int64)
    h = hval(order)

    τ = 1.0/(a * h * Δx) # τ=1/a H^{-1}
    return τ
end


"""
    hval(order::Int64)
Returns the value of ``h^{-1}`` for the penalties
"""
function hval(order::Int64)
    if order == 2
        return 0.5
    elseif order == 4
        return 17.0/48.0
    elseif order == 6
        return 13649.0/43200.0
    end
end

















#=
====================== Dirichlet boundary conditions ======================
=#
function SAT_Dirichlet1 end
function SAT_Dirichlet1(::NodeType{:Left},u::Vector{Float64},Δx::Float64,g;
        c=1.0,order::Int64=2,separate_forcing::Bool=false)
    # Penalty parameters
    α,τ = SATpenalties(Dirichlet,Δx,order)
    
    # Construct the SATs
    SAT = zeros(Float64,order)
    Dᵀu = boundary_Dₓᵀ(u,Δx,order)
    SAT[1] += τ * u[1] 
    SAT .+= α * c[1] * Dᵀu[1:order]
    
    # Forcing terms
    F = zeros(Float64,order)
    Dᵀf = boundary_Dₓᵀ(g,Δx,order)
    F[1] += -τ*g[1]
    F .+= -α * c[1] * Dᵀf[1:order]
    
    return SAT, F
end
function SAT_Dirichlet1(::NodeType{:Right},u::Vector{Float64},Δx,g;
    c=1.0,order::Int64=2,separate_forcing::Bool=false)
    # Penalty parameters
    α,τ = SATpenalties(Dirichlet,Δx,order)
    
    # Construct the SATs
    SAT = zeros(Float64,order)
    Dᵀu = boundary_Dₓᵀ(u,Δx,order)
    SAT[end] += τ*u[end]
    SAT .+= α * c[end] * Dᵀu[end-order+1:end]
    
    # Forcing terms
    F = zeros(Float64,order)
    Dᵀf = boundary_Dₓᵀ(g,Δx,order)
    F[end] += -τ*g[end]
    F .+= -α * c[end] * Dᵀf[end-order+1:end]
    
    return SAT, F
end