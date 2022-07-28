#=
    Second derivative boundary operators
=#


"""
    SAT_left(type::Symbol,u::Vector{Float64},Δx::Float64,g;order=2::Int,c::Union{Float64,Vector{Float64}}=1.0,αβ=[1.0,1.0],separate_forcing=false)

Returns the left SAT+Forcing term as a vector or SAT and Forcing term as separate vectors

Inputs: 
- type: :Dirichlet, :Neumann, :Robin
- u: Solution vector
- Δx: step size
- g: boundary value
Optional inputs:
- order: 2 (default), 4, 6
- coefficient: c in ∂ₓ(c∂ₓu)
    c=1.0 (default) integer for
- αβ: [α,β] Robin boundary coefficient αdₓu + βu = f
- separate_forcing: true, false (default) -- controls the return type.
    false; returns a single vector of the SAT
    true; returns SAT and F(orcing) vectors individually

For periodic conditions call Periodic(u,Δx,c;order), for a split domain call Split_domain(u⁻,u⁺,Δx⁻,Δx⁺,c⁻,c⁺;order=2,order⁻=2,order⁺=2)
"""



function SAT_left(type::BoundaryCondition,u::Vector{Float64},Δx::Float64,g;
    order=2::Int,c::Union{Float64,Vector{Float64}}=1.0,αβ::Vector{Float64}=[1.0,1.0],separate_forcing::Bool=false)
    if type == Dirichlet # Dirichlet boundaries
        # Ensure coefficient is correct type
        SAT, F = Dirichlet_left(u,Δx,g,c=c,order=order)

    elseif type == Neumann # Neumann boundary conditions
        # Ensure coefficient is correct type
        if typeof(c) <: Vector
            k = c[1]
        else
            k = c
        end
        SAT,F = Neumann_left(u,Δx,g,c=k,order=order)
        
    elseif type == Robin # Robin boundaries
        # Ensure coefficient is correct type
        SAT,F = Robin_left(u,Δx,g,order=order,a=αβ[1],b=αβ[2])
        
    else
        error("Must specify either :Dirichlet, :Neumann or :Robin. For periodic boundaries use the SAT_Periodic() function")
    end

    if !separate_forcing # If the forcing term is part of the SAT
        SAT += F
        return SAT
    else
        return SAT, F
    end
end


function SAT_left(type::Symbol,u::Matrix{Float64},Δx::Float64,Δy::Float64,nx::Int64,ny::Int64,g,dim;
    order::Int64=2,c::Union{Float64,Vector{Float64}}=1.0,αβ::Vector{Float64}=[1.0,1.0],separate_forcing::Bool=false)
    # Matrix variant of SAT_left
    if dim == 1
        SAT = zeros(Float64,ny,order)
        F = zeros(Float64,ny,order)
        for i = 1:ny
            SAT[i,:],F[i,:] = SAT_left(type=type,u[i,:],Δx,g,order=order,c=c,αβ=αβ,separate_forcing=separate_forcing)
        end
    elseif dim == 2
        SAT = zeros(Float64,nx,order)
        F = zeros(Float64,nx,order)
        for i = 1:nx
            SAT[i,:],F[i,:] = SAT_left(type=type,u[i,:],Δy,g,order=order,c=c,αβ=αβ,separate_forcing=separate_forcing)
        end
    end
    SAT += F
    return SAT
end



function SAT_right(type::Symbol,u::Vector{Float64},Δx::Float64,g;order::Int=2,c::Union{Float64,Vector{Float64}}=1.0,a::Float64=1.0,b::Float64=1.0,separate_forcing=false)

    # Right side SAT
    if type == :Dirichlet
        if typeof(c) <: Float64
            c = zeros(Float64,order) .+ c
        end
        # SAT, F = Dirichlet_right(u,Δx,g,c=c,order=order)
        SAT, F = Dirichlet_right(u,Δx,g,c=c,order=order)

    elseif type == :Neumann # Neumann boundary conditions

        # Ensure coefficient is correct type
        if typeof(c) <: Vector
            k = c[end]
        else
            k = c
        end

        SAT, F = Neumann_right(u,Δx,g,c=c,order=order)
        
    elseif type == :Robin

        SAT, F = Robin_right(u,Δx,g,order=order,a=a,b=b)

    else
        error("Must specify either :Dirichlet, :Neumann or :Robin. For :Periodic use the SAT_Periodic() function")
    end

    if !separate_forcing
        # If the forcing term is part of the SAT (i.e. not using GC method)
        SAT += F
        return SAT
    else
        return SAT, F
    end
end




#=
====================== Dirichlet boundary conditions ======================
=#
# function SAT_left end
# function SAT_left(::BoundaryCondition{:Dirichlet},u::AbstractVector,Δx::Float64,g;
#         c=1.0,order::Int64=2,penalty::Vector{Float64}=[-1.0,-1.0],separate_forcing::Bool=false)
#     α,τ = Dirichlet_penalties(Δx,order)

#     # Construct the SATs
#     SAT = zeros(Float64,order)
#     Dᵀu = boundary_Dₓᵀ(u,Δx,order)
#     SAT[1] += τ * u[1] 
#     SAT .+= α * c[1] * Dᵀu[1:order]

    
#     # Forcing terms
#     F = zeros(Float64,order)
#     Dᵀf = boundary_Dₓᵀ(g,Δx,order)
#     F[1] += -τ*g[1]
#     F .+= -α * c[1] * Dᵀf[1:order]
    
#     return SAT, F
# end

"""
    Dirichlet_left(u::Vector{Float64},Δx::Float64,g;c=1.0,order::Int64=2,penalty::Vector{Float64}=[-1.0,-1.0])
"""


function Dirichlet_left(u::Vector{Float64},Δx::Float64,g;
    c=1.0,order::Int64=2,seperate_forcing::Bool=false)

    # Penalty parameters
    α,τ = Dirichlet_penalties(Δx,order)
    
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
    
    if seperate_forcing
        return SAT, F
    else
        return SAT + F
    end
end

function Dirichlet_right(u::Vector{Float64},Δx,g;
    c=1.0,order::Int64=2,seperate_forcing::Bool=false)

    # Penalty parameters
    α,τ = Dirichlet_penalties(Δx,order)
    
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

    if seperate_forcing
        return SAT, F
    else
        return SAT + F
    end
end


function Dirichlet_penalties(Δx::Float64,order::Int64)
    # For reading in penalty parameters for Dirichlet SATs
    h = hval(order)

    α = 1.0 * (h * Δx)^-1

    τ = 1.0
    τ = -(1.0 + τ) * (h * Δx)^-2

    return α, τ
end



#=
====================== Neumann boundary conditions ======================
=#


function Neumann_left(u::Vector{Float64},Δx::Float64,g;
    c::Float64=1.0,order::Int64=2,seperate_forcing::Bool=false)

    # Penalties
    τ = Neumann_penalties(Δx,order)
    # SAT construction
    SAT = zeros(Float64,order)
    Du = boundary_Dₓ(u,Δx,order)
    SAT[1] += τ * c * Du[1]
    # Forcing term
    F = zeros(Float64,order)
    F[1] += -τ*g[1]

    if seperate_forcing
        return SAT, F
    else
        return SAT + F
    end
end

function Neumann_right(u::Vector{Float64},Δx::Float64,g;c=1.0,order=2,penalty=1.0)
    # Penalties
    τ = Neumann_penalties(Δx,order)
    # SAT construction
    SAT = zeros(Float64,order)
    Du = boundary_Dₓ(u,Δx,order)
    SAT[end] -= τ * c[end] * Du[end]
    # Forcing term
    F = zeros(Float64,order)
    F[end] -= -τ*g[end]

    return SAT, F
end


function Neumann_penalties(Δx::Float64,order::Int64;penalties::Float64=-1.0)
    # For reading in penalty parameters for Neumann SATs
    h = hval(order)
    τ = 1.0/(h * Δx)

    return τ
end

#=
====================== Robin boundary conditions ======================
=#

"""
    function Robin_left(u::Vector{Float64},Δx::Float64,g;order=2,a=1.0,b=1.0)

"""
function Robin_left(u::Vector{Float64},Δx::Float64,g;order=2,a=1.0,b=1.0)
    # Get penalties
    h = hval(order)
    τ = 1.0/(a * h * Δx)
    
    # Compute the SAT
    SAT = zeros(Float64,order)
    Du = boundary_Dₓ(u,Δx,order)
    SAT[1] = τ * b*u[1] + a*Du[1]
    
    # Forcing terms
    F = zeros(Float64,order)
    F[1] = -τ*g[1]
    
    return SAT, F
end

function Robin_right(u::Vector{Float64},Δx::Float64,g;order=2,a=1.0,b=1.0)
    # Get penalties
    h = hval(order)
    τ = 1.0/(a * h * Δx)

    # Compute the SAT
    SAT = zeros(Float64,order)
    Du = boundary_Dₓ(u,Δx,order)
    SAT[end] = τ * b*u[end] - a*Du[end] #-Du for directional derivative
    
    # Forcing terms
    F = zeros(Float64,order)
    F[end] = -τ * g[end]

    SAT, F
end


#=
====================== Periodic boundary conditions ======================
=#
function SAT_Periodic(u::Vector{Float64},Δx::Float64,c::Vector{Float64};order::Int64=2,separate_forcing::Bool=false)

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
    L₁u = zeros(Float64,2order)
    
    L₁u[1] = (u[1] - u[end])
    L₁u[end] = (u[1] - u[end])
    
    Dᵀu = boundary_Dₓᵀ(L₁u,Δx,order)

    SAT[1:order] += τ₁ * c[1] * Dᵀu[1:order]
    SAT[order+1:end] += -τ₁ * c[end] * Dᵀu[end-order+1:end] # negative for directional derivative

    # Neumann term
    Du = boundary_Dₓ(u,Δx,order)
    SAT[1] += α₀ * (c[1]*Du[1] - c[end]*Du[end]) # corrected for directional derivative
    SAT[end] += α₀ * (c[1]*Du[1] - c[end]*Du[end]) # corrected directional derivative

    return SAT[1:order], SAT[order+1:end]

end



#=
====================== Interface boundary ======================
=#
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




#=
====================== Supporting fns ======================
=#
function boundary_Dₓᵀ(u::Vector{Float64},Δx::Float64,order::Int64=2)
    # Implementation of the 

    uₓ = zeros(Float64,2*order)

    if order == 2
        #Boundary for the second order case

        BOp = [-1.0, 1.0]

        uₓ[1:2] = BOp*u[1]/Δx   # Left boundary
        uₓ[end-1:end] = BOp*u[end]/Δx   # Right boundary
        
    elseif order == 4
        #Boundary for the fourth order case

        BOp = [-24.0/17.0, 59.0/34.0, -4.0/17.0, -3.0/34.0]
        
        uₓ[1:4] = BOp*u[1]/Δx   # Left boundary
        uₓ[end-3:end] = -BOp*u[end]/Δx  # Right boundary
        
        
    elseif order == 6
        #Boundary for the sixth order case

        BOp = [-1.582533518939116, 2.033378678700676, -0.141512858744873, -0.450398306578272, 0.104488069284042, 0.036577936277544]
        
        uₓ[1:6] = BOp*u[1]/Δx   # Left boundary
        uₓ[end-5:end] = -BOp*u[end]/Δx  # Right boundary

    end

    return uₓ
end



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



function hval(order::Int64)
    if order == 2
        h = 0.5
    elseif order == 4
        h = 17.0/48.0
    elseif order == 6
        #TODO
        h = 13649.0/43200.0
    end
    return h
end