#=
    Second derivative boundary operators
=#



function SAT_left(type::Symbol,u::Vector{Float64},Δx::Float64,g;order=2,c::Union{Float64,Vector{Float64}}=1.0,a=1.0,b=1.0)
    # Left side SAT

    if type == :Dirichlet
        SAT = Dirichlet_left(u,Δx,g,c=c,order=order)
    elseif type == :Neumann
        SAT = Neumann_left(u,Δx,g,c=c,order=order)
    elseif type == :Robin
        SAT = Robin_left(u,Δx,g,order=order,a=a,b=b)
    else
        error("Must specify either :Dirichlet, :Neumann or :Robin. For :Periodic use the Periodic() function")
    end

    return SAT
end


function SAT_right(type::Symbol,u::Vector{Float64},Δx::Float64,g;order=2,c::Union{Float64,Vector{Float64}}=1.0,a=1.0,b=1.0)
    # Right side SAT
    if type == :Dirichlet
        SAT = Dirichlet_right(u,Δx,g,c=c,order=order)
    elseif type == :Neumann
        SAT = Neumann_right(u,Δx,g,c=c,order=order)
    elseif type == :Robin
        SAT = Robin_right(u,Δx,g,order=order,a=a,b=b)
    else
        error("Must specify either :Dirichlet, :Neumann or :Robin. For :Periodic use the SAT() function")
    end

    return SAT
end




#=
====================== Dirichletee boundary conditions ======================
=#


function Dirichlet_left(u::Vector{Float64},Δx,g;c=1.0,order=2,penalty=[-1.0,-1.0])

    # Fix up the coefficient 
    if typeof(c) <: Float64
        c = zeros(Float64,order) .+ c
    end
    
    # Penalty parameters
    α,τ = Dirichlet_penalties(Δx,order)

    # Construct the SATs
    SAT = zeros(Float64,order)

    Dᵀu = boundary_Dₓᵀ(u,Δx,order)
    Dᵀf = boundary_Dₓᵀ(g,Δx,order)

    SAT[1] += τ * (u[1] - g[1]) 
    SAT .+= α * c[1] * (Dᵀu[1:order] - Dᵀf[1:order])

    return SAT
end


function Dirichlet_right(u::Vector{Float64},Δx,g;c=1.0,order=2,penalty=[-1.0,-1.0])
    
    # Fix up the coefficient 
    if typeof(c) <: Float64
        c = zeros(Float64,order) .+ c
    end
    
    # Penalty parameters
    α,τ = Dirichlet_penalties(Δx,order)
    
    # Construct the SATs
    SAT = zeros(Float64,order)
    
    Dᵀu = boundary_Dₓᵀ(u,Δx,order)
    Dᵀf = boundary_Dₓᵀ(g,Δx,order)
    
    # Right SAT
    SAT[end] += τ*(u[end] - g[end])
    SAT .+= α * c[end] * (Dᵀu[end-order+1:end] - Dᵀf[end-order+1:end])
    return SAT
end


function Dirichlet_penalties(Δx,order;penalty=[-1.0,-1.0])
    # For reading in penalty parameters for Dirichlet SATs

    h = hval(order)

    α = penalty[1]
    τ = penalty[2]

    if α == -1.0
        α = 1.0 * (h * Δx)^-1
    end

    if τ == -1.0
        τ = 1.0
        τ = -(1.0 + τ) * (h * Δx)^-2
    end

    return α, τ
end



#=
====================== Neumann boundary conditions ======================
=#


function Neumann_left(u::Vector{Float64},Δx::Float64,g;c=1.0,order=2,penalty=1.0)

    τ = Neumann_penalties(Δx,order)

    SAT = zeros(Float64,order)
    # If k is a scalar, vectorise it
    if typeof(c) <: Float64
        c = zeros(2*order) .+ c
    end
    
    Du = boundary_Dₓ(u,Δx,order)
    SAT[1] += τ * c[1] * Du[1]

    return SAT
end

function Neumann_right(u::Vector{Float64},Δx::Float64,g;c=1.0,order=2,penalty=1.0)

    τ = Neumann_penalties(Δx,order)

    SAT = zeros(Float64,order)

    # If k is a scalar, vectorise it
    if typeof(c) <: Float64
        c = zeros(2*order) .+ c
    end

    Du = boundary_Dₓ(u,Δx,order)
    SAT[end] -= τ * c[end] * Du[end]

    return SAT
end


function Neumann_penalties(Δx,order;penalties=-1.0)

    h = hval(order)

    τ = 1.0/(h * Δx)

    return τ
end

#=
====================== Robin boundary conditions ======================
=#

function Robin_left(u::Vector{Float64},Δx::Float64,g;order=2,a=1.0,b=1.0)

    # Get penalties
    h = hval(order)
    if typeof(a) <: Vector
        a = a[1]
    end
    if typeof(b) <: Vector
        b = b[1]
    end
    τ = 1.0/(a * h * Δx)

    # Preallocate SAT
    SAT = zeros(Float64,order)

    # Compute the SAT
    Du = boundary_Dₓ(u,Δx,order)
    SAT[1] = τ * (b*u[1] + a*Du[1] - g[1])

    return SAT
end

function Robin_right(u::Vector{Float64},Δx::Float64,g;order=2,a=1.0,b=1.0)

    # Get penalties
    h = hval(order)
    if typeof(a) <: Vector
        a = a[end]
    end
    if typeof(b) <: Vector
        b = b[end]
    end
    τ = 1.0/(a * h * Δx)

    SAT = zeros(Float64,order)

    Du = boundary_Dₓ(u,Δx,order)
    SAT[end] = τ * (b*u[end] - a*Du[end] - g[end]) #-Du for directional derivative

    SAT
end


#=
====================== Periodic boundary conditions ======================
=#


function Periodic(u::Vector{Float64},Δx::Float64,c::Vector{Float64};order=2)

    # Get h value
    h = hval(order)

    # Penalties
    α₀ = 0.5/(h * Δx) # Derivative penatly
    τ₁ = -0.5/(h * Δx) # Symmeteriser penalty

    τ₀ = -max(c[1]/2(h*Δx),c[end]/2(h*Δx))/(h * Δx) # Function penalty


    SAT₀ = zeros(Float64,2order)
    SAT₁ = zeros(Float64,2order)
    SAT₂ = zeros(Float64,2order)


    SAT₀[1] = τ₀ * (u[1] - u[end])
    SAT₀[end] = τ₀ * (u[end] - u[1])


    L₁u = zeros(Float64,2order)
    L₁u[1] = (u[1] - u[end])
    L₁u[end] = (u[1] - u[end])
    Dᵀu = boundary_Dₓᵀ(L₁u,Δx,order)

    SAT₁[1:order] = τ₁ * c[1] * Dᵀu[1:order]
    SAT₁[order+1:end] = τ₁ * c[end] * Dᵀu[end-order+1:end] # negative for directional derivative

    Du = boundary_Dₓ(u,Δx,order)
    SAT₂[1] = α₀ * (c[1]*Du[1] - c[end]*Du[end]) # corrected for directional derivative
    SAT₂[end] = α₀ * (c[1]*Du[1] - c[end]*Du[end]) # corrected directional derivative


    SAT = SAT₀ + SAT₁ + SAT₂

    return SAT[1:order], SAT[order+1:end]
end



#=
====================== Interface boundary ======================
=#


function Split_domain(u⁻,u⁺,Δx⁻,Δx⁺,c⁻,c⁺;order=2,order_left=2,order_right=2)

    h⁻ = hval(order_left)
    h⁺ = hval(order_right)
    
    # Left side
    α₀ = -0.5
    τ₁ = 0.5
    τ₀ = max(c⁻[end]/(2*h⁻*Δx⁻), c⁺[1]/(2*h⁺*Δx⁺))

    SAT₀ = zeros(Float64,2order)
    SAT₁ = zeros(Float64,2order)
    SAT₂ = zeros(Float64,2order)
    
    # Function condition
    L₀u = zeros(Float64,2order)
    L₀u[order] += (u⁻[end] - u⁺[1])
    L₀u[order+1] += (u⁻[end] - u⁺[1])

    SAT₀[1:order] = τ₀/(h⁻ * Δx⁻) * L₀u[1:order]
    SAT₀[1:order] = τ₀/(h⁺ * Δx⁺) * L₀u[order+1:end]

    
    # Symmeteriser
    DₓᵀL₀u = zeros(Float64,2order)

    DₓᵀL₀u⁻ = boundary_Dₓᵀ(L₀u[1:order],Δx⁻,order)
    DₓᵀL₀u[1:order] = DₓᵀL₀u⁻[order+1:end]
    DₓᵀL₀u⁺ = boundary_Dₓᵀ(L₀u[order+1:end],Δx⁺,order)
    DₓᵀL₀u[order+1:end] = DₓᵀL₀u⁺[1:order]

    SAT₁[1:order] = τ₁/(h⁻ * Δx⁻) * c⁻[end] * DₓᵀL₀u[1:order]
    SAT₁[order+1:end] = τ₁/(h⁺ * Δx⁺) * c⁺[1] * DₓᵀL₀u[order+1:end]

    # Derivative condition
    Dₓu⁻ = boundary_Dₓ(u⁻,Δx⁻,order)
    Dₓu⁺ = boundary_Dₓ(u⁺,Δx⁺,order)

    SAT₂[1:order] = α₀/(h⁻ * Δx⁻) * c⁻[end] * Dₓu⁻[order+1:end]
    SAT₂[order+1:end] = -α₀/(h⁺ * Δx⁺) * c⁺[1] * Dₓu⁺[1:order]

    # println("Du+: ",Dₓu⁺[1:order])
    # println("Du+/Deltax: ",c⁺[1] .* Dₓu⁺[1:order] ./ Δx⁺)
    # println("SAT_0: ",SAT₀[order+1:end])
    # println("SAT_1: ",SAT₁[order+1:end])
    # println("SAT_2: ",SAT₂[order+1:end])

    SAT = SAT₀ + SAT₁ + SAT₂
    
    SAT⁻ = SAT[1:order]
    SAT⁺ = SAT[order+1:end]



    return SAT⁻, SAT⁺

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