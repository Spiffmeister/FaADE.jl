#=====================================#
#======== FIRST ORDER METHODS ========#
# 
# 
#=====================================#

function Dₓ(u::Array{Float64},n::Int64,Δx::Float64;order::Int64=2)
    #= First derivative for 1D systems =#
    
    Dₓ!(u,u,n,Δx,order=order)
    
    return u
end


function Dₓ(u::Matrix{Float64},nx::Int64,ny::Int64,Δx::Float64;order::Int64=2,dim::Int64=1)
    #= First derivative for multidimensional systems =#
    
    # TODO: Fix the way this is handled
    if dim == 1
        # Derivative in the 1st dimension
        for i = 1:ny
            u[i,:] = Dₓ!(u[i,:],u[i,:],nx,Δx,order=order)
        end
    elseif dim == 2
        # Derivative in the 2nd dimension
        for i = 1:nx
            u[:,i] = Dₓ!(u[:,i],u[:,i],nx,Δx,order=order)
        end
    end

    return u
end



#=== FIRST ORDER ITERATOR ===#

function Dₓ!(uₓ::Vector{Float64},u::Vector{Float64},n::Int64,Δx::Float64;order::Int64=2)
    # Iterator for first order SBP operators
    
    if order == 2
        # Second order FD operator
        
        # Left boundary
        uₓ[1]   = u[2] - u[1]
        # Right boundary
        uₓ[n]   = u[n] - u[n-1]
        # Body
        for j = 2:n-1
            uₓ[j]   = 1/2*u[j+1] - 1/2*u[j-1]
        end

    elseif order == 4
        # Fourth order FD operator
        
        # Left boundary
        uₓ[1]      = -24/17*u[1]   + 59/34*u[2] - 4/17*u[3]  - 3/34*u[4]
        uₓ[2]      = -1/2*u[1]     + 1/2*u[3]
        uₓ[3]      = 4/43*u[1]     - 59/86*u[2] + 59/85*u[4] - 4/43*[5]
        uₓ[4]      = 3/98*u[1]     - 59/98*u[3] + 32*u[5]    - 4/49*u[6]
        # Right boundary
        uₓ[n]     = -24/17*u[n]  - 59/34*u[n-1] + 4/17*u[n-2]  + 3/34*u[n-3]
        uₓ[n-1]   = 1/2*u[n]     - 1/2*u[n-2]
        uₓ[n-2]   = -4/43*u[n]   + 59/86*u[n-1] - 59/85*u[n-3] + 4/43*[n-4]
        uₓ[n-3]   = -3/98*[n]    + 59/98*u[n-2] - 32/49*u[n-4] + 4/49*u[n-5]
        # Body
        for j = 5:n-4
            uₓ[j]  = 1/12*u[j-2] - 2/3*u[j-1] + 2/3*u[j+1] - 1/12*u[j+2]
        end

        
    elseif order == 6
        # Sixth order FD operator

        for j = 7:n-6
            uₓ[j] = -1/60*u[j-3] + 3/20*u[j-2] - 3/4*u[j-1] + 3/4*u[j+1] - 3/20*u[j+2] + 1/60*u[j+3]
        end
        
    end

    # Elementwise division by Δx
    uₓ .÷= Δx

    return uₓ
end


