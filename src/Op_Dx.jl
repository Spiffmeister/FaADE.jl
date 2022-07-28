#=====================================#
#====== FIRST DERIVATIVE METHODS =====#
# 
# 
# TODO: Add 3D methods
#=====================================#



function Dₓ(u::Vector{Float64},n::Int64,Δx::Float64;order::Int64=2)
    #= First derivative for 1D systems =#
    # Computes the first derivative for 1D systems, the default order is 2
    uₓ = zeros(Float64,n)
    Dₓ!(uₓ,u,n,Δx,order=order)
    
    return uₓ
end


function Dₓ(u::Matrix{Float64},n::Vector{Int64},Δx::Float64;order::Int64=2,dims::Union{Int64,Vector{Int64}}=1)
    #= First derivative for multidimensional systems =#
    # TODO: Fix the way this is handled

    # try axes(u)[3]
    # catch
    # end

    uₓ = zeros(Float64,n)
    
    if 1 ∈ dims
        # Derivative in the 1st dimension
        # Dₓ_m!(uₓ,u) = Dₓ!(uₓ,u,nx,Δx,order=order)
        for i = 1:n[2]
            # u[i,:] = Dₓ_m!(u[i,:],u[i,:])
            Dₓ!(uₓ[:,i],u[:,i],nx,Δx,order=order)
        end
    end

    if 2 ∈ dims
        # Derivative in the 2nd dimension
        # Dₓ_m!(uₓ,u) = Dₓ!(uₓ,u,ny,Δx,order=order)
        for i = 1:n[1]
            # u[:,i] = Dₓ_m!(u[:,i],u[:,i])
            Dₓ!(uₓ[:,i],u[:,i],ny,Δx,order=order)
        end
    end

    return u
end



#=== FIRST ORDER ITERATOR ===#

function FirstDerivative end
### Internal Node
function FirstDerivative(u::AbstractVector{Float64},Δx::Float64,n::Int64,::NodeType{:Internal};order::Int64=2)

    if order == 2
        j = 1
        uₓ = (u[j+1] - u[j-1])/2.0
        return uₓ/Δx
    elseif order == 4
        j = 3
        uₓ = 1/12*u[j-2] - 2/3*u[j-1]  + 2/3*u[j+1] - 1/12*u[j+2]
        return uₓ/Δx
    elseif order == 6
        j = 4
        uₓ = -1/60*u[j-3] + 3/20*u[j-2] - 3/4*u[j-1] + 3/4*u[j+1] - 3/20*u[j+2] + 1/60*u[j+3]
        return uₓ/Δx
    end
end
### Left boundary
function FirstDerivative(u::AbstractVector{Float64},Δx::Float64,n::Int64,::NodeType{:Left};order::Int64=2)
    if order == 2
        return (u[2] - u[1])/Δx
    elseif order == 4
        uₓ = zeros(Float64,4)
        uₓ[1] = -24/17*u[1] + 59/34*u[2]  - 4/17*u[3] - 3/34*u[4]
        uₓ[2] = -1/2*u[1] + 1/2*u[3]
        uₓ[3] = 4/43*u[1] - 59/86*u[2]  + 59/86*u[4] - 4/43*u[5]
        uₓ[4] = 3/98*u[1] - 59/98*u[3]  + 32/49*u[5] - 4/49*u[6]
        return uₓ./Δx
    elseif order == 6
        uₓ = zeros(Float64,6)
        uₓ[1]   = -1.582533518939116*u[1]   + 2.033378678700676*u[2] - 0.141512858744873*u[3] +
            -0.450398306578272*u[4] + 0.104488069284042*u[5] + 0.036577936277544*u[6]
        uₓ[2]   = -0.462059195631158*u[1] + 0.287258622978251*u[3] + 0.258816087376832*u[4] +
            -0.069112065532624*u[5] - 0.014903449191300*u[6]
        uₓ[3]   = 0.071247104721830*u[1] - 0.636451095137907*u[2] + 0.606235523609147*u[4] +
            -0.022902190275815*u[5] - 0.018129342917256*u[6]
        uₓ[4]   = 0.114713313798970*u[1] - 0.290087484386815*u[2] - 0.306681191361148*u[3] +
            0.520262285050482*u[5] - 0.051642265516119*u[6] + 0.013435342414630*u[7]
        uₓ[5]   = -0.036210680656541*u[1] + 0.105400944933782*u[2] + 0.015764336127392*u[3] +
            -0.707905442575989*u[4] + 0.769199413962647*u[6] - 0.164529643265203*u[7] + 0.018281071473911*u[8]
        uₓ[6]   = -0.011398193015050*u[1] + 0.020437334208704*u[2] + 0.011220896474665*u[3] + 
            0.063183694641876*u[4] - 0.691649024426814*u[5] + 0.739709139060752*u[7] +
            - 0.147941827812150*u[8] + 0.016437980868017*u[9]
        return uₓ./Δx
    end
end
### Right boundary
function FirstDerivative(u::AbstractVector{Float64},Δx::Float64,n::Int64,::NodeType{:Right};order::Int64=2)
    if order == 2
        return (u[2] - u[1])/Δx
    elseif order == 4
        m = 4
        n = 6
        uₓ = zeros(Float64,4)
        uₓ[m] = 24/17*u[n] - 59/34*u[n-1]  + 4/17*u[n-2] + 3/34*u[n-3]
        uₓ[m-1] = 1/2*u[n] - 1/2*u[n-2]
        uₓ[m-2] = -4/43*u[n] + 59/86*u[n-1]  - 59/86*u[n-3] + 4/43*u[n-4]
        uₓ[m-3] = -3/98*u[n] + 59/98*u[n-2]  - 32/49*u[n-4] + 4/49*u[n-5]
        return uₓ./Δx
    elseif order == 6
        m = 6
        n = 9
        uₓ = zeros(Float64,6)
        uₓ[m]      = 1.582533518939116*u[n]   - 2.033378678700676*u[n-1] + 0.141512858744873*u[n-2] +
            0.450398306578272*u[n-3] - 0.104488069284042*u[n-4] - 0.036577936277544*u[n-5]
        uₓ[m-1]    = 0.462059195631158*u[n] - 0.287258622978251*u[n-2] - 0.258816087376832*u[n-3] +
            0.069112065532624*u[n-4] + 0.014903449191300*u[n-5]
        uₓ[m-2]    = -0.071247104721830*u[n] + 0.636451095137907*u[n-1] - 0.606235523609147*u[n-3] +
            0.022902190275815*u[n-4] + 0.018129342917256*u[n-5]
        uₓ[m-3]    = -0.114713313798970*u[n] + 0.290087484386815*u[n-1] + 0.306681191361148*u[n-2] +
            -0.520262285050482*u[n-4] + 0.051642265516119*u[n-5] - 0.013435342414630*u[n-6]
        uₓ[m-4]    = 0.036210680656541*u[n] - 0.105400944933782*u[n-1] + 0.015764336127392*u[n-2] +
            0.707905442575989*u[n-3] - 0.769199413962647*u[n-5] + 0.164529643265203*u[n-6] - 0.018281071473911*u[n-7]
        uₓ[m-5]    = 0.011398193015050*u[n] - 0.020437334208704*u[n-1] - 0.011220896474665*u[n-2] + 
            -0.063183694641876*u[n-3] + 0.691649024426814*u[n-4] - 0.739709139060752*u[n-6] +
            0.147941827812150*u[n-7] - 0.016437980868017*u[n-8]
        return uₓ./Δx
    end
end
### OLD
function Dₓ!(uₓ::Vector{Float64},u::Vector{Float64},n::Int64,Δx::Float64;order::Int64=2)
    # Iterator for first derivative SBP operators
    
    if order == 2
        # Second order FD operator
        
        # Left boundary
        uₓ[1]   = u[2] - u[1]
        # Right boundary
        uₓ[n]   = u[n] - u[n-1]
        # Body
        for j = 2:n-1
            uₓ[j]   = 1/2*(u[j+1] - u[j-1])
        end

    elseif order == 4
        # Fourth order FD operator
       
        uₓ[1] = -24/17*u[1] + 59/34*u[2]  - 4/17*u[3] - 3/34*u[4]
        uₓ[2] = -1/2*u[1] + 1/2*u[3]
        uₓ[3] = 4/43*u[1] - 59/86*u[2]  + 59/86*u[4] - 4/43*u[5]
        uₓ[4] = 3/98*u[1] - 59/98*u[3]  + 32/49*u[5] - 4/49*u[6]

        uₓ[n] = 24/17*u[n] - 59/34*u[n-1]  + 4/17*u[n-2] + 3/34*u[n-3]
        uₓ[n-1] = 1/2*u[n] - 1/2*u[n-2]
        uₓ[n-2] = -4/43*u[n] + 59/86*u[n-1]  - 59/86*u[n-3] + 4/43*u[n-4]
        uₓ[n-3] = -3/98*u[n] + 59/98*u[n-2]  - 32/49*u[n-4] + 4/49*u[n-5]

        for j = 5:n-4
        
            uₓ[j] = 1/12*u[j-2] - 2/3*u[j-1]  + 2/3*u[j+1] - 1/12*u[j+2];
            
        end

    
    elseif order == 6
        # Sixth order FD operator

        uₓ[1]   = -1.582533518939116*u[1]   + 2.033378678700676*u[2] - 0.141512858744873*u[3] +
            -0.450398306578272*u[4] + 0.104488069284042*u[5] + 0.036577936277544*u[6]
        uₓ[2]   = -0.462059195631158*u[1] + 0.287258622978251*u[3] + 0.258816087376832*u[4] +
            -0.069112065532624*u[5] - 0.014903449191300*u[6]
        uₓ[3]   = 0.071247104721830*u[1] - 0.636451095137907*u[2] + 0.606235523609147*u[4] +
            -0.022902190275815*u[5] - 0.018129342917256*u[6]
        uₓ[4]   = 0.114713313798970*u[1] - 0.290087484386815*u[2] - 0.306681191361148*u[3] +
            0.520262285050482*u[5] - 0.051642265516119*u[6] + 0.013435342414630*u[7]
        uₓ[5]   = -0.036210680656541*u[1] + 0.105400944933782*u[2] + 0.015764336127392*u[3] +
            -0.707905442575989*u[4] + 0.769199413962647*u[6] - 0.164529643265203*u[7] + 0.018281071473911*u[8]
        uₓ[6]   = -0.011398193015050*u[1] + 0.020437334208704*u[2] + 0.011220896474665*u[3] + 
            0.063183694641876*u[4] - 0.691649024426814*u[5] + 0.739709139060752*u[7] +
            - 0.147941827812150*u[8] + 0.016437980868017*u[9]
        
        uₓ[n]      = 1.582533518939116*u[n]   - 2.033378678700676*u[n-1] + 0.141512858744873*u[n-2] +
            0.450398306578272*u[n-3] - 0.104488069284042*u[n-4] - 0.036577936277544*u[n-5]
        uₓ[n-1]    = 0.462059195631158*u[n] - 0.287258622978251*u[n-2] - 0.258816087376832*u[n-3] +
            0.069112065532624*u[n-4] + 0.014903449191300*u[n-5]
        uₓ[n-2]    = -0.071247104721830*u[n] + 0.636451095137907*u[n-1] - 0.606235523609147*u[n-3] +
            0.022902190275815*u[n-4] + 0.018129342917256*u[n-5]
        uₓ[n-3]    = -0.114713313798970*u[n] + 0.290087484386815*u[n-1] + 0.306681191361148*u[n-2] +
            -0.520262285050482*u[n-4] + 0.051642265516119*u[n-5] - 0.013435342414630*u[n-6]
        uₓ[n-4]    = 0.036210680656541*u[n] - 0.105400944933782*u[n-1] + 0.015764336127392*u[n-2] +
            0.707905442575989*u[n-3] - 0.769199413962647*u[n-5] + 0.164529643265203*u[n-6] - 0.018281071473911*u[n-7]
        uₓ[n-5]    = 0.011398193015050*u[n] - 0.020437334208704*u[n-1] - 0.011220896474665*u[n-2] + 
            -0.063183694641876*u[n-3] + 0.691649024426814*u[n-4] - 0.739709139060752*u[n-6] +
            0.147941827812150*u[n-7] - 0.016437980868017*u[n-8]


        for j = 7:n-6
            uₓ[j] = -1/60*u[j-3] + 3/20*u[j-2] - 3/4*u[j-1] + 3/4*u[j+1] - 3/20*u[j+2] + 1/60*u[j+3]
        end
        
    end

    # Elementwise division by Δx
    uₓ ./= Δx

    return uₓ
end

