





#======================================#
#======== SECOND ORDER METHODS ========#
#======================================#


function Dₓₓ(u::Vector{Float64},n::Int64,Δx::Float64,coeff::Vector{Float64};order::Int64=2,boundary::Symbol=:Dirichlet)
    # Function call for the first order SBP operator

    Dₓₓ!(u,u,coeff,n,Δx,order=order)

    return u
end


function Dₓₓ(u::Matrix{Float64},nx::Int64,ny::Int64,Δ::Float64;order::Int64=2,boundary::Symbol=:Dirichlet)

    if dim == 1
        for i = 1:ny
            u[i,:] = Dₓₓ!(u[i,:],u[i,:],coeff,nx,Δx)
        end
    elseif dim == 2
        for i = 1:ny
            u[:,i] = Dₓₓ!(u[:,i],u[:,i],coeff,nx,Δx)
        end
    end
    
end    



#=== SECOND ORDER ITERATOR ===#

function Dₓₓ!(uₓₓ::Vector{Float64},u::Vector{Float64},coeff::Vector{Float64},n::Int64,Δx::Float64;order::Int64=2,boundary::Symbol=:Dirichlet)
    # Iterator for second order SBP operators
    
    if order == 2
        for j = 2:n-1
            uₓₓ[j] = 1/2*( (coeff[j] + coeff[j-1])*u[j-1] - (c[j-1] + 2coeff[j-1])*u[j] + (c[j] + c[j+1])*u[j+1])/Δx^2
        end
    elseif order == 4
    end
end







