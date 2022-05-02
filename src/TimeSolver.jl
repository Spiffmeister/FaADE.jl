

struct prob
    RHS             :: Function
    boundary_left   :: Symbol
    boundary_right  :: Symbol
    boundary        :: Union{Vector{Float64},Function}
    order           :: Int64
    k               :: Union{Vector{Float64},Float64}
end

mutable struct solution
    u       :: Array{Float64}
    x       :: Array{Float64}
    Δx      :: Float64
    n       :: Int64
    t_f     :: Float64
    Δt      :: Float64
    solver  :: Symbol

    function solution(u₀,x,Δx,t_f,Δt,solver)
        N = ceil(Int64,t_f/Δt)
        n = length(x)

        u = zeros(Float64,n,N)
        u[:,1] = u₀(x)

        new(u,x,Δx,n,t_f,Δt,solver)
    end
end



"""
"""
function time_solver(PDE::Function,u₀::Function,n::Int64,x::Vector{Float64},Δx::Float64,t_f::Float64,Δt::Float64,k::Vector{Float64},boundary::Function,boundary_left::Symbol;boundary_right::Symbol=boundary_left,method::Symbol=:euler,order::Int64=2)
    #=
    Expects two functions PDE and SAT from the SBP_operators package

    =#

    # Initialise solution
    soln = solution(u₀,x,Δx,t_f,Δt,method)

    # Get the length of the time array
    N = ceil(Int64,t_f/Δt)

    if method != :cgie
        function RHS(uₓₓ,u,n,x,Δx,t,Δt,k,g)
            # Combines the PDE and SATs (forcing term included)
            uₓₓ = PDE(uₓₓ,u,n,x,Δx,t,Δt,k,order=order)

            uₓₓ[1:order] .+= SAT_left(boundary_left,u,Δx,g(t),c=k,order=order)
            uₓₓ[end-order+1:end] .+= SAT_right(boundary_right,u,Δx,g(t),c=k,order=order)
            
            return uₓₓ
        end
    # else
        # function RHS(uₓₓ,u,n,Δx,k,t,x,g)

        # end
    end

    if method == :euler
        for i = 1:N-1
            t = i*Δt
            soln.u[:,i+1] = forward_euler(soln.u[:,i+1],soln.u[:,i],RHS,n,Δx,Δt,k,t,x,boundary)
        end
    elseif method == :rk4
        for i = 1:N-1
            t = i*Δt
            u[:,i+1] = RK4(u[:,i+1],u[:,i],RHS,n,Δx,Δt,k,t,x,boundary,order=order)
        end
    elseif method == :impliciteuler
        for i = 1:N-1
            t = i*Δt
            u[:,i+1] = implicit_euler(u[:,i+1],u[:,i],RHS,n,Δx,Δt,t,k,x,order=order)
        end
    elseif method == :cgie

        function RHS()
            uₓₓ = PDE(uₓₓ,u,n,Δx,Δt,k,t,x,order=order)
            SATₗ, = SAT_left(boundary_left,u,Δx,g,order=order,seperate_forcing=true)
            SATᵣ, = SAT_right(boundary_left,u,Δx,g,order=order,seperate_forcing=true)
            uₓₓ[1:order] = SAT(boundary_left,boundary_right,uₓₓ,u,Δx,g,t,k,order=order,seperate_forcing=true)
            return uₓₓ
        end

        if order == 2
            H = diagm(ones(length(u)))
            H[1,1] = H[end,end] = 0.5
        elseif order == 4
            H = diagm(ones(length(u)))
            H[1,1] = H[end,end] = 0.5
            H[2,2] = H[end-1,end-1] = 0.5
            H[3,3] = H[end-2,end-2] = 0.5
        elseif order == 6
            H = diagm(ones(length(u)))
        end
        H = Δx*H

        for i = 1:N-1
            t = i*Δt
            uⱼ = u[:,i]
            u[:,i+1] = conj_grad(uⱼ,uⱼ,RHS,n,Δx,Δt,k,t,x,boundary;order=2,tol=1e-5,maxIT=20)
            SATₗ,Fₗ = SAT_left(boundary_left,u,Δx,g,order=order,seperate_forcing=true)
            SATᵣ,Fᵣ = SAT_right(boundary_right,u,Δx,g,order=order,seperate_forcing=true)
            u[1:order,i+1] += Fₗ
            u[order,i+1] += Fₗ
        end
    else
        error("Method must be :euler, :rk4, :impliciteuler or :cgie")
    end
    return soln
end


# function time_solver(u,PDE,SAT,nx,ny,Δx,Δy,x,y,Δt,k,boundary_x,boundary_y;method=:euler,order=2,maxIT=5)
# end



#= EXPLICIT METHODS =#

function forward_euler(uₙ::Vector,uₒ::Vector,RHS::Function,n::Int,Δx::Float64,Δt::Float64,k::Vector,t,x::Vector,g)
    # Simple forward euler method
    uₙ = uₒ + Δt*RHS(uₙ,uₒ,n,x,Δx,t,Δt,k,g)
    return uₙ
end

function RK4(uₙ,uₒ,RHS,n,Δx,Δt,k,t,x,boundary;order=2)
    k1 = RHS(uₙ,uₒ        ,n,Δx,Δt,k,t        ,x,boundary,order=order)
    k2 = RHS(uₙ,uₒ+Δt/2*k1,n,Δx,Δt,k,t+0.5Δt  ,x,boundary,order=order)
    k3 = RHS(uₙ,uₒ+Δt/2*k2,n,Δx,Δt,k,t+0.5Δt  ,x,boundary,order=order)
    k4 = RHS(uₙ,uₒ+Δt*k3  ,n,Δx,Δt,k,t+Δt     ,x,boundary,order=order)
    uₙ = uₒ + Δt/6 * (k1 + 2k2 + 2k3 + k4)
    return uₙ
end


#= IMPLICIT METHODS =#

function implicit_euler(uₙ,uₒ,RHS,n,Δx,Δt,k,t,x;order=2,maxIT=100,α=1.5)
    uⱼ = uₒ
    for j = 1:maxIT
        uⱼ = uⱼ - α * Δt * (uⱼ - uₒ - Δt*RHS(uₙ,uⱼ,n,Δx,Δt,k,t,x,order=order))
    end
    uⱼ = uₒ + Δt*RHS(uₙ,uⱼ,n,Δx,Δt,k,t,x,order=order)
    return uⱼ
end




#= SUPPORT =#
function A(uⱼ,PDE::Function,n,Δx,x,Δt,t,k;order=2)
    # tmp can be any vector of length(uⱼ)
    tmp = uⱼ
    tmp -= Δt*PDE(tmp,uⱼ,n,Δx,Δt,k,t,x,order=order)
    return tmp
end

function A(uⱼ,PDE::Function,nx,ny,Δx,x,Δy,y,Δt,t,k,boundary_x,boundary_y;order=2)
    # A for 2D arrays
    tmp = uⱼ
    tmp -= Δt*PDE(tmp,uⱼ,nx,ny,Δx,x,Δy,y,Δt,t,k,boundary_x,boundary_y,order=order)
    return tmp
end


function innerH(u::Vector,H::Matrix,v::Vector)
    # H inner product for 1D problems
    return dot(u,H*v)
end

function innerH(u::Matrix,Hx::Vector,Hy::Vector,v::Matrix)
    # H inner product for 2D problems
    nx,ny = size(u)
    tmp = 0
    for i = 1:nx
        for j = 1:ny
            tmp += u[i,j]*v[i,j]*Hx[i]*Hy[j]
        end
    end
    return tmp
end




function conj_grad(b,uⱼ,RHS,n,Δx,Δt,k,t,x,H;order=2,tol=1e-5,maxIT=10)
    xₖ = zeros(length(b)) #Initial guess
    rₖ = A(uⱼ,RHS,n,Δx,x,Δt,t,k,order=order) - b
    dₖ = -rₖ
    i = 0
    rnorm = sqrt(innerH(rₖ,H,rₖ))
    while (rnorm > tol) & (i < maxIT)
        Adₖ = A(dₖ,RHS,n,Δx,Δt,k,t,x,order=order)
        dₖAdₖ = innerH(dₖ,H,Adₖ)
        αₖ = -innerH(rₖ,H,dₖ)/dₖAdₖ
        xₖ = xₖ + αₖ*dₖ
        rₖ = A(xₖ,RHS,n,Δx,x,Δt,t,k,order=order) - b
        βₖ = innerH(rₖ,H,A(rₖ,RHS,n,Δx,Δt,k,t,x,order=order))/dₖAdₖ
        dₖ = - rₖ + βₖ*dₖ
        rnorm = sqrt(innerH(rₖ,H,rₖ))
        i += 1
    end
    if (norm(rₖ)>tol)
        warnstr = string("CG did not converge at t=",t)
        @warn warnstr
    end
    return xₖ
end

function conj_grad(b,uⱼ,RHS,nx,ny,Δx,x,Δy,y,Δt,t,k,Hx,Hy;order=2,tol=1e-5,maxIT=10)
    xₖ = zeros(length(b)) #Initial guess
    rₖ = A(uⱼ,RHS,n,Δx,Δt,k,t,x,order=order) - b
    dₖ = -rₖ
    i = 0
    rnorm = sqrt(innerH(rₖ,Hx,Hy,rₖ))
    while (rnorm > tol) & (i < maxIT)
        Adₖ = A(dₖ,RHS,nx,ny,Δx,x,Δy,y,Δt,t,k,boundary_x,boundary_y,order=order)
        dₖAdₖ = innerH(dₖ,H,Adₖ)
        αₖ = -innerH(rₖ,H,dₖ)/dₖAdₖ
        xₖ = xₖ + αₖ*dₖ
        rₖ = A(xₖ,RHS,n,Δx,Δt,k,t,x,order=order) - b
        βₖ = innerH(rₖ,H,A(rₖ,RHS,n,Δx,Δt,k,t,x,order=order))/dₖAdₖ
        dₖ = - rₖ + βₖ*dₖ
        rnorm = sqrt(innerH(rₖ,Hx,Hy,rₖ))
        i += 1
    end
    if (norm(rₖ)>tol)
        warnstr = string("CG did not converge at t=",t)
        @warn warnstr
    end
    return xₖ
end


