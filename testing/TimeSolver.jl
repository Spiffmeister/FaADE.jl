module TimeSolver




function time_solver(u,PDE,SAT,n,x,Δt,Δx,k,t_f,boundary,boundary_left,boundary_right;method=:euler,order=2)

    N = ceil(Int64,t_f/Δt)

    function RHS(uₓₓ,u,n,Δx,Δt,k,t,x,g,boundary_left,boundary_right;order=2)
        # Combines the PDE and SATs
        uₓₓ = PDE(uₓₓ,u,n,Δx,Δt,k,t,x,order=order)
        uₓₓ = SAT(boundary_left,boundary_right,uₓₓ,u,Δx,g,t,k,order=order)
        return uₓₓ
    end
    
    if method == :euler
        for i = 1:N-1
            t = i*Δt
            u[:,i+1] = forward_euler(u[:,i+1],u[:,i],RHS,n,Δx,Δt,k,t,x,order=order)
        end
    elseif method == :rk4
        for i = 1:N-1
            t = i*Δt
            u[:,i+1] = RK4(u[:,i+1],u[:,i],RHS,n,Δx,Δt,k,t,x,order=order)
        end
    elseif method == :impliciteuler
        for i = 1:N-1
            t = i*Δt
            u[:,i+1] = implicit_euler(u[:,i+1],u[:,i],RHS,n,Δx,Δt,t,k,x,order=order)
        end
    elseif method == :cgie
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
            u[:,i+1] = cg(uⱼ,uⱼ,RHS,n,Δx,Δt,k,t,x,boundary;order=2,tol=1e-5,maxIT=20)
        end
    else
        error("Method must be :euler, :rk4, :impliciteuler or :cgie")
    end
    return u
end

#= EXPLICIT METHODS =#

function forward_euler(uₙ,uₒ,RHS,n,Δx,Δt,k,t,x;order=2)
    uₙ = uₒ + Δt*RHS(uₙ,uₒ,n,Δx,Δt,k,t,x,order=order)
    return uₙ
end

function RK4(uₙ,uₒ,RHS,n,Δx,Δt,k,t,x;order=2)
    k1 = RHS(uₙ,uₒ        ,n,Δx,Δt,k,t        ,x,order=order)
    k2 = RHS(uₙ,uₒ+Δt/2*k1,n,Δx,Δt,k,t+0.5Δt  ,x,order=order)
    k3 = RHS(uₙ,uₒ+Δt/2*k2,n,Δx,Δt,k,t+0.5Δt  ,x,order=order)
    k4 = RHS(uₙ,uₒ+Δt*k3  ,n,Δx,Δt,k,t+Δt     ,x,order=order)
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

function CG_implicit_euler()
    
end


#= SUPPORT =#
function A(uⱼ,RHS,n,Δx,Δt,k,t,x;order=2)
    # tmp can be any vector of length(uⱼ)
    tmp = uⱼ
    tmp -= Δt*RHS(tmp,uⱼ,n,Δx,Δt,k,t,x,order=order)
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
    rₖ = A(uⱼ,RHS,n,Δx,Δt,k,t,x,order=order) - b
    dₖ = -rₖ
    i = 0
    rnorm = sqrt(innerH(rₖ,H,rₖ))
    while (norm(rₖ) > tol) & (i < maxIT)
        Adₖ = A(dₖ,RHS,n,Δx,Δt,k,t,x,order=order)
        dₖAdₖ = innerH(dₖ,H,Adₖ)
        αₖ = -innerH(rₖ,H,dₖ)/dₖAdₖ
        xₖ = xₖ + αₖ*dₖ
        rₖ = A(xₖ,RHS,n,Δx,Δt,k,t,x,order=order) - b
        βₖ = innerH(rₖ,H,A(rₖ,RHS,n,Δx,Δt,k,t,x,order=order))/dₖAdₖ
        dₖ = - rₖ + βₖ*dₖ
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
    rnorm = sqrt(innerH(rₖ,H,rₖ))
    while (norm(rₖ) > tol) & (i < maxIT)
        Adₖ = A(dₖ,RHS,n,Δx,Δt,k,t,x,order=order)
        dₖAdₖ = innerH(dₖ,H,Adₖ)
        αₖ = -innerH(rₖ,H,dₖ)/dₖAdₖ
        xₖ = xₖ + αₖ*dₖ
        rₖ = A(xₖ,RHS,n,Δx,Δt,k,t,x,order=order) - b
        βₖ = innerH(rₖ,H,A(rₖ,RHS,n,Δx,Δt,k,t,x,order=order))/dₖAdₖ
        dₖ = - rₖ + βₖ*dₖ
        i += 1
    end
    if (norm(rₖ)>tol)
        warnstr = string("CG did not converge at t=",t)
        @warn warnstr
    end
    return xₖ
end




end