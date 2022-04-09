module TimeSolver

    function innerH(u::Vector{Float64},H::Matrix{Float64},v::Vector{Float64})
        return dot(u,H*v)
    end

    function A(uⱼ,PDE::Function,n,Δx,Δy,k,t,x,H,boundary;order=2)

    end

end