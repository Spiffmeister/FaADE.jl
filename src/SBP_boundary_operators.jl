

function boundary_Dₓ(u::Vector{Float64},Δx::Float64,order=2)
    # Implementation of the 

    uₓ = zeros(Float64,2*order)

    if order == 2
        #Boundary for the second order case
        
        uₓ[1] = u[1]/Δx
        uₓ[2] = u[1]/Δx
        
        uₓ[end-1] = u[1]/Δx
        uₓ[end] = u[1]/Δx
        
    elseif order == 4
        #Boundary for the fourth order case

        BOp = [-24.0/17.0, 59.0/34.0, -4.0/17.0, -3.0/34.0]
        
        uₓ[1:4] = BOp*u[1]/Δx
        
        uₓ[end-3] = -BOp*u[1]/Δx
        
        
    elseif order == 6
        #Boundary for the sixth order case

        BOp = [-1.582533518939116, 2.033378678700676, -0.141512858744873, -0.450398306578272, 0.104488069284042, 0.036577936277544]

        uₓ[1:6] = BOp*u[1]/Δx

        uₓ[end-5:end] = -BOp*u[1]/Δx

    end

    return uₓ
end