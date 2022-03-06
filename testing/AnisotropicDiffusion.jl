#=
Testing


∂ₜu = ∂ₓ(k ∂ₓu) = ∂ₓk ∂ₓu + k ∂ₓ∂ₓu
=#


# Initalise the things
n = 100

domain = [0,2π]

x = collect(range(domain[1],domain[2],length=n))

u = InitialCondition(x)

struct soln
    x :: Vector{Float64}    # Grid
    u :: Matrix{Float64}    # Solution
    Δt :: Float64           # Step size
end





function InitialCondition(x)
    return -cos(x)
end


function k(x)
    return x
end




function timesolve()
end