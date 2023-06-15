
push!(LOAD_PATH,"../FaADE")
using FaADE


𝒟x = [0.0,1.0]
𝒟y = [-π,π]
nx = 21
ny = 21
Dom = Grid2D(𝒟x,𝒟y,nx,ny)

H = FaADE.Helpers.innerH(Dom,2)


u = rand(Dom.nx,Dom.ny);
v = rand(Dom.nx,Dom.ny);
c = rand(Dom.nx,Dom.ny);

H(u,v)

@code_warntype H(u,v)





χₘₙ = 2.1e-3 + 5.0e-3
params = (ϵₘₙ = [χₘₙ/2.0,χₘₙ/3.0 ],m=[2.0,3.0],n=[1.0,2.0])
function χ_h!(χ,x::Array{Float64},p,t)
    # Hamiltons equations for the field-line Hamiltonian
    # H = ψ²/2 - ∑ₘₙ ϵₘₙ(cos(mθ - nζ))
    χ[2] = x[1] #p_1            qdot        θ
    χ[1] = -sum(p.ϵₘₙ .*(sin.(p.m*x[2] - p.n*t) .* p.m)) #q_1        pdot        ψ
end

dH(X,x,p,t) = χ_h!(X,x,params,t)
PGrid = FaADE.construct_grid(dH,Dom,[-2π,2π])
Pfn = FaADE.generate_parallel_penalty(PGrid,Dom,2)






@code_warntype Pfn(u,v,0.1)




#=

=#
BoundaryLeft = Boundary(Dirichlet,(y,t) -> 0.0,Left,1)
BoundaryRight = Boundary(Dirichlet,(y,t) -> 1.0,Right,1)
BoundaryUpDown = PeriodicBoundary(2)
u₀(x,y) = x
order = 2

P = VariableCoefficientPDE2D(u₀,(x,y)->1e-8,(x,y)->1e-8,order,BoundaryLeft,BoundaryRight,BoundaryUpDown)

typeof(P.BoundaryConditions)







Diff! = FaADE.Derivatives.generate_SecondDerivative(Dom.nx,Dom.ny,Dom.Δx,Dom.Δy,2)

SAT_P = FaADE.SATs.SAT_Periodic(Dom.Δx,2,2)
SAT_P_Fn! = FaADE.SATs.generate_Periodic(SAT_P,:cgie)


SAT_DL = FaADE.SATs.SAT_Dirichlet((x,y)->0.0,Dom.Δx,Left,1,2)
SAT_DR = FaADE.SATs.SAT_Dirichlet((x,y)->0.0,Dom.Δx,Right,1,2)

SAT_DL_Fn! = FaADE.SATs.generate_Dirichlet(SAT_DL,:cgie)
SAT_DR_Fn! = FaADE.SATs.generate_Dirichlet(SAT_DR,:cgie)



# SAT_P_Fn(u,v,c)


function test_SAT_P_Fn!(u::AT,v,c) where AT
    SAT_P_Fn!(u,v,c)
end


# @code_warntype SAT_P_Fn!(u,v,c)
SAT_P(u,v,c)
@code_warntype SAT_P(u,v,c)
# @code_warntype test_SAT_P_Fn!(u,v,c)


mode = FaADE.Helpers.SolutionMode


CGRHS! = let btype1 = P.BoundaryConditions.Left.type,
    btype2 = P.BoundaryConditions.Up.type,
    mode = mode
    function CGRHS!(cache::AT,u::AT,K::KT) where {AT,KT}
        Diff!(cache,u,K,K)
        if (btype1 != Periodic) #Left/Right boundaries
            SAT_DL_Fn!(cache,u,K,mode)
            SAT_DR_Fn!(cache,u,K,mode)
        else
            # SAT_LR!(cache,u,K) #Periodic SAT
        end
        if btype2 != Periodic #Up/Down boundaries
            # SAT_Up!(cache,u,K,mode)
            # SAT_Down!(cache,u,K,mode)
        else
            SAT_P_Fn!(cache,u,K) #Periodic SAT
            # SAT_P(cache,u,K) #Periodic SAT
        end
        cache
    end
end

@code_warntype CGRHS!(u,v,c)






solve(P,Dom,Dom.Δx^2,10Dom.Δx^2,:cgie,adaptive=false,penalty_func=Pfn);

@code_warntype solve(P,Dom,Dom.Δx^2,10Dom.Δx^2,:cgie,adaptive=false,penalty_func=Pfn)




