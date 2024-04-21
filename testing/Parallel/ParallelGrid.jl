
using DifferentialEquations
using Test

using Pkg
Pkg.activate(".")
using FaADE


k_para = 1.0
order = 2


@testset "HUDSON AND BRESLAU" begin
    𝒟x = [0.0,1.0]
    𝒟y = [-π,π]
    nx = 21
    ny = 31
    Dom = Grid2D(𝒟x,𝒟y,nx,ny)
    
    
    
    # HUDSON AND BRESLAU
    χₘₙ = 2.1e-3
    params = (ϵₘₙ = [χₘₙ/2., χₘₙ/3.],m=[2.0, 3.0],n=[1.0, 2.0])
    function χ_h!(χ,x::Array{Float64},p,t)
        # Hamiltons equations for the field-line Hamiltonian
        # H = ψ²/2 - ∑ₘₙ ϵₘₙ(cos(mθ - nζ))
        χ[2] = x[1] #p_1            qdot        θ
        χ[1] = -sum(p.ϵₘₙ .*(sin.(p.m*x[2] - p.n*t) .* p.m)) #q_1        pdot        ψ
    end
    
    dH(X,x,p,t) = χ_h!(X,x,params,t)
    gdata = FaADE.construct_grid(dH,Dom,[-2π,2π])
    PData = FaADE.ParallelData(gdata,Dom,order,κ=k_para)
    
    for N in eachindex(Dom)
        # N = 25;
        
        pt = [Dom[N][1],Dom[N][2]];
        test_prob = ODEProblem(dH,pt,(0.0,2π));
        test_soln = DifferentialEquations.solve(test_prob,Tsit5(),reltol=1e-6);
        
        @test isapprox(test_soln.u[end][1],gdata.Fplane.x[N])
        @test isapprox(rem2pi(test_soln.u[end][2],RoundNearest),gdata.Fplane.y[N])
    end
end



# SINGLE ISLAND
@testset "SINGLE ISLAND" begin

    𝒟x = [0.0,1.0]
    𝒟y = [0,2π]
    nx = 21
    ny = 31
    Dom = Grid2D(𝒟x,𝒟y,nx,ny)
    
    
    
    δ=0.005
    xₛ=0.5
    function B(X,x::Array{Float64},p,t)
        X[1] = δ*x[1]^2 * (1.0-x[1]^4) * sin(x[2])
        X[2] = 2.0*x[1] - 2xₛ + 2*δ*(1-x[1]^4)*cos(x[2]) - 4*δ*x[1]^5*cos(x[2])
    end
    dH(X,x,p,t) = B(X,x,params,t)
    gdata = FaADE.construct_grid(dH,Dom,[-2π,2π],ymode=:period)
    PData = FaADE.ParallelData(gdata2,Dom,order,κ=k_para)
    
    
    for N in eachindex(Dom)
        # N = 25;
        
        pt = [Dom[N][1],Dom[N][2]];
        test_prob = ODEProblem(dH,pt,(0.0,2π));
        test_soln = DifferentialEquations.solve(test_prob,Tsit5(),reltol=1e-6);
        
        @test isapprox(test_soln.u[end][1],gdata.Fplane.x[N])
        @test isapprox(mod(test_soln.u[end][2],2π),gdata.Fplane.y[N])
    end

end



