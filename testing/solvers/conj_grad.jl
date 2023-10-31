
using LinearAlgebra
using FaADE
using Test


using Random
Random.seed!(1234321)

@warn "testA! must be swtiched on in conj_grad.jl to test"

b = rand(10)

Dom = Grid1D([0.0,1.0],10)

order = 2

K = 1.0

u₀(x) = 0.0

Dl = FaADE.SATs.SAT_Dirichlet(t->0,Dom.Δx,Left,1,order)
Dr = FaADE.SATs.SAT_Dirichlet(t->0,Dom.Δx,Right,1,order)
BD = FaADE.Inputs.SATBoundaries(Dl,Dr)
P = newProblem1D(order,u₀,K,Dom,BD,nothing,nothing)

SC = FaADE.solvers.StepConfig{Float64}()


println("TESTS")

@testset "small identity map" begin
    K = diagm(10,10,ones(10))
    D = FaADE.solvers.newLocalDataBlock(P,Dom,SC,K)
    DB = FaADE.solvers.DataMultiBlock(D)

    DB[1].b = b

    FaADE.solvers.conj_grad!(DB)

    @test all(DB[1].uₙ₊₁ .== DB[1].b)
    @test norm(DB[1].rₖ) < 1e-14
    @test norm(DB[1].K * DB[1].uₙ₊₁ - DB[1].b) < 1e-14
end
# println("b ",DB[1].b)
# println("un ",DB[1].K*DB[1].uₙ₊₁)
# println("resid ",DB[1].K*DB[1].uₙ₊₁ .- DB[1].b," relerr ",norm(DB[1].K*DB[1].uₙ₊₁ .- DB[1].b)/norm(DB[1].b))

@testset "small random matrix" begin
    K = rand(10,10)
    K = K*K' + I
    @testset "small random matrix with exact solution" begin
        D = FaADE.solvers.newLocalDataBlock(P,Dom,SC,K)
        DB = FaADE.solvers.DataMultiBlock(D)

        DB[1].b = b
    
        DB[1].u .= DB[1].K \ DB[1].b
        DB[1].uₙ₊₁ .= DB[1].u
        @test norm(DB[1].K * DB[1].u - DB[1].b) < 1e-14
    
        FaADE.solvers.conj_grad!(DB)
        @test norm(DB[1].rₖ) < 1e-14
        @test norm(DB[1].K * DB[1].uₙ₊₁ - DB[1].b) < 1e-14
    end

    @testset "small random matrix zero guess" begin
        D = FaADE.solvers.newLocalDataBlock(P,Dom,SC,K)
        DB = FaADE.solvers.DataMultiBlock(D)

        DB[1].b = b

        FaADE.solvers.conj_grad!(DB)

        @test norm(DB[1].K*DB[1].uₙ₊₁ - DB[1].b) < 1e-14
        println(DB[1].K*DB[1].uₙ₊₁ - DB[1].b)
    end
end




