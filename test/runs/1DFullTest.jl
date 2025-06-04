using FaADE


###
𝒟 = [0.0,1.0]
n = 101
Dom = Grid1D(𝒟,n)

K(x) = 1.0

Δt = 0.01
t_f = 200.0

u₀(x) = sin.(2π*x*2 .+ 1.0)

order = 2



@testset "Dirichlet boundaries" begin
    BoundaryLeft =  SAT_Dirichlet(t->sin(1.0),       Dom.Δx, Left, order)
    BoundaryRight = SAT_Dirichlet(t->sin(4π + 1.0),  Dom.Δx, Right, order)
    
    BD = (BoundaryLeft,BoundaryRight)
    P = Problem1D(order,u₀,K,Dom,BD)
    soln = solve(P,Dom,Δt,t_f)
    
    @test all(isapprox.(soln.u[2],sin(1.0),atol=1e-6))
end


@testset "Periodic boundaries" begin
    BoundaryLeft = SAT_Periodic(Dom.Δx,Left,order)
    BoundaryRight = SAT_Periodic(Dom.Δx,Right,order)
    
    BD = (BoundaryLeft,BoundaryRight)
    P = Problem1D(order,u₀,K,Dom,BD)
    soln = solve(P,Dom,Δt,t_f)
   
    @test all(isapprox.(soln.u[2],0.0,atol=1e-12))
end


# New initial condition
u₀(x) = sin.(2π*x)

@testset "Neumann boundaries" begin
    BoundaryLeft = FaADE.SATs.SAT_Neumann(t->0.0, Dom.Δx,Left,order)
    BoundaryRight = FaADE.SATs.SAT_Neumann(t->0.0, Dom.Δx,Right,order)
    
    BD = (BoundaryLeft,BoundaryRight)
    P = Problem1D(order,u₀,K,Dom,BD)
    soln = solve(P,Dom,Δt,t_f)
   
    @test all(isapprox.(soln.u[2],0.0,atol=1e-12))
end


@testset "Robin boundaries" begin
    BoundaryLeft = SAT_Robin(t->0.0,Dom.Δx,Left, order)
    BoundaryRight = SAT_Robin(t->0.0,Dom.Δx,Right,order)

    BD = (BoundaryLeft,BoundaryRight)
    P = Problem1D(order,u₀,K,Dom,BD)
    soln = solve(P,Dom,Δt,t_f)
   
    @test all(isapprox.(soln.u[2],0.0,atol=1e-12))
end


