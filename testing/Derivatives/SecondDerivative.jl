# == TESTING FIRST DERIVATIVE ==#
using Test

using LinearAlgebra

using FaADE





function buildgrid(n)
    x = collect(range(0.0,stop=1.0,length=n+1))
    Δx = x[2]-x[1]
    return n+1, x, Δx
end

function buildprob1D(D,f,c,e)
    ue = e.(D.grid)
    u = f.(D.grid)
    k = c.(D.grid)
    return u, k, ue
end

function buildprob2D(D,f,cx,cy,e)
    ue  = zeros(eltype(D.gridx),size(D))
    u   = zeros(eltype(D.gridx),size(D))
    kx  = zeros(eltype(D.gridx),size(D))
    ky  = zeros(eltype(D.gridx),size(D))
    for I in eachindex(D)
        ue[I] = e(D[I]...)
        u[I] = f(D[I]...)
        kx[I] = cx(D[I]...)
        ky[I] = cy(D[I]...)
    end
    return u, kx, ky, ue
end


##======##
# SECOND ORDER
##======##
@testset "1D second derivative second order" begin
    order = 2
    n = 21
    Dom = FaADE.Grid1D([-1.0,1.0],n)
    # x = Dom.grid
    @testset "constant coefficient" begin
        # Linear Function u=x, ∂ₓ(1 ∂ₓx) = 0
        u, k, ue = buildprob1D(Dom,x->x,x->1.0,x->0.0)
        uxx = D₂(u,k,Dom.Δx,order)
        @test norm(uxx[2:end-1] .- ue[2:end-1]) .≤ 1e-12
    
        # Quadratic function, ∂ₓ(1 ∂ₓx²) = 2
        u, k, ue = buildprob1D(Dom,x->x^2,x->1.0,x->2.0)
        uxx = D₂(u,k,Dom.Δx,order)
        @test norm(uxx[2:end-1] .- ue[2:end-1]) ≤ 1.0e-12
        @test (uxx[1] == 0) & (uxx[end] == 0)
    
        # Cubic function, ∂ₓ(1 ∂ₓx³) = 6x
        u, k, ue = buildprob1D(Dom,x->x^3,x->1.0,x->6x)
        uxx = D₂(u,k,Dom.Δx,order)
        @test norm(uxx[2:end-1] .- ue[2:end-1]) ≤ 1.0e-12
        @test (uxx[1] == 0) & (uxx[end] == 0)
    end

    @testset "variable coefficient" begin

        u, k, ue = buildprob1D(Dom,x->x,x->x,x->1.0)
        uxx = D₂(u,k,Dom.Δx,order)
        @test norm(uxx[2:end-1] .- ue[2:end-1]) ≤ 1.0e-12
        @test (uxx[1] == 0) & (uxx[end] == 0)
    
        # Quadratic Function, ∂ₓ(x ∂ₓx²) = 4x
        u, k, ue = buildprob1D(Dom,x->x^2,x->x,x->4x)
        uxx = D₂(u,k,Dom.Δx,order)
        @test norm(uxx[2:end-1] .- ue[2:end-1]) ≤ 1.0e-12
        @test (uxx[1] == 0) & (uxx[end] == 0)

    
        # Cubic Function, ∂ₓ(x ∂ₓx³) = 9x²
        u, k, ue = buildprob1D(Dom,x->x^3,x->x,x->9x^2)
        uxx = D₂(u,k,Dom.Δx,order)
        @test norm(uxx[2:end-1] .- ue[2:end-1]) ≤ 1.0e-12 broken = true
        @test (uxx[1] == 0) & (uxx[end] == 0)
    end
end

@testset "2D second derivative second order" begin
    nx = 21
    ny = 31
    order = 2
    Dom = FaADE.Grid2D([-1.0,1.0],[0.0,1.0],nx,ny)

    @testset "constant coefficient" begin
        # Linear Function u=x, ∇(1 ∇xy) = 0
        u,kx,ky,ue = buildprob2D(Dom,(x,y)->x*y,(x,y)->1.0,(x,y)->1.0,(x,y)->0.0)
        uxx = D₂(u,kx,ky,Dom.Δx,Dom.Δy,order)
        @test norm(uxx[2:end-1,2:end-1] .- ue[2:end-1,2:end-1]) ≤ 1.0e-12
        @test all(uxx[[1,Dom.nx],:] .≤ 1e-12) & all(uxx[:,[1,Dom.ny]] .≤ 1e-12)
        
        # Quadratic function, ∇(1 ∇x²y) = 2
        u,kx,ky,ue = buildprob2D(Dom,(x,y)->x^2*y,(x,y)->1.0,(x,y)->1.0,(x,y)->2.0y)
        uxx = D₂(u,kx,ky,Dom.Δx,Dom.Δy,order)
        @test norm(uxx[2:end-1,2:end-1] .- ue[2:end-1,2:end-1]) ≤ 1.0e-12
        @test all(uxx[[1,Dom.nx],:] .≤ 1e-12)
        @test all(uxx[:,1] .≤ 1e-12) & (norm(uxx[2:end-1,end] .- 2) ≤ 1e-12)

        # Cubic Function u=x, ∇(1 ∇x²y³) = 6y + 2y^3
        u,kx,ky,ue = buildprob2D(Dom,(x,y)->x^2*y^3,(x,y)->1.0,(x,y)->1.0,(x,y)->6y*x^2+2y^3)
        uxx = D₂(u,kx,ky,Dom.Δx,Dom.Δy,order)
        @test norm(uxx[2:end-1,2:end-1] .- ue[2:end-1,2:end-1]) ≤ 1.0e-12 # CHECK THIS TEST
    end

    # @testset "variable coefficient" begin
    # end

end




##======##
# FOURTH ORDER
##======##
@testset "1D second derivative fourth order" begin
    order = 4
    n = 21
    Dom = FaADE.Grid1D([-1.0,1.0],n)
    # x = Dom.grid
    @testset "constant coefficient" begin
        # Linear Function u=x, ∂ₓ(1 ∂ₓx) = 0
        u, k, ue = buildprob1D(Dom,x->x,x->1.0,x->0.0)
        uxx = D₂(u,k,Dom.Δx,order)
        @test norm(uxx .- ue) .≤ 1e-12
    
        # Quadratic function, ∂ₓ(1 ∂ₓx²) = 2
        u, k, ue = buildprob1D(Dom,x->x^2,x->1.0,x->2.0)
        uxx = D₂(u,k,Dom.Δx,order)
        @test norm(uxx .- ue) ≤ 1.0e-12
    
        # Cubic function, ∂ₓ(1 ∂ₓx³) = 6x
        u, k, ue = buildprob1D(Dom,x->x^3,x->1.0,x->6x)
        uxx = D₂(u,k,Dom.Δx,order)
        @test norm(uxx .- ue) ≤ 1.0e-12
    end

    @testset "variable coefficient" begin

        u, k, ue = buildprob1D(Dom,x->x,x->x,x->1.0)
        uxx = D₂(u,k,Dom.Δx,order)
        @test norm(uxx[2:end-1] .- ue[2:end-1]) ≤ 1.0e-12
    
        # Quadratic Function, ∂ₓ(x ∂ₓx²) = 4x
        u, k, ue = buildprob1D(Dom,x->x^2,x->x,x->4x)
        uxx = D₂(u,k,Dom.Δx,order)
        @test norm(uxx[2:end-1] .- ue[2:end-1]) ≤ 1.0e-12

    
        # Cubic Function, ∂ₓ(x ∂ₓx³) = 9x²
        u, k, ue = buildprob1D(Dom,x->x^3,x->x,x->9x^2)
        uxx = D₂(u,k,Dom.Δx,order)
        @test norm(uxx[2:end-1] .- ue[2:end-1]) ≤ 1.0e-12 broken = true
    end
end

#=
@testset "2D second derivative second order" begin
    nx = 21
    ny = 31
    order = 4
    Dom = FaADE.Grid2D([-1.0,1.0],[0.0,1.0],nx,ny)

    @testset "constant coefficient" begin
        # Linear Function u=x, ∇(1 ∇xy) = 0
        u,kx,ky,ue = buildprob2D(Dom,(x,y)->x*y,(x,y)->1.0,(x,y)->1.0,(x,y)->0.0)
        uxx = D₂(u,kx,ky,Dom.Δx,Dom.Δy,order)
        @test norm(uxx[2:end-1,2:end-1] .- ue[2:end-1,2:end-1]) ≤ 1.0e-12
        @test all(uxx[[1,Dom.nx],:] .≤ 1e-12) & all(uxx[:,[1,Dom.ny]] .≤ 1e-12)
        
        # Quadratic function, ∇(1 ∇x²y) = 2
        u,kx,ky,ue = buildprob2D(Dom,(x,y)->x^2*y,(x,y)->1.0,(x,y)->1.0,(x,y)->2.0y)
        uxx = D₂(u,kx,ky,Dom.Δx,Dom.Δy,order)
        @test norm(uxx[2:end-1,2:end-1] .- ue[2:end-1,2:end-1]) ≤ 1.0e-12
        @test all(uxx[[1,Dom.nx],:] .≤ 1e-12)
        @test all(uxx[:,1] .≤ 1e-12) & (norm(uxx[2:end-1,end] .- 2) ≤ 1e-12)

        # Cubic Function u=x, ∇(1 ∇x²y³) = 6y + 2y^3
        u,kx,ky,ue = buildprob2D(Dom,(x,y)->x^2*y^3,(x,y)->1.0,(x,y)->1.0,(x,y)->6y*x^2+2y^3)
        uxx = D₂(u,kx,ky,Dom.Δx,Dom.Δy,order)
        @test norm(uxx[2:end-1,2:end-1] .- ue[2:end-1,2:end-1]) ≤ 1.0e-12 # CHECK THIS TEST
    end

    # @testset "variable coefficient" begin
    # end

end
=#








#=
nx = ny = 21
order = FaADE.Derivatives.DerivativeOrder{2}()
Dom = FaADE.Grid2D([-1.0,1.0],[-1.0,1.0],nx,ny)
D = FaADE.Derivatives.DerivativeOperator{Float64,2,typeof(order),:Constant}(order,nx,ny,Dom.Δx,Dom.Δy,true,true)

uxx, u, kx, ky, ue = buildprob2D(Dom,(x,y)->sin(2π*x)*sin(2π*y),(x,y)->1.0,(x,y)->1.0,(x,y)->-8π^2*sin(2π*x)*sin(2π*y))
FaADE.Derivatives.mul!(uxx,u,[kx,ky],D)
=#




