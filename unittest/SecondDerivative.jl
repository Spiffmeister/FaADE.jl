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
    uxx = zeros(eltype(D.grid),D.n)
    k = c.(D.grid)
    return uxx, u, k, ue
end

function buildprob2D(D,f,cx,cy,e)
    ue  = zeros(eltype(D.gridx),(D.nx,D.ny))
    u   = zeros(eltype(D.gridx),(D.nx,D.ny))
    kx  = zeros(eltype(D.gridx),(D.nx,D.ny))
    ky  = zeros(eltype(D.gridx),(D.nx,D.ny))
    for i in eachindex(D.gridx)
        for j in eachindex(D.gridy)
            ue[i,j] = e(D.gridx[i],D.gridy[j])
            u[i,j]  = f(D.gridx[i],D.gridy[j])
            kx[i,j] = cx(D.gridx[i],D.gridy[j])
            ky[i,j] = cy(D.gridx[i],D.gridy[j])
        end
    end
    uxx = zeros(eltype(D.gridx),(D.nx,D.ny))
    return uxx, u, kx, ky, ue
end

# function T1DD2(n,D,f,c,e)
#     Dom = FaADE.Grid1D([-1.0,1.0],n)
#     buildprob1D(D,f,c,e)
#     D₂!(uxx,u,k,Dom.n,Dom.Δx,order,0.0)
#     chk1 = norm(∂ₓₓu[2:end-1] .- ue[2:end-1])
#     chk2 = (uxx[1] == 0) & (uxx[end] == 0)
#     return chk1 & chk2
# end

##======##
# SECOND ORDER
##======##
@testset "1D second derivative second order" begin
    order = FaADE.Derivatives.DerivativeOrder{2}()
    n = 21
    Dom = FaADE.Grid1D([-1.0,1.0],n)
    # x = Dom.grid
    @testset "constant coefficient" begin
        # Linear Function u=x, ∂ₓ(1 ∂ₓx) = 0
        uxx, u, k, ue = buildprob1D(Dom,x->x,x->1.0,x->0.0)
        D₂!(uxx,u,k,Dom.n,Dom.Δx,order,0.0)
        @test norm(uxx[2:end-1] .- ue[2:end-1]) .≤ 1e-12
    
        # Quadratic function, ∂ₓ(1 ∂ₓx²) = 2
        uxx, u, k, ue = buildprob1D(Dom,x->x^2,x->1.0,x->2.0)
        D₂!(uxx,u,k,Dom.n,Dom.Δx,order,0.0)
        @test norm(uxx[2:end-1] .- ue[2:end-1]) ≤ 1.0e-12
        @test (uxx[1] == 0) & (uxx[end] == 0)
    
        # Cubic function, ∂ₓ(1 ∂ₓx³) = 6x
        uxx, u, k, ue = buildprob1D(Dom,x->x^3,x->1.0,x->6x)
        D₂!(uxx,u,k,Dom.n,Dom.Δx,order,0.0)
        @test norm(uxx[2:end-1] .- ue[2:end-1]) ≤ 1.0e-12
        @test (uxx[1] == 0) & (uxx[end] == 0)
    end

    @testset "variable coefficient" begin

        uxx, u, k, ue = buildprob1D(Dom,x->x,x->x,x->1.0)
        D₂!(uxx,u,k,Dom.n,Dom.Δx,order,0.0)
        @test norm(uxx[2:end-1] .- ue[2:end-1]) ≤ 1.0e-12
        @test (uxx[1] == 0) & (uxx[end] == 0)
    
        # Quadratic Function, ∂ₓ(x ∂ₓx²) = 4x
        uxx, u, k, ue = buildprob1D(Dom,x->x^2,x->x,x->4x)
        D₂!(uxx,u,k,Dom.n,Dom.Δx,order,0.0)
        @test norm(uxx[2:end-1] .- ue[2:end-1]) ≤ 1.0e-12
        @test (uxx[1] == 0) & (uxx[end] == 0)

    
        # Cubic Function, ∂ₓ(x ∂ₓx³) = 9x²
        uxx, u, k, ue = buildprob1D(Dom,x->x^3,x->x,x->9x^2)
        D₂!(uxx,u,k,Dom.n,Dom.Δx,order,0.0)
        @test norm(uxx[2:end-1] .- ue[2:end-1]) ≤ 1.0e-12 broken = true
        @test (uxx[1] == 0) & (uxx[end] == 0)
    end
end

@testset "2D second derivative second order" begin
    nx = 21
    ny = 31
    order = FaADE.Derivatives.DerivativeOrder{2}()
    Dom = FaADE.Grid2D([-1.0,1.0],[0.0,1.0],nx,ny)

    @testset "constant coefficient" begin
        # Linear Function u=x, ∇(1 ∇xy) = 0
        uxx,u,kx,ky,ue = buildprob2D(Dom,(x,y)->x*y,(x,y)->1.0,(x,y)->1.0,(x,y)->0.0)
        D₂!(uxx,u,kx,ky,Dom.nx,Dom.ny,Dom.Δx,Dom.Δy,order,order)
        @test norm(uxx[2:end-1,2:end-1] .- ue[2:end-1,2:end-1]) ≤ 1.0e-12
        @test all(uxx[[1,Dom.nx],:] .≤ 1e-12) & all(uxx[:,[1,Dom.ny]] .≤ 1e-12)
        
        # Quadratic function, ∇(1 ∇x²y) = 2
        uxx,u,kx,ky,ue = buildprob2D(Dom,(x,y)->x^2*y,(x,y)->1.0,(x,y)->1.0,(x,y)->2.0y)
        D₂!(uxx,u,kx,ky,Dom.nx,Dom.ny,Dom.Δx,Dom.Δy,order,order)
        @test norm(uxx[2:end-1,2:end-1] .- ue[2:end-1,2:end-1]) ≤ 1.0e-12
        @test all(uxx[[1,Dom.nx],:] .≤ 1e-12)
        @test all(uxx[:,1] .≤ 1e-12) & (norm(uxx[2:end-1,end] .- 2) ≤ 1e-12)

        # Cubic Function u=x, ∇(1 ∇x²y³) = 6y + 2y^3
        uxx,u,kx,ky,ue = buildprob2D(Dom,(x,y)->x^2*y^3,(x,y)->1.0,(x,y)->1.0,(x,y)->6y*x^2+2y^3)
        D₂!(uxx,u,kx,ky,Dom.nx,Dom.ny,Dom.Δx,Dom.Δy,order,order)
        @test norm(uxx[2:end-1,2:end-1] .- ue[2:end-1,2:end-1]) ≤ 1.0e-12 broken = true
        @test all(uxx[[1,Dom.nx],:] .≤ 1e-12) & all(uxx[:,[1,Dom.ny]] .≤ 1e-12)
    end

    @testset "variable coefficient" begin
    end

end



#=

# # Quartic function, ∂ₓ(1 ∂ₓx⁴) = 12x²
# u = x.^4
# ∂ₓₓuₑ = 12x.^2
# ∂ₓₓu = D₂(u,n,Δx,c)

# lowres = norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1])

# n, x, Δx = buildgrid(20)
# c = ones(n)
# u = x.^4
# ∂ₓₓuₑ = 12x.^2
# ∂ₓₓu = D₂(u,n,Δx,c)

# highres = norm(∂ₓₓuₑ[2:n-1] .- ∂ₓₓu[2:n-1])

# (lowres-highres)/highres



#= 2D methods =#
@testset "Dxx second order 2D" begin
    nx = ny = 21
    Dx = Dy = [-1.0,1.0]
    Dom = FaADE.Grid2D(Dx,Dy,nx,ny)


end








##======##
# FOURTH ORDER
##======##

### Constant coefficient
# Linear Function
n, x, Δx = buildgrid(15)
c = ones(n)
u = x
∂ₓₓuₑ = zeros(n)
∂ₓₓu = D₂(u,n,Δx,c,order=4)
@test norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-13

# Quadratic function
n, x, Δx = buildgrid(15)
c = ones(n)
u = x.^2
∂ₓₓuₑ = 2*ones(n)
∂ₓₓu = D₂(u,n,Δx,c,order=4)
@test norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-12

# Cubic function
n, x, Δx = buildgrid(15)
c = ones(n)
u = x.^3
∂ₓₓuₑ = 6x
∂ₓₓu = D₂(u,n,Δx,c,order=4)
@test norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-13

# Quartic function
n, x, Δx = buildgrid(15)
c = ones(n)
u = x.^4
∂ₓₓuₑ = 12x.^2
∂ₓₓu = D₂(u,n,Δx,c,order=4)
@test norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-13

# Quintic function
n, x, Δx = buildgrid(15)
c = ones(n)
u = x.^5
∂ₓₓuₑ = 20x.^3
∂ₓₓu = D₂(u,n,Δx,c,order=4)
@test norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-13

# Hextic function
n, x, Δx = buildgrid(15)
c = ones(n)
u = x.^6
∂ₓₓuₑ = 30x.^4
∂ₓₓu = D₂(u,n,Δx,c,order=4)
@test_broken norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-12


### Variable coefficient
# Linear function
n, x, Δx = buildgrid(15)
c = x
u = x
∂ₓₓuₑ = ones(n)
∂ₓₓu = D₂(u,n,Δx,c,order=4)
@test norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-13

# Quadratic function
n, x, Δx = buildgrid(15)
c = x
u = x.^2
∂ₓₓuₑ = 4x
∂ₓₓu = D₂(u,n,Δx,c,order=4)
@test norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-13

# Cubic function
n, x, Δx = buildgrid(15)
c = x
u = x.^3
∂ₓₓuₑ = 9x.^2
∂ₓₓu = D₂(u,n,Δx,c,order=4)
@test norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-13

# BROKEN TESTS
# Quartic function
n, x, Δx = buildgrid(15)
c = x
u = x.^4
∂ₓₓuₑ = 16x.^3
∂ₓₓu = D₂(u,n,Δx,c)
@test_broken norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-14


# Quadratic coefficient
n, x, Δx = buildgrid(15)
c = x.^2
u = x
∂ₓₓuₑ = ones(n)
∂ₓₓu = D₂(u,n,Δx,c,order=4)
@test_broken norm(∂ₓₓuₑ[5:n-4] .- ∂ₓₓu[5:n-4]) ≤ 1.0e-13

#=
##======##
# SIXTH ORDER
##======##

### Constant coefficient
# Linear Function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x
∂ₓₓuₑ = zeros(n)
∂ₓₓu = D₂(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ[6:n-5] .- ∂ₓₓu[6:n-5]) ≤ 1.0e-12

# Quadratic function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^2
∂ₓₓuₑ = 2*ones(n)
∂ₓₓu = D₂(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ .- ∂ₓₓu) ≤ 1.0e-10

# Cubic function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^3
∂ₓₓuₑ = 6x
∂ₓₓu = D₂(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ .- ∂ₓₓu) ≤ 1.0e-10

# Quartic function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^4
∂ₓₓuₑ = 12x.^2
∂ₓₓu = D₂(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ .- ∂ₓₓu) ≤ 1.0e-12

# Quintic function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^5
∂ₓₓuₑ = 20x.^3
∂ₓₓu = D₂(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ .- ∂ₓₓu) ≤ 1.0e-12

# Hextic function - only body is 6th order
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^6
∂ₓₓuₑ = 30x.^4
∂ₓₓu = D₂(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ .- ∂ₓₓu) ≤ 1.0e-13

# Heptic coefficient
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^7
∂ₓₓuₑ = 42x.^5
∂ₓₓu = D₂(u,n,Δx,c,order=6)

@test_broken norm(∂ₓₓuₑ .- ∂ₓₓu) ≤ 1.0e-13



### VARIABLE COEFFICIENT
# Linear Function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x
∂ₓₓuₑ = zeros(n)
∂ₓₓu = D₂(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ[6:n-5] .- ∂ₓₓu[6:n-5]) ≤ 1.0e-12

# Quadratic function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^2
∂ₓₓuₑ = 2*ones(n)
∂ₓₓu = D₂(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ[6:n-5] .- ∂ₓₓu[6:n-5]) ≤ 1.0e-12

# Cubic function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^3
∂ₓₓuₑ = 6x
∂ₓₓu = D₂(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ[6:n-5] .- ∂ₓₓu[6:n-5]) ≤ 1.0e-12

# Quartic function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^4
∂ₓₓuₑ = 12x.^2
∂ₓₓu = D₂(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ[6:n-5] .- ∂ₓₓu[6:n-5]) ≤ 1.0e-12

# Quintic function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^5
∂ₓₓuₑ = 20x.^3
∂ₓₓu = D₂(u,n,Δx,c,order=6)

@test norm(∂ₓₓuₑ[7:n-6] .- ∂ₓₓu[7:n-6]) ≤ 1.0e-12

# Hextic function
n, x, Δx = buildgrid(20)
c = ones(n)
u = x.^6
∂ₓₓuₑ = 30x.^4
∂ₓₓu = D₂(u,n,Δx,c,order=6)
@test norm(∂ₓₓuₑ[7:n-6] .- ∂ₓₓu[7:n-6]) ≤ 1.0e-12

# Heptic coefficient
n, x, Δx = buildgrid(20)
c = x
u = x.^7
∂ₓₓuₑ = 42x.^5
∂ₓₓu = D₂(u,n,Δx,c,order=6)
@test_broken norm(∂ₓₓuₑ[7:n-6] .- ∂ₓₓu[7:n-6]) ≤ 1.0e-14
=#
=#