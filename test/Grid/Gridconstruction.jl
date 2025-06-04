using Test
using FaADE



##======##
# Grid1D
##======##
@testset "1D grid tests" begin

    𝒟 = [0.0,1.0]
    n = 5
    grid = range(𝒟[1], 𝒟[2], n)
    Δx = (𝒟[2] - 𝒟[1]) / (n - 1)
    g = Grid1D(𝒟, n)
    @test g.grid ≈ grid rtol=1e-12
    @test g.Δx ≈ Δx rtol=1e-12
    @test g.n == n

end



# Test Grid2D struct with CartesianMetric
@testset "2D cartesian grid construction" begin
    𝒟x = [0.0, 1.0]
    𝒟y = [0.0, 2.0]
    nx, ny = 5, 3
    gridx = range(𝒟x[1], 𝒟x[2], nx)
    gridy = range(𝒟y[1], 𝒟y[2], ny)
    Δx, Δy = (𝒟x[2] - 𝒟x[1]) / (nx - 1), (𝒟y[2] - 𝒟y[1]) / (ny - 1)
    g = Grid2D(𝒟x,𝒟y, nx, ny)
    @testset "gridx" begin
        @test all(g.gridx[:,1] .== gridx)
        @test all(g.gridx[:,2] .== gridx)
        @test all(g.gridx[:,3] .== gridx)
    end
    @testset "gridy" begin
        @test all(g.gridy[1,:] .== gridy)
        @test all(g.gridy[2,:] .== gridy)
        @test all(g.gridy[3,:] .== gridy)
        @test all(g.gridy[4,:] .== gridy)
        @test all(g.gridy[5,:] .== gridy)
    end
    @test g.Δx ≈ Δx
    @test g.Δy ≈ Δy
    @test g.nx == nx
    @test g.ny == ny
end




# Test Grid2D struct with curvilinear coordinates
@testset "2D curvilinear grid construction" begin

    nx = 5
    ny = 3

    cbottom(u)  = u * [1.0, 0.0]
    cleft(v)    = v * [0.0, 1.0]
    cright(v)   = v * [0.0, 1.0] + [1.0, 0.0]
    ctop(u)     = u * [1.0, 0.0] + [0.0, 1.0]

    Dom = Grid2D(cbottom, cleft, cright, ctop, nx, ny)

    gridx = range(0.0,1.0,nx)
    gridy = range(0.0,1.0,ny)

    Δx, Δy = 1/(nx - 1), 1/(ny - 1)

    g = Grid2D(cbottom,cleft,cright,ctop, nx, ny)
    @testset "gridx" begin
        @test all(g.gridx[:,1] .== gridx)
        @test all(g.gridx[:,2] .== gridx)
        @test all(g.gridx[:,3] .== gridx)
    end
    @testset "gridy" begin
        @test all(g.gridy[1,:] .== gridy)
        @test all(g.gridy[2,:] .== gridy)
        @test all(g.gridy[3,:] .== gridy)
        @test all(g.gridy[4,:] .== gridy)
        @test all(g.gridy[5,:] .== gridy)
    end
    @test all(g.J .≈ 1.0)
    @test all(g.qx .≈ 1.0)
    @test all(g.qy .≈ 0.0)
    @test all(g.rx .≈ 0.0)
    @test all(g.ry .≈ 1.0)
    @test g.Δx ≈ Δx
    @test g.Δy ≈ Δy
    @test g.nx == nx
    @test g.ny == ny
end

