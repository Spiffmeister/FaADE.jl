using Test
using FaADE



##======##
# Grid1D
##======##
@testset "1D grid tests" begin

    𝒟 = [0.0,1.0]
    n = 5
    grid = LinRange(𝒟[1], 𝒟[2], n)
    Δx = (𝒟[2] - 𝒟[1]) / (n - 1)
    g = FaADE.Grid1D(collect(grid), Δx, n)
    @test g.grid ≈ grid rtol=1e-6
    @test g.Δx ≈ Δx rtol=1e-6
    @test g.n == n

end



# Test Grid2D struct with CartesianMetric
function test_Grid2D_CartesianMetric()
    𝒟x = [0.0, 1.0]
    𝒟y = [0.0, 2.0]
    nx, ny = 5, 3
    grid = (linspace(𝒟x[1], 𝒟x[2], nx), linspace(𝒟y[1], 𝒟y[2], ny))
    Δx, Δy = (𝒟x[2] - 𝒟x[1]) / (nx - 1), (𝒟y[2] - 𝒟y[1]) / (ny - 1)
    g = Grid2D{Float64, CartesianMetric, Float64, Tuple{Vector{Float64}, Vector{Float64}}}(grid, Δx, Δy, nx, ny)
    @test g.gridx ≈ grid[1] rtol=1e-6
    @test g.gridy ≈ grid[2] rtol=1e-6
    @test g.Δx ≈ Δx rtol=1e-6
    @test g.Δy ≈ Δy rtol=1e-6
    @test g.nx == nx
    @test g.ny == ny
end

# Run all tests
@testset "grid.jl tests" begin
    @testset "Grid2D tests" begin
        test_Grid2D_CartesianMetric()
    end
end




# Test Grid2D struct with curvilinear coordinates
function test_Grid2D_Curvilinear()
    𝒟x = [0.0, 1.0]
    𝒟y = [0.0, 2.0]
    nx, ny = 5, 3
    grid = (linspace(𝒟x[1], 𝒟x[2], nx), linspace(𝒟y[1], 𝒟y[2], ny))
    Δx, Δy = (𝒟x[2] - 𝒟x[1]) / (nx - 1), (𝒟y[2] - 𝒟y[1]) / (ny - 1)
    metric = CurvilinearMetric((x, y) -> (x^2, y^2), (x, y) -> (2x, 2y))
    g = Grid2D{Float64, CurvilinearMetric, Float64, Tuple{Vector{Float64}, Vector{Float64}}}(grid, Δx, Δy, nx, ny, metric)
    @test g.gridx ≈ grid[1] .^ 2 rtol=1e-6
    @test g.gridy ≈ grid[2] .^ 2 rtol=1e-6
    @test g.Δx ≈ Δx rtol=1e-6
    @test g.Δy ≈ Δy rtol=1e-6
    @test g.nx == nx
    @test g.ny == ny
end

# Run all tests
@testset "grid.jl tests" begin
    @testset "Grid2D tests" begin
        test_Grid2D_Curvilinear()
    end
end