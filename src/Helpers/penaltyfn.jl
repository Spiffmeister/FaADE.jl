


"""
parallel_storage{T,N}
"""

# struct parallel_storage{T,N}
#     F_plane :: AbstractArray{T}
#     B_plane :: AbstractArray{T}
# end

# function parallel_storage()
#     if typeof(grid) == Grid1D
#         F = zeros(T,grid.n)
#         B = zeros(T,grid.n)
#         new(F,B)

#     elseif typeof(grid) == Grid2D
#     end
# end


"""
    parallel_grid
Storage for parallel grids
"""
struct ParallelGrid{T,N}
    z               :: Vector{T}
    nz              :: Int
    FowardPlane     :: AbstractArray{T} #Need x and y for 2D grid
    BackwardPlane   :: AbstractArray{T}

    function ParallelGrid(ForwardGrids,BackwardGrids,z::T) where T
        if typeof(z) == Float64
            zplanes = [z]
        elseif typeof(z) == AbstractArray
            zplanes = z
        end

        if ndims(ForwardGrids) == 1
            new{T,1}(zplanes,length(zplanes),ForwardGrids,BackwardGrids)
        elseif ndims(ForwardGrids) == 2
            new{T,2}(zplanes,length(zplanes),ForwardGrids,BackwardGrids)
        end
    end
end


struct PGrid1D{T} <: ParallelGridStorage{T,1}
    z               :: Vector{T}
    nz              :: Int
    ForwardPlane    :: AbstractArray{T}
    BackwardPlane   :: AbstractArray{T}
end

struct PGrid2D{T} <: ParallelGridStorage{T,2}
    z               :: Vector{T}
    nz              :: Int
    ForwardPlaneX   :: AbstractArray{T}
    ForwardPlaneY   :: AbstractArray{T}
    BackwardPlaneX  :: AbstractArray{T}
    BackwardPlaneY  :: AbstractArray{T}
end




# struct P∥{T,N}

# end



# function choose_interpmode(;interpmode::Symbol=:linear)
#     if interpmode == :linear
#         return Interpolations.LinearInterpolation
#     else
#         error("interpmode must be :linear")
#     end
# end

function generate_parallel_penalty(planes::ParallelGrid,grid::Grid1D,order::Int;κ::T=1.0,τ::T=-1.0,interpmode::Symbol=:cust,interpfn::Union{Nothing,Function}=nothing) where T
    # interp= choose_interpmode(interpmode=interpmode)
    if typeof(interpfn) == Nothing
        error("No interpolation function specified.")
    end

    H = build_H(order,grid.n)
    H = H.^-1.0
    
    α = -2.0
    τ = α/2.0

    let H=H, grid=grid, τ=τ, κ=κ, planes=planes

        ParPen(u,u₀,Δt) = ParallelPenalty1D!(interpfn,u,u₀,planes,Δt,grid,τ,κ,H)
    end
end

function generate_parallel_penalty(planes::ParallelGrid,grid::Grid2D,order::Int;κ::T=1.0,τ::T=-1.0,interpmode::Symbol=:linear,interpfn::Union{Nothing,Function}=nothing) where T

    # interp = choose_interpmode(interpmode=interpmode)
    if typeof(interpfn) == Nothing
        error("No interpolation function specified.")
    end

    H_x = SBP_operators.build_H(grid.ny,order)
    H_y = SBP_operators.build_H(grid.nx,order)

    H_x = 1.0 ./H_x.^2
    H_y = 1.0 ./H_y.^2

    # α = -2.0
    # τ = α/2.0 #TODO MIGHT BE A SQUARED HERE

    let H_x=H_x, H_y=H_y,
            grid=grid,
            τ=τ, κ=κ

        ParPen(u,u₀,Δt) = ParallelPenalty(interp,u,u₀,Δt,grid,τ,κ,H_x,H_y)
    end
    return ParPen

end



function ParallelPenalty1D!(interp::Function,u::AbstractArray{T},u₀::AbstractArray{T},planes::ParallelGrid{T,1},Δt::T,grid::Grid1D,τ::T,κ::T,H::AbstractArray{T}) where T
    I = interp(grid.grid,u₀)
    for i = 1:grid.n
        u[i] = 1.0/(1.0 - κ*τ/2.0 * Δt * H[i]) * 
            (u[i] - Δt*κ*τ/4.0 * H[i] * (I(planes.FowardPlane[i]) + I(planes.BackwardPlane[i])))
        
    end
end

function ParallelPenalty2D!(interp::Function,u::AbstractArray{T},u₀::AbstractArray{T},Δt::T,grid::Grid2D,τ::T,κ::T,H_x::AbstractArray{T},H_y::AbstractArray{T}) where T

    I = interp((grid.gridx,grid.gridy),u₀)

    for j = 1:ny
        for i = 1:nx
            u[i,j] = 1.0/(1.0 - κ/2.0 * Δt * (H_y[i] + H_x[j])) * (u[i,j] - Δt*τ/4.0 * (H_y[i]+H_x[j])*(I(F_plane[i,j]) + I(B_plane([i,j]))))
        end
    end
end




# function forward_map(interp,foward_points)
# end

# function backward_map(interp,backward_points)
# end

