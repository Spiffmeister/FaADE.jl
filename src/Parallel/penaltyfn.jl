


function generate_parallel_penalty(planes::ParallelGrid{1},order::Int;κ::T=1.0,τ::T=-1.0,interpmode::Symbol=:cust,interpfn::Union{Nothing,Function}=nothing) where T
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

function generate_parallel_penalty(pargrid::ParallelGrid{2},grid::Grid2D,order::Int;κ::T=1.0,τ::T=-1.0,interpfn::Union{Nothing,Function}=nothing) where T

    # interp = choose_interpmode(interpmode=interpmode)
    if typeof(interpfn) == Nothing
        interp = LinearInterpolation
    end

    H_x = build_H(order,grid.ny)
    H_y = build_H(order,grid.nx)

    H_x = 1.0 ./H_x.^2
    H_y = 1.0 ./H_y.^2

    # α = -2.0
    # τ = α/2.0 #TODO MIGHT BE A SQUARED HERE

    let H_x=H_x, H_y=H_y,
            grid=grid,
            τ=τ, κ=κ,
            pargrid=pargrid
        ParPen(u,u₀,Δt) = ParallelPenalty2D!(interp,u,u₀,Δt,pargrid,grid,τ,κ,H_x,H_y)
        return ParPen
    end

end


"""
"""
function ParallelPenalty1D!(interp::Function,u::AbstractArray{T},u₀::AbstractArray{T},planes::ParallelGrid{1},Δt::T,grid::Grid1D,τ::T,κ::T,H::AbstractArray{T}) where T
    I = interp(grid.grid,u₀)
    for i = 1:grid.n
        u[i] = 1.0/(1.0 - κ*τ/2.0 * Δt * H[i]) * 
            (u[i] - Δt*κ*τ/4.0 * H[i] * (I(planes.FowardPlane[i]) + I(planes.BackwardPlane[i])))
    end
end
"""
    ParallelPenalty2D!(interp::Function,u::AbstractArray{T},u₀::AbstractArray{T},Δt::T,pargrid::ParallelGrid{2},grid::Grid2D,τ::T,κ::T,H_x::AbstractArray{T},H_y::AbstractArray{T}) where T
"""
function ParallelPenalty2D!(interp::Function,u::AbstractArray{T},u₀::AbstractArray{T},Δt::T,pargrid::ParallelGrid{2},grid::Grid2D,τ::T,κ::T,H_x::AbstractArray{T},H_y::AbstractArray{T}) where T

    I = interp((grid.gridx,grid.gridy),u₀)

    for j = 1:grid.ny
        for i = 1:grid.nx
            u[i,j] = 1.0/(1.0 - κ/2.0 * Δt * (H_y[i] + H_x[j])) * (u[i,j] - Δt*τ/4.0 * (H_y[i]+H_x[j])*(I(pargrid.Fplane[i,j]) + I(pargrid.Bplane[i,j])))
        end
    end
end


