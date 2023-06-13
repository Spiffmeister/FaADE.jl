
"""
    generate_parallel_penalty
Generates a 2D or 3D parallel penalty function given a parallel grid mapping.

Operator in the 1D case is found [here](@ref https://arxiv.org/abs/2303.15447). In the 2D case details can be found [here](@ref https://arxiv.org/abs/2306.00423).
"""
function generate_parallel_penalty end
function generate_parallel_penalty(planes::ParallelGrid{T,1,AT},order::Int;κ::T=1.0,τ::T=-1.0,interpmode::Symbol=:cust,interpfn::Union{Nothing,Function}=nothing) where {T,AT}
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
function generate_parallel_penalty(pargrid::ParallelGrid{T,2,AT},grid::Grid2D,order::Int;κ::T=1.0,τ::T=-1.0,perp::T=1.0,interpfn::Union{Nothing,Function}=nothing) where {T,AT}

    # interp = choose_interpmode(interpmode=interpmode)
    # if typeof(interpfn) == Nothing
    #     interp = LinearInterpolation
    # end

    H_x = build_H(order,grid.ny)
    H_y = build_H(order,grid.nx)

    H_x = 1.0 ./(H_x*grid.Δx)
    H_y = 1.0 ./(H_y*grid.Δy)

    # α = -2.0
    # τ = α/2.0 #TODO MIGHT BE A SQUARED HERE

    # itp(u) = LinearInterpolation((grid.gridx,grid.gridy),u)

    let H_x=H_x, H_y=H_y,
            τ=τ, κ=κ,
            nx=grid.nx, ny=grid.ny,
            D=grid,
            PGrid = pargrid,
            perp = perp
        ParPen(u,u₀,Δt) = ParallelPenalty2D!(u,u₀,Δt,PGrid,D,nx,ny,τ,κ,H_x,H_y,perp)
        # ParPen(u,u₀,Δt) = ParallelPenalty2D!(interp,u,u₀,Δt,PGrid,D,nx,ny,τ,κ,H_x,H_y)
        return ParPen
    end

end


"""
    ParallelPenalty1D
The default generated parallel penalty function for 2D (1D+parallel) problems
"""
function ParallelPenalty1D!(interp::Function,u::AbstractArray{T},u₀::AbstractArray{T},planes::ParallelGrid{1},Δt::T,grid::Grid1D,τ::T,κ::T,H::AbstractArray{T}) where T
    I = interp(grid.grid,u₀)
    for i = 1:grid.n
        u[i] = 1.0/(1.0 - κ*τ/2.0 * Δt * H[i]) * 
            (u[i] - Δt*κ*τ/4.0 * H[i] * (I(planes.FowardPlane[i]) + I(planes.BackwardPlane[i])))
    end
end
"""
    ParallelPenalty2D!(u::AbstractArray{T},u₀::AbstractArray{T},Δt::T,pargrid::ParallelGrid{2},grid::Grid2D,τ::T,κ::T,H_x::AbstractArray{T},H_y::AbstractArray{T}) where T
Default generated parallel penalty function for 3D (2D+parallel) problems
"""
function ParallelPenalty2D!(u::AbstractArray{T},u₀::AbstractArray{T},Δt::T,
        PGrid::ParallelGrid,D::Grid2D,
        nx::Integer,ny::Integer,τ::T,κ::T,H_x::AbstractArray{T},H_y::AbstractArray{T},perp::T) where T
    local wf :: T
    local wb :: T
    I = linear_interpolation((D.gridx,D.gridy),u)
    # I = scale(interpolate(u,BSpline(Linear())),(D.gridx,D.gridy))
    # I = scale(interpolate(u, BSpline(Quadratic(Line(OnGrid()))),),
        # range(D.gridx[1],D.gridx[end],length=D.nx),range(D.gridy[1],D.gridy[end],length=D.ny))
    # I = LinearInterpolation((D.gridx,D.gridy),u₀)

    # τ = -(1.e-14 + 1.0/κ)
    # τ = -(perp/κ)
    # τ = -sqrt(D.Δx * D.Δy/( (D.gridx[end]-D.gridx[1]) * (D.gridy[end] - D.gridy[1]) ))
    τ = -sqrt((D.gridx[end]-D.gridx[1])*(D.gridy[end]-D.gridy[1])/(D.Δx*D.Δy))  
    # -1.0/sqrt(D.Δx * D.Δy/( (D.gridx[end]-D.gridx[1]) * (D.gridy[end] - D.gridy[1]) ))
    # τ = -T(1)
    # τ = -D.Δx
    
    for j = 1:ny
        for i = 1:nx
            # H = H_y[i]*H_x[j]
            # H = 1.0
            wf = I(PGrid.Fplane.x[i,j],PGrid.Fplane.y[i,j])
            wb = I(PGrid.Bplane.x[i,j],PGrid.Bplane.y[i,j])
            u[i,j] = 1.0/(1.0 - κ* τ * Δt) * 
                ( u[i,j] -  κ*Δt*τ/2.0 * (wf + wb) )
            
            # u[i,j] = 1.0/(1.0 - κ* τ/2.0 * Δt * H) * 
            #     ( u[i,j] -  κ*Δt*τ/4.0 * H * 
            #         (I(PGrid.Fplane.x[i,j],PGrid.Fplane.y[i,j]) + I(PGrid.Bplane.x[i,j],PGrid.Bplane.y[i,j])) )
    
            
            # u[i,j] = u[i,j] + Δt * κ * τ/2.0 * H * (u[i,j] - 0.5*(I(PGrid.Fplane.x[i,j],PGrid.Fplane.y[i,j]) + I(PGrid.Bplane.x[i,j],PGrid.Bplane.y[i,j])))
            



        end
    end
end


