
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
function ParallelPenalty1D!(interp::Function,u::AbstractArray{TT},u₀::AbstractArray{TT},planes::ParallelGrid{TT,1},Δt::TT,grid::Grid1D,τ::TT,κ::TT,H::AbstractArray{TT}) where TT
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
function ParallelPenalty2D!(u::AbstractArray{TT},u₀::AbstractArray{TT},Δt::TT,
        PGrid::ParallelGrid{TT,2},D::Grid2D,
        nx::Integer,ny::Integer,τ::TT,κ::TT,H_x::AbstractArray{TT},H_y::AbstractArray{TT},perp::TT) where TT
    local wf :: TT
    local wb :: TT
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



"""
    For multiblock
"""
function applyParallelPenalty! end

function applyParallelPenalty!(u::VT,u₀::VT,Δt::TT,θ::TT,P::ParallelData{TT,1},grid::Grid1D{TT,MET}) where {TT,VT,MET}
    κ = P.κ
    w_f = P.w_f
    H = P.H

    I = LinearInterpolator(P.gridx,u)

    Fplane = P.PGrid.Fplane
    Bplane = P.PGrid.Bplane

    # w_f ← P_f u + P_b u
    for i in eachindex(P.gridx)
        x_f = Fplane.x[i]
        x_b = Bplane.x[i]
        if grid[1] ≤ x_f ≤ grid[grid.n]
            w_f[i] = I(x_f)
        else
            w_f[i] = TT(0)
        end
        if grid[1] ≤ x_b ≤ grid[grid.n]
            w_f[i] += I(x_b)
        else
            w_f[i] += TT(0)
        end
        # w_f[i] = w_fi + w_bi
        # w_f[i] = w_f[i]/2
    end

    # Tune the parallel penalty
    τ = P.τ/grid.Δx
    # P.τ_i[1] = τ * κ * 1/(maximum(J) * P.Δx*P.Δy * maximum(H.H[1].Boundary) * maximum(H.H[2].Boundary))

    P.w_b[1,1] = maximum(u - w_f)

    # u^{n+1} = (1+θ)Δt κ τ
    # _compute_u!(u,w_f,κ,τ,Δt,H,length(P.gridx),length(P.gridy))

    for i in 1:grid.n
        u[i] = TT(1)/(TT(1) - τ*κ/TT(2) * Δt * H.H[1][i]) *
            (u[i] - Δt * κ * τ / TT(4) * H.H[1][i] * w_f[i])
    end

    u
end

function applyParallelPenalty!(u::AbstractArray{TT},u₀::AbstractArray{TT},Δt::TT,θ::TT,
    P::ParallelData{TT,2},grid::Grid2D{TT,MET}) where {TT,MET}

    κ = P.κ
    w_f = P.w_f
    H = P.H
    J = grid.J

    I   = BicubicInterpolator(P.gridx,P.gridy,u)

    # w_f ← P_f u + P_b u
    _compute_w!(I,w_f,P.PGrid.Fplane,P.PGrid.Bplane,length(P.gridx),length(P.gridy))
    
    # Tune the parallel penalty
    τ = P.τ*0.1*(maximum(abs.(u - w_f))/ maximum(abs.(w_f)))^2.0
    P.τ_i[1] = τ * κ * 1/(maximum(J) * P.Δx*P.Δy * maximum(H.H[1].Boundary) * maximum(H.H[2].Boundary))

    P.w_b[1,1] = maximum(u - w_f)

    # u^{n+1} = (1+θ)Δt κ τ
    if MET == CartesianMetric
        _compute_u!(u,w_f,κ,τ,Δt,H,length(P.gridx),length(P.gridy))
    else
        _compute_u_curvilinear!(u,w_f,κ,τ,Δt,H,length(P.gridx),length(P.gridy),J)
    end

    # Replace the old parallel diffusion with the new one
    # TODO : FIX THIS FOR ADAPTIVE
    # @. w_b = u - w_f

    # u = w_f
    u
end

function _compute_u!(u::AT,w_f::AT,κ::TT,τ::TT,Δt::TT,H::CompositeH,nx::Int,ny::Int) where {TT,AT}
    # println(norm(u))
    # um = maximum(abs.(u .- w_f))
    for j in 1:ny
        for i in 1:nx
            ### CURRENTLY IN THE PAPER

            u[i,j] = 1/(1 + Δt * κ * τ / (H[i,j])) * (
                u[i,j] + 
                +Δt*κ*τ*w_f[i,j]/(H[i,j])
                    )

        end
    end
    u
end

function _compute_u_curvilinear!(u::AT,w_f::AT,κ::TT,τ::TT,Δt::TT,H::CompositeH,nx::Int,ny::Int,J::AT) where {TT,AT}
    for j in 1:ny
        for i in 1:nx

            u[i,j] = 1/(1 + Δt * κ * τ / (J[i,j]*H[i,j])) * (
                u[i,j] + 
                +Δt*κ*τ*w_f[i,j]/(J[i,j]*H[i,j])
                )

        end
    end
    u
end

"""
    _compute_w!
Computes ``P_parallel u`` and stores it in `dest`.
"""
function _compute_w!(itp,dest::AT,Fplane::ParGrid,Bplane::ParGrid,nx::Int,ny::Int) where {AT}
    for j in 1:ny
        for i in 1:nx
            dest[i,j] = itp(Fplane.x[i,j],Fplane.y[i,j]) + itp(Bplane.x[i,j],Bplane.y[i,j])
            dest[i,j] = dest[i,j]/2
        end
    end
    dest
end
function _compute_w!(itp,dest::AT,plane::ParGrid,nx::Int,ny::Int) where {AT}
    for j in 1:ny
        for i in 1:nx
            dest[i,j] = itp(plane.x[i,j],plane.y[i,j])
        end
    end
    dest
end

"""
    compute_parallel_operator
"""
function compute_parallel_operator(dest::AT,u::AT,P::ParallelData) where {AT}
    # I = scale( interpolate(u,BSpline(Cubic())),(P.gridx,P.gridy) )
    # I = scale( interpolate(u,BSpline(Quadratic())),(P.gridx,P.gridy) )
    # I = scale( interpolate(u,BSpline(Linear())),(P.gridx,P.gridy) )
    I = BicubicInterpolator(P.gridx,P.gridy,u)
    # _compute_w!(I,dest,P.PGrid.Fplane,P.PGrid.Bplane,P.gridx.len,P.gridy.len)
    _compute_w!(I,dest,P.PGrid.Fplane,P.PGrid.Bplane,length(P.gridx),length(P.gridy))
    @. dest = u - dest
end


"""
    applyParallelPenalty!(u::AbstractArray{TT},uglobal::Vector{Matrix{TT}},Δt::TT,P::ParallelData,grid::Grid2D{TT,MET}) where {TT,MET}
Applies the parallel penalty for a multiblock problem.
"""
function applyParallelPenalty!(u::AbstractArray{TT},uglobal::Vector{Matrix{TT}},Δt::TT,P::ParallelData,grid::Grid2D{TT,MET}) where {TT,MET}
    
    # Ipt = [BicubicInterpolator(P.gridx,P.PGrid.Fplane.y,uglobal[I]) for I in eachindex(uglobal)]

    κ = P.κ
    w_f = P.w_f
    H = P.H
    J = grid.J

    sgiF = P.PGrid.Fplane.subgrid
    sgiB = P.PGrid.Bplane.subgrid
    nnF = P.PGrid.Fplane.x
    nnB = P.PGrid.Bplane.x

    for I in eachindex(w_f)
        w_f[I] = uglobal[sgiF[I]][nnF[I]]
        w_f[I] += uglobal[sgiB[I]][nnB[I]]
        w_f[I] = w_f[I]/2
    end

    τ = P.τ * 0.1 * (maximum(abs.(u - w_f))/ maximum(abs.(w_f)))^2.0

    for j in 1:grid.ny
        for i in 1:grid.nx
            u[i,j] = 1/(1 + Δt * κ * τ / H[i,j]) * (
                u[i,j] + Δt * κ * τ * w_f[i,j] / H[i,j]
            )
        end
    end

end
