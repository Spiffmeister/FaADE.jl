"""
    applyParallelPenalty!
The default generated parallel penalty function for 2D (1D+parallel) problems

Operator in the 1D case is found [here](@ref https://arxiv.org/abs/2303.15447). In the 2D case details can be found [here](@ref https://arxiv.org/abs/2306.00423).    
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
# function applyParallelPenalty!(u::AbstractArray{TT},uglobal::Vector{Matrix{TT}},Δt::TT,P::ParallelData,grid::Grid2D{TT,MET}) where {TT,MET}
function applyParallelPenalty!(u::AbstractArray{TT},uglobal::Vector{Matrix{TT}},Δt::TT,P::Vector{ParallelData{TT,DIM,PGT,GT,BT,IT}},grid::Grid2D{TT,MET},I) where {TT,MET,DIM,PGT,GT,BT,IT}
    
    Ipt = [BicubicInterpolator(P[I].gridx,P[I].gridy,uglobal[I]) for I in eachindex(uglobal)]
    # Ipt = [interpolate(P[I].gridx,P[I].gridy,uglobal[I]) for I in eachindex(uglobal)]

    κ = P[I].κ
    w_f = P[I].w_f
    H = P[I].H
    J = grid.J

    sgiF = P[I].PGrid.Fplane.subgrid
    sgiB = P[I].PGrid.Bplane.subgrid
    # nnF = P.PGrid.Fplane.x
    # nnB = P.PGrid.Bplane.x

    nnFx = P[I].PGrid.Fplane.x
    nnFy = P[I].PGrid.Fplane.y
    nnBx = P[I].PGrid.Bplane.x
    nnBy = P[I].PGrid.Bplane.y

    for J in eachindex(w_f)
        # w_f[I] = uglobal[sgiF[I]][nnF[I]]
        # w_f[I] += uglobal[sgiB[I]][nnB[I]]
        w_f[J] = Ipt[sgiF[J]](nnFx[J],nnFy[J])
        w_f[J] += Ipt[sgiB[J]](nnBx[J],nnBy[J])
        w_f[J] = w_f[J]/2
    end

    τ = P[I].τ * 0.1 * (maximum(abs.(u - w_f))/ maximum(abs.(w_f)))^2.0

    for j in 1:grid.ny
        for i in 1:grid.nx
            u[i,j] = 1/(1 + Δt * κ * τ / H[i,j]) * (
                u[i,j] + Δt * κ * τ * w_f[i,j] / H[i,j]
            )
        end
    end

end
function applyParallelPenalty!(u::AbstractArray{TT},uglobal::Vector{Matrix{TT}},Δt::TT,P::Vector{ParallelData{TT,DIM,ParallelGrid{TT,DIM,PMT,AT},GT,BT,IT}},grid::Grid2D{TT,MET},I) where {TT,MET,DIM,AT,GT,BT,IT, PMT<:ParGridLinear}
    
    κ = P[I].κ
    w_f = P[I].w_f
    H = P[I].H
    J = grid.J

    sgiF = P[I].PGrid.Fplane.subgrid
    sgiB = P[I].PGrid.Bplane.subgrid


    nnFx = P[I].PGrid.Fplane.x
    nnFy = P[I].PGrid.Fplane.y
    nnBx = P[I].PGrid.Bplane.x
    nnBy = P[I].PGrid.Bplane.y

    for J in eachindex(w_f)
        # w_f[I] = uglobal[sgiF[I]][nnF[I]]
        # w_f[I] += uglobal[sgiB[I]][nnB[I]]
        w_f[J] = Ipt[sgiF[J]](nnFx[J],nnFy[J])
        w_f[J] += Ipt[sgiB[J]](nnBx[J],nnBy[J])
        w_f[J] = w_f[J]/2
    end

    τ = P[I].τ * 0.1 * (maximum(abs.(u - w_f))/ maximum(abs.(w_f)))^2.0

    for j in 1:grid.ny
        for i in 1:grid.nx
            u[i,j] = 1/(1 + Δt * κ * τ / H[i,j]) * (
                u[i,j] + Δt * κ * τ * w_f[i,j] / H[i,j]
            )
        end
    end

end



function _linear_interpolation(u::AbstractArray{TT},wx::TT,wy::TT,i::Int,j::Int) where TT

    w = (u[i,j] * w11) + (u[i+1,j] * w12) + (u[i,j+1] * w21) + (u[i+1,j+1] * w22)

    return w
end


