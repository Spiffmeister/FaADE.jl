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
"""
Parallel penalty for single block problems
"""
function applyParallelPenalty!(u::AbstractArray{TT},t::TT,Δt::TT,
    P::ParallelData{TT,2,PGT,GT,BT,IT},grid::Grid2D{TT,MET}) where {TT,MET,PGT,GT,BT,IT}

    κ = P.κ
    w_f = P.w_f
    H = P.H
    J = grid.J

    if IT <: Nothing
        I   = BicubicInterpolator(P.gridx,P.gridy,u)
    elseif IT <: BivariateCHSInterpolation
        I = P.Interpolant
        IC = P.Intercept
    end

    # w_f ← P_f u + P_b u
    if isnothing(IC)
        _compute_w!(I,w_f,P.PGrid.Fplane,P.PGrid.Bplane,length(P.gridx),length(P.gridy))
    else
        _compute_w!(I,IC,w_f,P.PGrid.Fplane,P.PGrid.Bplane,t,length(P.gridx),length(P.gridy))
    end
    
    # Tune the parallel penalty
    τ = P.τ*0.1*(maximum(abs.(u - w_f))/ maximum(abs.(w_f)))^2.0
    isinf(τ) ? τ = TT(1) : nothing
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
"""
Applies the parallel penalty for a multiblock problem.
"""
function applyParallelPenalty!(u::AbstractArray{TT},τ::TT,Δt::TT,P::Dict{Int64,ParallelData{TT,DIM}},grid::Grid2D{TT,MET},I) where {TT,MET,DIM}

    κ = P[I].κ
    w_f = P[I].w_f
    H = P[I].H
    Jac = grid.J

    if length(Jac) > 1
        for j in 1:grid.ny
            for i in 1:grid.nx
                u[i,j] = 1/(1 + Δt * κ * τ / (Jac[i,j] * H[i,j])) * (
                    u[i,j] + Δt * κ * τ * w_f[i,j] / (Jac[i,j] * H[i,j])
                )
            end
        end
    else
        for j in 1:grid.ny
            for i in 1:grid.nx
                u[i,j] = 1/(1 + Δt * κ * τ / H[i,j]) * (
                    u[i,j] + Δt * κ * τ * w_f[i,j] / H[i,j]
                )
            end
        end
    end

    # τ = P[I].τ * 0.1 * (maximum(abs.(u - w_f))/ maximum(abs.(w_f)))^2.0

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
function _compute_w! end
"""
Compute `w = (w_f + w_b)/2` when an interpolating function is provided
"""
function _compute_w!(itp,dest::AT,Fplane::ParGrid,Bplane::ParGrid,nx::Int,ny::Int) where {AT}
    # @show size(dest), size(Fplane.x)
    for j in 1:ny
        for i in 1:nx
            dest[i,j] = itp(Fplane.x[i,j],Fplane.y[i,j]) + itp(Bplane.x[i,j],Bplane.y[i,j])
            # if isnan(dest[i,j])
            #     @show i,j, (Fplane.x[i,j], Fplane.y[i,j]), (Bplane.x[i,j], Bplane.y[i,j])
            #     dest[i,j] = zero(eltype(dest))
            # end
            dest[i,j] = dest[i,j]/2
        end
    end
    dest
end
"""
If an `intercept(x,y,t)` function is provided
"""
function _compute_w!(itp,itc,dest::AT,Fplane::ParGrid,Bplane::ParGrid,t::TT,nx::Int,ny::Int) where {AT,TT}
    for j in 1:ny
        for i in 1:nx
            dest[i,j] = itp(Fplane.x[i,j],Fplane.y[i,j]) + itp(Bplane.x[i,j],Bplane.y[i,j])
            # if isnan(dest[i,j])
                dest[i,j] = itc(dest[i,j],Fplane.x[i,j],Fplane.y[i,j],t) + itc(dest[i,j],Bplane.x[i,j],Bplane.y[i,j],t)
            # end
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
function _compute_w!(itp,dest::AT,plane::ParGrid) where {AT}
    for j in 1:plane.y.len
        for i in 1:plane.x.len
            dest[i,j] = itp(plane.x[i,j],plane.y[i,j])
        end
    end
    dest
end
# function _compute_w!(itp::BivariateCHSInterpolation,dest::AT,FPlane::ParGrid,BPlane::ParGrid,nx::Int,ny::Int) where AT
#     for j in 1:ny
#         for i in 1:nx
#             dest[i,j] = itp()
#         end
#     end
# end

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
    computeglobalw!
"""
function computeglobalw! end

function computeglobalw!(PD::ParallelMultiBlock{TT,DIM,IT,CT},uglobal::VAT,t::TT,Δt::TT) where {TT,VAT,DIM,IT,CT}

    Interpolant = PD.Interpolant
    Intercept = PD.Intercept
    τglobal = PD.τ

    for I in eachindex(Interpolant)
        w_f = PD.PData[I].w_f
        PData = PD.PData[I].PGrid
        computeglobalw!(Interpolant,Intercept,w_f,uglobal[I],PData,τglobal,t,I)
    end
    
    τglobal .= τglobal * PD.PData[1].τ
end
function computeglobalw!(interpolant::IT,intercept::CT,w_f::AT,u::AT,PGrid::ParallelGrid,τglobal::VT,t::TT,I::Int) where {AT,VT,TT,IT,CT}

    sgiF = PGrid.Fplane.subgrid
    sgiB = PGrid.Bplane.subgrid

    nnFx = PGrid.Fplane.x
    nnFy = PGrid.Fplane.y
    nnBx = PGrid.Bplane.x
    nnBy = PGrid.Bplane.y

    for J in eachindex(w_f)
        tmpf = interpolant[sgiF[J]](nnFx[J],nnFy[J])
        tmpb = interpolant[sgiB[J]](nnBx[J],nnBy[J])
        # if isnan(tmpf)
        #     @show "aaaaaa"
        #     @show nnFx[J], nnFy[J], sgiF[J], J
        #     @show nnBx[J], nnBy[J], sgiB[J]
        #     @show tmpf, tmpb
        #     # @show interpolant[sgiF[J]].x
        #     error("aa")
        # end
        if !(CT == Nothing)
            tmpf = intercept[sgiF[J]](tmpf,nnFx[J],nnFy[J],t)
            tmpb = intercept[sgiB[J]](tmpb,nnBx[J],nnBy[J],t)
        end

        w_f[J] = (tmpf + tmpb)/2
    end
    τglobal[I] = 0.1 * (maximum(abs.(u - w_f))/ maximum(abs.(w_f)))^2.0

end


"""
Compute the globalw for parallel blocks where a package is used for bicubic interpolation
"""
function computeglobalw!(u::AbstractArray{TT},uglobal::Vector{Matrix{TT}},τglobal::Vector{TT},Δt::TT,P::Vector{ParallelData{TT,DIM,ParallelGrid{TT,DIM,PMT,AT},GT,BT,IT}},grid::Grid2D{TT,MET},I) where {TT,MET,DIM,AT,GT,BT,IT, PMT<:ParGrid}
    
    Ipt = [BicubicInterpolator(P[I].gridx,P[I].gridy,uglobal[I]) for I in eachindex(uglobal)]

    w_f = P[I].w_f
    # J = grid.J

    sgiF = P[I].PGrid.Fplane.subgrid
    sgiB = P[I].PGrid.Bplane.subgrid

    nnFx = P[I].PGrid.Fplane.x
    nnFy = P[I].PGrid.Fplane.y
    nnBx = P[I].PGrid.Bplane.x
    nnBy = P[I].PGrid.Bplane.y

    for J in eachindex(w_f)
        w_f[J] = Ipt[sgiF[J]](nnFx[J],nnFy[J])
        w_f[J] += Ipt[sgiB[J]](nnBx[J],nnBy[J])
        w_f[J] = w_f[J]/2
    end

    τglobal[I] = P[I].τ * 0.1 * (maximum(abs.(u - w_f))/ maximum(abs.(w_f)))^2.0

end
"""
Compute the globalw for parallel blocks where a custom interpolation scheme is used
"""
function computeglobalw!(u::AbstractArray{TT},uglobal::Vector{Matrix{TT}},τglobal::Vector{TT},Δt::TT,P::Vector{ParallelData{TT,DIM,ParallelGrid{TT,DIM,PMT,AT},GT,BT,IT}},grid::Grid2D{TT,MET},I) where {TT,MET,DIM,AT,GT,BT,IT, PMT<:ParGridLinear{TT,AT,METHOD}} where METHOD
    
    w_f = P[I].w_f

    sgiF = P[I].PGrid.Fplane.subgrid
    sgiB = P[I].PGrid.Bplane.subgrid


    nnFx = P[I].PGrid.Fplane.i
    nnFy = P[I].PGrid.Fplane.j
    nnBx = P[I].PGrid.Bplane.i
    nnBy = P[I].PGrid.Bplane.j

    w11F = P[I].PGrid.Fplane.weight11
    w12F = P[I].PGrid.Fplane.weight12
    w21F = P[I].PGrid.Fplane.weight21
    w22F = P[I].PGrid.Fplane.weight22

    w11B = P[I].PGrid.Bplane.weight11
    w12B = P[I].PGrid.Bplane.weight12
    w21B = P[I].PGrid.Bplane.weight21
    w22B = P[I].PGrid.Bplane.weight22

    if METHOD == :NearestNeighbour
        for J in eachindex(u)
            w_f[J] = uglobal[sgiB[J]][nnBx[J]]
            w_f[J] += uglobal[sgiF[J]][nnFx[J]]
            w_f[J] = w_f[J]/2
        end
    elseif METHOD == :CHS
        
    else
        for J in eachindex(w_f)
            # w_f[I] = uglobal[sgiF[I]][nnF[I]]
            # w_f[I] += uglobal[sgiB[I]][nnB[I]]
            # w_f[J] = Ipt[sgiF[J]](nnFx[J],nnFy[J])
            # w_f[J] += Ipt[sgiB[J]](nnBx[J],nnBy[J])
            w_f[J] = _interpolation(uglobal[sgiB[J]], w11B, w12B, w21B, w22B, nnBx[J], nnBy[J], J)
            w_f[J] += _interpolation(uglobal[sgiF[J]],w11F, w12F, w21F, w22F, nnFx[J], nnFy[J], J)
            w_f[J] = w_f[J]/2
        end
    end
    
    τglobal[I] = P[I].τ * 0.1 * (maximum(abs.(u - w_f))/ maximum(abs.(w_f)))^2.0

end





"""
    _linear_interpolation
"""
# function _linear_interpolation(u::AbstractArray{TT},wx::TT,wy::TT,i::Int,j::Int) where TT

#     w = (u[i,j] * w11) + (u[i+1,j] * w12) + (u[i,j+1] * w21) + (u[i+1,j+1] * w22)

#     return w
# end

"""
    _interpolation
"""
function _interpolation(u::AT,w11::AT,w12::AT,w21::AT,w22::AT,i::Int,j::Int,I::Int) where {AT}
    w = (u[i,j] * w11[I])
    w += (u[i+1,j] * w12[I])
    w += (u[i,j+1] * w21[I])
    w += (u[i+1,j+1] * w22[I])
    return w
end