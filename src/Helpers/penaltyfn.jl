





function generate_parallel_penalty(w_f::Function,w_b::Function,grid::Grid2D,order::Int;κ::T=1,τ::T=-1,interpmode::Symbol=:linear) where T

    H_x = SBP_operators.build_H(ny,order)
    H_y = SBP_operators.build_H(nx,order)

    H_x = 1.0 ./H_x.^2
    H_y = 1.0 ./H_y.^2

    let H_x=H_x, H_y=H_y,
            w_f=w_f, w_b=w_b,
            grid=grid,
            τ=τ, κ=κ

        ParPen(u,u₀,Δt) = ParallelPenalty(w_f,w_b,interp,u,u₀,Δt,grid,τ,κ,H_x,H_y)
    end


end


struct parallel_storage
    F_plane
    B_plane
end

function ParallelPenalty(w_f::Function,w_b::Function,interp::Function,u::AbstractArray{T},u₀::AbstractArray{T},Δt::T,grid::Grid2D,τ::T,κ::T,H_x::AbstractArray{T},H_y::AbstractArray{T}) where T
    u_f = w_f(grid)
    u_b = w_b(grid)
    for j = 1:ny
        for i = 1:nx
            u[i,j] = 1.0/(1.0 - κ/2.0 * Δt * (H_y[i] + H_x[j])) * (u₀[i,j] - Δt*τ/4.0 * (H_y[i]+H_x[j])*(u_f[i,j] + u_b[i,j]))
            # u[i,j] = (0.5Δt*τ*(H_y[i]+H_x[j])*(u_f[i,j]+u_b[i,j]) - 2u[i,j])/(Δt*κ*(H_y[i]+H_x[j]) - 2.0)
        end
    end
end


# function forward_map(interp,foward_points)
# end

# function backward_map(interp,backward_points)
# end

