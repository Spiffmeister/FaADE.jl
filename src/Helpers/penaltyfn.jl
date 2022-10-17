





function generate_parallel_penalty(w_f::Function,w_b::Function,order::Int;κ::T=1,τ::T=-1) where T

    H_x = SBP_operators.build_H(ny,order)
    H_y = SBP_operators.build_H(nx,order)

    H_x = 1.0 ./H_x.^2
    H_y = 1.0 ./H_y.^2




end


struct parallel_storage
    F_plane
    B_plane
end

function P_∥(u::AbstractArray{T},u₀::AbstractArray{T},Δt::T,w_f::Function,w_b::Function,nx::Int,ny::Int,τ::T,κ::T) where T
    u_f = w_f()
    u_b = w_b()
    for j = 1:ny
        for i = 1:nx
            u[i,j] = 1.0/(1.0 - κ/2.0 * Δt * (H_y[i] + H_x[j])) * (u₀[i,j] - Δt*τ/4.0 * (H_y[i]+H_x[j])*(u_f[i,j] + u_b[i,j]))
            # u[i,j] = (0.5Δt*τ*(H_y[i]+H_x[j])*(u_f[i,j]+u_b[i,j]) - 2u[i,j])/(Δt*κ*(H_y[i]+H_x[j]) - 2.0)
        end
    end
end