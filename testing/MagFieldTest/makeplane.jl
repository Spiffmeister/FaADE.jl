module makeplane
    

    
    struct z_plane
        # plane_data      :: Array{Float64}
        z       :: Float64
        x       :: Matrix{Float64}
        y       :: Matrix{Float64}
        B       :: Matrix{Float64}
        xproj    :: Matrix{Int64}
        yproj    :: Matrix{Int64}
    end

    struct grid_data
        #=
        Holding grid and interpolation data
        =#
        # Arrays for ψ and θ grid points
        x                   :: Array{Float64}
        y                   :: Array{Float64}
        # Grid size
        Δx                  :: Float64
        Δy                  :: Float64
        Δz                  :: Float64
        # Size of the grid
        nx                   :: Int64 #x
        ny                   :: Int64 #y
        nz                   :: Int64 #z
        # x and y on the z_{k+1} and z_{k-1} planes
        z_planes            :: Vector{z_plane}
    end


    function construct_grid(𝒟x::Vector{Float64},𝒟y::Vector{Float64},nx::Int64,ny::Int64,nz::Int64,H::Function,params;
        method=:nearestneighbours)
        #=
        Take ψ₀,θ₀ on ζ=0 plane and form backward and forward interpolation grid points
        ψ₀ and θ₀ can be non-equidistant arrays of grid of points
            
            Requires: the field line Hamiltonian χ -- see FieldLines.jl for example
        =#

        # Set step sizes
        Δx = (𝒟x[end]-𝒟x[1])/(nx-1)
        Δy = (𝒟y[end]-𝒟y[1])/(ny-1)
        # Build arrays
        x = collect(range(𝒟x[1],stop=𝒟x[2],length=nx))
        y = collect(range(𝒟y[1],stop=𝒟y[2],length=ny))

        # Perform field line tracing
        if nz == 1
            # Forward plane
            z = 2π
            planex,planey = trace(H,x,y,z,params,nx,ny)
            xproj,yproj = nearest_neighbours(x,y,planex,planey,nx,ny)
            zf_plane = z_plane(z,planex,planey,ones(nx,ny),xproj,yproj)
            # Backward plane
            z = -2π
            planex,planey = trace(H,x,y,z,params,nx,ny)
            xproj,yproj = nearest_neighbours(x,y,planex,planey,nx,ny)
            zb_plane = z_plane(z,planex,planey,ones(nx,ny),xproj,yproj)

            gdata = grid_data(x,y,Δx,Δy,2π,nx,ny,nz,[zf_plane,zb_plane])
        end

        return gdata
    end



    #= FIELD LINE TRACING =#

    function trace(H::Function,x::Vector{Float64},y::Vector{Float64},z::Float64,params::H_params,nx::Int64,ny::Int64)
        # Returns an array of [ψ,θ] points on the desired ζ plane
        plane = zeros(2,nx*ny)
        x₀ = [[θ,ψ] for θ in y for ψ in x]

        function prob_fn(prob,i,repeat)
            remake(prob,u0=x₀[i])
        end
        P = ODEProblem(H,x₀[1],(0.0,z),params)
        EP = EnsembleProblem(P,prob_func=prob_fn)

        # dispatch_size = floor(Int64,m*n/nworkers())
        dispatch_size = nx*ny/nworkers()
        isinteger(dispatch_size) ? dispatch_size += 1 : nothing

        sim = solve(EP,Tsit5(),EnsembleDistributed(),trajectories=ny*nx,batch_size=floor(Int64,ny*nx/nworkers()),save_on=false,save_end=true)

        for i = 1:length(sim.u)
            plane[:,i] = mod.(sim.u[i][2:-1:1,2],2π)
        end

        planex = zeros(Float64,nx,ny)
        planey = zeros(Float64,nx,ny)

        for i = 1:length(sim.u)
            planex[i] = mod.(sim.u[i][2,2],2π)
            planey[i] = mod.(sim.u[i][1,2],2π)
        end

        return planex,planey
    end


    #= Weighting functions =#

    function nearest_neighbours(x,y,trace_x,trace_y,nx,ny)
        # Find the nearest neighbours for a set of points
        maptox = zeros(Int64,nx,ny)
        maptoy = zeros(Int64,nx,ny)
        
        for i = 1:nx
            for j = 1:ny
                ii = argmin(abs.(x.-trace_x[i,j]))
                jj = argmin(abs.(y.-trace_y[i,j]))
                maptox[i,j] = ii
                maptoy[i,j] = jj
            end
        end

        return maptox, maptoy
    end




end