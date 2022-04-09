module makeplane
    

    
    function construct_grid(ψ,θ,nx,ny,nz,Δψ,Δθ,H,Bfield;bc=[0.,1.])
        #=
            Domain is [0,1]×[0,2π), ψ and θ inputs zoom in on region
        =#

        ψ = range(ψ[1],stop=ψ[2],length=nx)[2:end-1]
        θ = range(θ[1],stop=θ[2],length=ny)
        Δψ = (ψ[2]-ψ[1])
        Δθ = (θ[2]-θ[1])
        
        if mnk[3] == 0
            error("Number of planes cannot be zero.")
        end
        
        Δζ = 2π #HARD CODED

        #=
            Load fns to avoid variable passing
        =#
        ini(pt) = internal_node_interp(pt,ψ,θ,Δψ,Δθ,m,n)
        
        B = zeros(m*n)
        for i = 1:m
            for j = 1:n
                ind = i + (j-1)*m
                B[ind] = norm(Bfield([ψ[i],θ[j],1.0]))
            end
        end
        # Preallocate planes

        ζ = 2π
        plane = trace(H,ψ,θ,ζ,m,n)
        f_plane = grid_interpolate(plane,ψ,bc,ζ,ini)
        
        ζ = -2π
        plane = trace(H,ψ,θ,ζ,m,n)
        b_plane = grid_interpolate(plane,ψ,bc,ζ,ini)

        gdata = grid_data(ψ,θ,Δψ,Δθ,Δζ,m,n,f_plane,b_plane,B)

        # Form an interpolation grid for recoving T data on
        # Pack grid data
        # f_plane = plane_parallel(Δζ,f_points,f_plane[2:-1:1,:],f_weights,f_modB)
        # b_plane = plane_parallel(-Δζ,b_points,b_plane[2:-1:1,:],b_weights,b_modB)
        
        # rgrid = RectangleGrid(Δψ:Δψ:ψ[1]-Δψ,θ[1]:Δθ:θ[2])
        
        # gdata = grid_data(ψ,θ,Δψ,Δθ,m,n,f_plane,b_plane,rgrid,modB)
        # Return points and interpolation data
        return gdata
    end


    function trace(H,ψ,θ,ζ,m,n)
        # Returns an array of [ψ,θ] points on the desired ζ plane
        plane = zeros(2,m*n)
        for i = 1:m
            for j = 1:n
                P = ODEProblem(H,[θ[j],ψ[i]],(0.0,ζ))
                sol = solve(P)
                plane[:,i+(j-1)*m] = mod.(sol.u[end],2π)[[2,1]] #flip to ψ,θ
            end
        end
        return plane
    end


    

    function grid_interpolate(plane,ψ,bc,ζ,ini)
        ind_inner = findall(x-> ψ[1] < x < ψ[end], plane[1,:])
        ind_exterior = collect(1:1:size(plane)[2])[filter(x->!(x in ind_inner),eachindex(ind_inner))]

        # boarder_inner = 
        # boarder_exterior = 

        i_points = zeros(Int64,5,length(ind_inner)) #Need 5 to store the row# and col#'s
        i_weights = zeros(Float64,4,length(ind_inner))

        e_points = zeros(Int64,length(ind_exterior)) #only need 1 for RHS storage
        e_weights = zeros(Float64,length(ind_exterior))
        for i = 1:length(ind_inner)
            i_points[1,i] = ind_inner[i]
            i_points[2:5,i],i_weights[:,i] = ini(plane[:,ind_inner[i]])
        end

        for i = 1:length(ind_exterior)
            # Need to convert from ind back to i,j
            e_points[i] = ind_exterior[i]
            # e_points[2,i] = mod1(ind_exterior[i],m)
            if ψ[end] > plane[1,ind_exterior[i]]
                e_weights[i] = bc[1]
            elseif ψ[1] < plane[1,ind_exterior[i]]
                e_weights[i] = bc[2]
            end
        end

        return plane_parallel(ζ,plane,i_points,i_weights,e_points,e_weights,ones(size(plane)[2]))
    end


        
    function internal_node_locate(pt,x,y,n)
        # pt = [ψ,θ]
        i₁ = findall(x-> x>=0.0, pt[1].-x)
        i₂ = i₁ + 1
        # Loopback on periodic
        j₁ = mod1(findall(x->x>=0.0, pt[2].-y)[end],n)
        j₂ = mod1(j₁ + 1,n)
        return i₁,i₂,j₁,j₂
    end

    function internal_node_interp(pt::Vector{Float64},x::StepRangeLen,y::StepRangeLen,Δx::Float64,Δy::Float64,m::Int64,n::Int64)
        # pt = [ψ,θ]
        li(i,j) = LinInd(i,j,m,n)

        Δ = Δx*Δy
        i₁ = findall(x-> x>=0.0, pt[1].-x)[end]
        i₂ = i₁ + 1
        # Loopback on periodic
        j₁ = mod1(findall(x->x>=0.0, pt[2].-y)[end],n)
        j₂ = mod1(j₁ + 1,n)
        
        # BUILD WEIGHTS
        x₁ = (pt[1] - x[i₁])

        if i₂ <= m
            x₂ = (x[i₂] - pt[1])
        else
            x₂ = (x[i₁]+Δx)-pt[1]
        end
        y₁ = (pt[2] - y[j₁])
        y₂ = (y[j₂] - pt[2])
        w₁₁ = x₂*y₂/Δ #Bottom left
        w₂₁ = x₁*y₂/Δ #Bottom right
        w₂₂ = x₁*y₁/Δ #Top right
        w₁₂ = x₂*y₁/Δ #Top left

        pts = [li(i₁,j₁),li(i₂,j₁),li(i₂,j₂),li(i₁,j₂)]
        weights = [w₁₁,w₂₁,w₂₂,w₁₂]

        return pts,weights
    end

    #=
        EXTERIOR NODES
    =#
    function external_node_interp(ij,bc)
        pts = [LinInd(ij[1],ij[2])]
        weights = [bc[2]]
        return pts, weights
    end




    #=
        LINEAR INDICES
    =#


    function LinInd(i,j,m,n)
        i < n ? j : j = 1
        return i + (j-1)*m
    end





end