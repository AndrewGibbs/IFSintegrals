# provides a sketch of an parent_measure in N topological dimensions
function sketch_attractor(Γ::SelfSimilarFractal; mem_const = 100000, start_count = 10)
    if isa(Γ,InvariantMeasure)
        N = Γ.spatial_dimension
    elseif isa(Γ,SubInvariantMeasure)
        N = Γ.parent_measure.spatial_dimension
    end
    X = [Γ.barycentre]
    num_its = floor(log(mem_const/(N*start_count))/log(length(Γ.IFS)))
    for _ = 1:num_its
        X = full_map(Γ.IFS, X)
    end
    return X
end

function matrix_of_vectors(Nx::Int64, Ny::Int64)
    array = Array{Vector{Float64}}(undef, Nx,Ny)
    for i in eachindex(array)
        array[i] = Vector{Float64}()
    end
    return array
end

function box(c1::Vector{Float64},c2::Vector{Float64},Nx::Int64,Ny::Int64)
    xmin = min(c1[1],c2[1])
    xmax = max(c1[1],c2[1])
    ymin = min(c1[2],c2[2])
    ymax = max(c1[2],c2[2])
    xh = (xmax - xmin)/(Nx-1)
    yh = (ymax - ymin)/(Ny-1)
    x = [j for j=-xmin:xh:xmax]
    y = [j for j=-ymin:yh:ymax]
    M = [[Vector{Float64}[] for _=1:Nx] for _=1:Ny ]
    for nx=1:Nx
        for ny=1:Ny
            M[nx,ny] = [x[nx],y[Ny-ny+1]]
        end
    end
    return M
end

function slice(c1::Vector{Float64},c2::Vector{Float64},z::Float64,Nx::Int64,Ny::Int64)
    Nx = max(Nx,2)
    Ny = max(Ny,2)
    xmin = min(c1[1],c2[1])
    xmax = max(c1[1],c2[1])
    ymin = min(c1[2],c2[2])
    ymax = max(c1[2],c2[2])
    xh = (xmax - xmin)/(Nx-1)
    yh = (ymax - ymin)/(Ny-1)
    x = [j for j=xmin:xh:xmax]
    y = [j for j=ymin:yh:ymax]
    Nx = length(x)
    Ny = length(y)
    M = matrix_of_vectors(Ny,Nx)
    for nx=1:Nx
        for ny=1:Ny
            M[ny,nx] = [x[nx],y[ny],z]
        end
    end
    return M,x,y
end

function fractal_pre_plot(Γ::SelfSimilarFractal, mem_const)
    if isa(Γ,InvariantMeasure)
        n = Γ.spatial_dimension
    else
        n = Γ.parent_measure.spatial_dimension
    end
    X = sketch_attractor(Γ,mem_const=mem_const)
    z = Float64[]
    if n == 2
        x = [X[j][1] for j=1:length(X)]
        y = [X[j][2] for j=1:length(X)]
        xyz = (x,y)
    elseif n==3 
        x = [X[j][1] for j=1:length(X)]
        y = [X[j][2] for j=1:length(X)]
        z = [X[j][3] for j=1:length(X)]
        xyz = (x,y,z)
    elseif n == 1
        x = X
        y = zeros(length(X))
        xyz = (x,y)
    else
        error("Can only plot in one and two spatial dimensions")
    end
    return xyz
end

function fractal_pre_plot(Γ::UnionInvariantMeasure, mem_const)
    X = Vector{Float64}[]
    Y = Vector{Float64}[]
    # p = plot!()
    for γ ∈ Γ.invariant_measures
        x,y = fractal_pre_plot(γ,mem_const)
        push!(X,x)
        push!(Y,y)
    end
    return X,Y
end

function adjust_defaults(kwargs)
    kwargs_dict = Dict(kwargs)
    !haskey(kwargs,:markersize) ? merge(kwargs_dict,Dict("markersize"=>0.1)) : nothing
    !haskey(kwargs,:color) ? merge(kwargs_dict,Dict("color"=>"black")) : nothing
    return kwargs_dict
end

"""
    plot(Γ::SelfSimilarFractal; markersize=0.1, color="black")
Provides a simple sketch of the parent_measure Γ, by repeatedly applying the IFS.

See also: [`plot!`](@ref)
"""
# function plot(Γ::FractalMeasure; kwargs...)#mem_const = 100000, kwargs...)
#     # plot()
#     p = plot!(Γ; kwargs...)
#     # # println(kwargs)
#     # x,y = fractal_pre_plot(Γ,mem_const)
#     # # scatter(x,y;kwargs...)
#     # # kwargs = adjust_defaults(kwargs)
#     # scatter(x,y;adjust_defaults(kwargs)...)
#     return p
# end

"""
    plot!(Γ::SelfSimilarFractal; markersize=0.1, color="black")
Similar to [`draw`](@ref), except it will draw on the current image.
"""
function plot(Γ::FractalMeasure; mem_const = 100000, mswidth=0, kwargs...)
    xyz = fractal_pre_plot(Γ,mem_const)
    p = scatter(xyz;
                markerstrokewidth=mswidth,
                kwargs...)
    return p
end

function plot!(Γ::FractalMeasure; mem_const = 100000, mswidth=0, kwargs...)
    xyz = fractal_pre_plot(Γ,mem_const)
    p = scatter!(xyz;
                markerstrokewidth=mswidth,
                kwargs...)
    return p
end

# function plot!(Γ::UnionInvariantMeasure; mem_const = 100000, kwargs...)
#     X = Vector{Float64}[]
#     Y = Vector{Float64}[]
#     # p = plot!()
#     for γ ∈ Γ.invariant_measures
#         x,y = fractal_pre_plot(γ,mem_const)
#         push!(X,x)
#         push!(Y,y)
#     end
#     p = scatter!(X,Y;kwargs...)
#     return p
# end

# the big boy plotting algorithm

function check_for_intersection(x₁::Vector{Float64},x₂::Vector{Float64},y₁::Vector{Float64},y₂::Vector{Float64})

    a = LineSegment(x₁,x₂)
    b = LineSegment(y₁,y₂)
    is_disjoint, intersection_point = isdisjoint(a, b, true)
    return !is_disjoint, intersection_point
end


function get_clockwise_angle_between_vectors(v₁::Vector{Float64},v₂::Vector{Float64})
    θ = acos(clamp(dot(v₁,v₂) /(norm(v₁)*norm(v₂)), -1, 1))
    if cross(vcat(v₁,0),vcat(v₂,0))[3] < 0 #
        ϑ = θ
    else
        ϑ = 2π - θ
    end
end

function unify_polygons(P₁::Vector{Vector{Float64}},P₂::Vector{Vector{Float64}})

    # cyclic index map
    cycdex(n,P) = mod(n-1,length(P))+1

    # create seperate graphs for each polygon, defined by their edges:
    P1_edges = Vector{Int64}[]
    for n=1:length(P₁)
        push!(P1_edges,[cycdex(n,P₁),cycdex(n+1,P₁)])
    end
    P2_edges = Vector{Int64}[]
    for n=1:length(P₂)
        push!(P2_edges,[cycdex(n,P₂),cycdex(n+1,P₂)])
    end

    # part one - find all intersections and store them
    P₁_intersection_edge_indices = [Int64[] for _=1:length(P₁)]
    P₂_intersection_edge_indices = [Int64[] for _=1:length(P₂)]
    P₁_intersection_point_indices = [Int64[] for _=1:length(P₁)]
    P₂_intersection_point_indices = [Int64[] for _=1:length(P₂)]
    intersection_points = Vector{Float64}[]
    intersection_point_index = 0#length(P₁)+length(P₂);

    for (index1,e₁) ∈ enumerate(P1_edges)
        # global intersection_point_index, P₁_intersection_edge_indices, P₂_intersection_edge_indices, P₁_intersection_point_indices, P₂_intersection_point_indices, intersection_points
        for (index2,e₂) ∈ enumerate(P2_edges)
            #yn_edges, P1_pts, P2_pts = are_lines_overlapping(P₁[e₁[1]],P₁[e₁[2]],P₂[e₂[1]],P₂[e₂[2]])
            yn_points, x = check_for_intersection(P₁[e₁[1]],P₁[e₁[2]],P₂[e₂[1]],P₂[e₂[2]])
            if yn_points # two edges cross each other in the usual way, at just one point
                # check for similar intersection points which already existing
                matching_intersection_index = 0
                reuse_intersection = false
                for (n,y) ∈ enumerate(intersection_points)
                    if isapprox(y,x,atol=1e-6)
                        matching_intersection_index = n
                        reuse_intersection = true
                        break
                    end
                end
                if reuse_intersection
                    push!(P₁_intersection_point_indices[index1],matching_intersection_index)
                    push!(P₂_intersection_point_indices[index2],matching_intersection_index)
                    union!(P₁_intersection_point_indices[index1])
                    union!(P₂_intersection_point_indices[index2])
                else
                    intersection_point_index +=1 
                    push!(intersection_points,Vector{Float64}(x))
                    push!(P₁_intersection_point_indices[index1],intersection_point_index)
                    push!(P₂_intersection_point_indices[index2],intersection_point_index)
                end
            end
            if yn_points ## if there's been any intersection for these two edges
                push!(P₁_intersection_edge_indices[index1],index2)
                push!(P₂_intersection_edge_indices[index2],index1)
            end
        end
    end

    if intersection_point_index == 0 # no edge intersections
        inpoly2(P₁,P₂)[:,1]  == zeros(Bool,length(P₁)) ? P1_not_in_P2 = true : P1_not_in_P2 = false
        inpoly2(P₂,P₁)[:,1]  == zeros(Bool,length(P₂)) ? P2_not_in_P1 = true : P2_not_in_P1 = false
        if P1_not_in_P2 && P2_not_in_P1
            polygons_intersect = false
            union_of_polygons_vertices = Vector{Float64}[]
        elseif P1_not_in_P2 && !P2_not_in_P1
            polygons_intersect = true
            union_of_polygons_vertices = P₁
        elseif P2_not_in_P1 && !P1_not_in_P2
            polygons_intersect = true
            union_of_polygons_vertices = P₂
        else
            error("exhausted all cases, cant find one")
        end
    else # there are edge intersections, so taking the union is difficult
        polygons_intersect = true

        # part two

        combined_edges = Vector{Int64}[]
        combined_verticles = vcat(P₁,P₂,intersection_points)
        intersection_point_shift = length(P₁)+length(P₂)

        for (index1,e₁) ∈ enumerate(P1_edges)
            # global combined_edges, combined_verticles
            if P₁_intersection_point_indices[index1] == []
                push!(combined_edges,e₁)
            else
                intersection_points_ordered = sortperm([norm(P₁[e₁[1]]-y) for y ∈ intersection_points[P₁_intersection_point_indices[index1]]])
                
                # check if intersection points are the same as anything we've seen in either P. If so, make a new vector of indices for them.


                local_intersection_indices = Int64[]
                for n ∈ intersection_points_ordered
                    vertex_match = false
                    for (m,v) ∈ enumerate(vcat(P₁,P₂))
                        if isapprox(intersection_points[P₁_intersection_point_indices[index1][n]],v)
                            vertex_match = true
                            push!(local_intersection_indices,m)
                            break
                        end
                    end
                    if !vertex_match
                        push!(local_intersection_indices,intersection_point_shift+P₁_intersection_point_indices[index1][n])
                    end
                end
                #connect first existing vertex to first intersection point
                push!(combined_edges,[e₁[1],local_intersection_indices[1]])
                
                # connect intermediate vertices
                if length(P₁_intersection_point_indices[index1])>1
                    for n=1:(length(P₁_intersection_point_indices[index1])-1)
                        push!(combined_edges,[local_intersection_indices[n],local_intersection_indices[n+1]])
                    end
                end
                
                # connect final intersection point to second vertex
                push!(combined_edges,[local_intersection_indices[end],e₁[2]])

            end
        end

        num_P1_edges = length(P₁)
        for (index2,e₂) ∈ enumerate(P2_edges)
            if P₂_intersection_point_indices[index2] == []
                push!(combined_edges,[e₂[1]+num_P1_edges,e₂[2]+num_P1_edges])
            else
                intersection_points_ordered = sortperm([norm(P₂[e₂[1]]-y) for y ∈ intersection_points[P₂_intersection_point_indices[index2]]])
                
                local_intersection_indices = Int64[]
                for n ∈ intersection_points_ordered
                    vertex_match = false
                    for (m,v) ∈ enumerate(vcat(P₁,P₂))
                        if isapprox(intersection_points[P₂_intersection_point_indices[index2][n]],v)
                            vertex_match = true
                            push!(local_intersection_indices,m)
                            break
                        end
                    end
                    if !vertex_match
                        push!(local_intersection_indices,intersection_point_shift+P₂_intersection_point_indices[index2][n])
                    end
                end

                #connect first existing vertex to first intersection point
                push!(combined_edges,[e₂[1]+num_P1_edges,local_intersection_indices[1]])
                
                # connect intermediate vertices
                if length(P₂_intersection_point_indices[index2])>1
                    for n=1:(length(P₂_intersection_point_indices[index2])-1)
                        push!(combined_edges,[local_intersection_indices[n],local_intersection_indices[n+1]])
                    end
                end
                
                # connect final intersection point to second vertex
                push!(combined_edges,[local_intersection_indices[end],e₂[2]+num_P1_edges])

            end
        end

        # part two and a half. delete edges of zero length
        delete_index_list = Int64[]

        for (n,v) ∈ enumerate(combined_verticles)
            for (m,w) ∈ enumerate(combined_verticles)
                if m>n && isapprox(v,w,atol=1e-6)
                    # replace all m's with n's
                    for e ∈ combined_edges
                        e[e.==m].=n
                    end
                end
            end
        end


        for (n,e) ∈ enumerate(combined_edges)
            if e[1] == e[2]
                push!(delete_index_list,n)
            end
        end
        # end

        deleteat!(combined_edges, delete_index_list)

        # part three. traverse the graph in a roughly clockwise sense.

        # choose the starting point
        P1_points_in_P2_shape = inpoly2(P₁,P₂)[:,1]
        if P1_points_in_P2_shape == zeros(Bool,length(P₁))
            P2_points_in_P1_shape = inpoly2(P₂,P₁)[:,1]
            if P2_points_in_P1_shape == zeros(Bool,length(P₁))
                error("Both polygons are inside of each other :/")
            else
                start_node_index = findfirst(P2_points_in_P1_shape .== 0) + length(P₁)
            end
        else
            start_node_index = findfirst(P1_points_in_P2_shape .== 0)
        end

        # set things up for while loop
        current_edge = [0,0]
        for e ∈ combined_edges
            # global current_edge
            if e[1] == start_node_index
                current_edge = e
                break # can only be one of these
            end
        end
        current_node_index = current_edge[2]

        union_vertex_indices = [start_node_index,current_node_index]

        counter = 0
        while current_node_index != start_node_index
            counter +=1
            if counter>length(combined_verticles)
                error("stuck in a never ending loop")
            end
            # global current_edge, current_node_index, union_vertex_indices
            adjacent_edges = Vector{Int64}[] # edges adjacent to current node
            adjacent_edge_indices = Int64[]
            for (n,e) ∈ enumerate(combined_edges)
                if (e != current_edge) # ignore self case
                    if e[1] == current_node_index
                        push!(adjacent_edges,e)
                        push!(adjacent_edge_indices,n)
                    elseif e[2] == current_node_index
                        push!(adjacent_edges,[e[2],e[1]])
                        push!(adjacent_edge_indices,n)
                    end
                end
            end

            v₁ = combined_verticles[current_edge[2]]-combined_verticles[current_edge[1]]
            angles = [get_clockwise_angle_between_vectors(v₁,combined_verticles[e[1]]-combined_verticles[e[2]]) for e ∈ adjacent_edges]
            
            chosen_adjacent_edge_index = argmin(angles)
            push!(union_vertex_indices, adjacent_edges[chosen_adjacent_edge_index][2])
            # move to next the best adjacent edge
            current_edge = combined_edges[adjacent_edge_indices[chosen_adjacent_edge_index]]
            current_node_index = current_edge[2]

        end
        # step 4. clean up. some vertices might be very close together

        union_of_polygons_vertices = combined_verticles[union_vertex_indices]
        similar_points = true
        while similar_points
            similar_points = false
            for n=1:(length(union_of_polygons_vertices)-1)
                if isapprox(union_of_polygons_vertices[n],union_of_polygons_vertices[n+1],atol=1e-4)
                    similar_points = true
                    deleteat!(union_of_polygons_vertices,n)
                    break
                end
            end
        end
        pop!(union_of_polygons_vertices)

        # union_of_polygons_vertices = combined_verticles[union_vertex_indices[1:(end-1)]]
    end
    return polygons_intersect, union_of_polygons_vertices
end

function sketch_attractor_boundary(Γ::SelfSimilarFractal, levels::Int64; mem_const = 100000)

    Γ.spatial_dimension != 2 ? error("plotting only defined for fractals with two spatial dimensions") : nothing

    unit_square = [[-1/2,-1/2],[-1/2,1/2],[1/2,1/2],[1/2,-1/2]]
    # stretched_square = (1.0+0.1*rand())*Γ.diameter.*unit_square
    rand_vals = rand(length(Γ.IFS))
    # K = [x + Γ.barycentre for x∈ stretched_square]
    K = [[x + Γ.barycentre for x∈ ((1.0+0.1*rand_vals[m])*Γ.diameter.*unit_square)] for m=1:length(Γ.IFS)]

    for ℓ_=1:levels
        # S□ = [Vector{Vector{Float64}}([sim_map(s,x) for x∈K ]) for s∈Γ.IFS]
        if ℓ_==1
            S□ = [Vector{Vector{Float64}}([sim_map(Γ.IFS[m],x) for x∈K[m] ]) for m=1:length(Γ.IFS)]
        else
            S□ = [Vector{Vector{Float64}}([sim_map(s,x) for x∈K ]) for s∈Γ.IFS]
        end
        # println(ℓ_)
        count = 1
        while length(S□) > 1
            count +=1
            # println("\t",count)
            for n=2:length(S□)
                is_intersection, US□ = unify_polygons(S□[1],S□[n])
                if is_intersection
                    S□[1] = US□
                    deleteat!(S□,n)
                    break
                end
            end
        end
        K = S□[1]
        length(K)>mem_const ? break : nothing
    end

    return push!(K,K[1])
end

function polygonise_mesh(mesh::Vector{SubInvariantMeasure{V,M}}, prefractal_guess::Vector{Vector{Float64}}) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    
    Γ = mesh[1].parent_measure

    # prefractal_guess = sketch_attractor_boundary(Γ, levels, mem_const=mem_const)

    N = length(mesh)

    mesh_shapes = [Shape([(0.0,0.0)]) for _=1:N]

    if N == 1
        mesh_shapes[1] = Shape([y[1] for y ∈ prefractal_guess], [y[2] for y ∈ prefractal_guess])
    else
        for (n,mesh_el) ∈ enumerate(mesh)
            Y = [x for x∈prefractal_guess]
            for mᵢ ∈ reverse(mesh_el.index)
                Y = [sim_map(Γ.IFS[mᵢ],y) for y∈Y]
            end
            mesh_shapes[n] = Shape([y[1] for y ∈ Y], [y[2] for y ∈ Y])
        end
    end
    return mesh_shapes
end

"""
    plot(mesh::Vector{SubInvariantMeasure{V,M}}, vals::Vector{Float64};
            colour_map = :jet, linewidth =0.0,
            levels::Int64 = 3, mem_const = 100000, 
            prefractal_guess::Vector{Vector{Float64}} = sketch_attractor_boundary(mesh[1].parent_measure::SelfSimilarFractal, levels, mem_const=mem_const),
            kwargs...) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}

Provides a colourmap-style plot of a fractal ```mesh````,
where the colours are determined by ```vals```
"""
function plot(mesh::Vector{SubInvariantMeasure{V,M}}, vals::Vector{Float64};
        colour_map = :jet, linewidth =0.0,
        levels::Int64 = 3, mem_const = 100000, 
        prefractal_guess::Vector{Vector{Float64}} = sketch_attractor_boundary(mesh[1].parent_measure::SelfSimilarFractal, levels, mem_const=mem_const),
        kwargs...) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}

    length(mesh) != length(vals) ? error("values vector and mesh need to be same size") : nothing

    plot(polygonise_mesh(mesh, prefractal_guess);
        c=colour_map, mc=colour_map, fill_z=permutedims(vals), labels=:none, linewidth=linewidth, kwargs...)

end

function plot!(mesh::Vector{SubInvariantMeasure{V,M}}, vals::Vector{Float64};
    colour_map = :jet, linewidth =0.0,
    levels::Int64 = 3, mem_const = 100000, 
    prefractal_guess::Vector{Vector{Float64}} = sketch_attractor_boundary(mesh[1].parent_measure::SelfSimilarFractal, levels, mem_const=mem_const),
    kwargs...) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}

length(mesh) != length(vals) ? error("values vector and mesh need to be same size") : nothing

plot!(polygonise_mesh(mesh, prefractal_guess);
    c=colour_map, mc=colour_map, fill_z=permutedims(vals), labels=:none, linewidth=linewidth, kwargs...)

end