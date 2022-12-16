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

"""
    draw(Γ::SelfSimilarFractal; markersize=0.1, color="black")
Provides a simple sketch of the parent_measure Γ, by repeatedly applying the IFS.

See also: [`draw!`](@ref)
"""
function draw(Γ::SelfSimilarFractal; markersize=1.0, color="black", grid=true, mem_const = 100000)
    if isa(Γ,InvariantMeasure)
        n = Γ.spatial_dimension
    else
        n = Γ.parent_measure.spatial_dimension
    end
    X = sketch_attractor(Γ,mem_const=mem_const)
    if n == 2
        scatter([X[j][1] for j=1:length(X)],[X[j][2] for j=1:length(X)],legend=:false,markerstrokewidth=0, markersize=markersize, markercolor=color, grid=grid)
    elseif n == 1
        scatter(X,zeros(length(X)),legend=:false,markerstrokewidth=0, markersize=markersize, markercolor=color, grid=grid)
    else
        error("Can only plot in one and two spatial dimensions")
    end
end

"""
    draw!(Γ::SelfSimilarFractal; markersize=0.1, color="black")
Similar to [`draw`](@ref), except it will draw on the current image.
"""
function draw!(Γ::SelfSimilarFractal; markersize=1.0, color="black", grid=true, mem_const = 100000)
    if isa(Γ,InvariantMeasure)
        n = Γ.spatial_dimension
    else
        n = Γ.parent_measure.spatial_dimension
    end
    X = sketch_attractor(Γ,mem_const=mem_const)
    if n == 2
        scatter!([X[j][1] for j=1:length(X)],[X[j][2] for j=1:length(X)],legend=:false,markerstrokewidth=0, markersize=markersize, markercolor=color, grid=grid)
    elseif n == 1
        scatter!(X,zeros(length(X)),legend=:false,markerstrokewidth=0, markersize=markersize, markercolor=color, grid=grid)
    else
        error("Can only plot in one and two spatial dimensions")
    end
end