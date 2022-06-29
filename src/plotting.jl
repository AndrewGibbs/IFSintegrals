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

function draw(Γ::SelfSimilarFractal; markersize=0.1, color="black")
    X = sketch_attractor(Γ)
    scatter([X[j][1] for j=1:length(X)],[X[j][2] for j=1:length(X)],legend=:false,markerstrokewidth=0, markersize=markersize, markercolor=color)
end

function draw!(Γ::SelfSimilarFractal; markersize=0.1, color="black")
    X = sketch_attractor(Γ)
    scatter!([X[j][1] for j=1:length(X)],[X[j][2] for j=1:length(X)],legend=:false,markerstrokewidth=0, markersize=markersize, markercolor=color)
end