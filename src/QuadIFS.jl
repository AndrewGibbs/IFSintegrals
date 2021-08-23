module QuadIFS

using Base: Float64
export Similarity, Attractor, SubAttractor, Fractal, sketch_attractor, CantorSet,
        Sierpinski, sim_map,
        barycentre_rule, subdivide_indices, eval_green_double_integral,
        CantorDust, eval_green_single_integral_fixed_point, fixed_point,
        HelhmoltzGreen2D, dist, singular_elliptic_double_integral,
        BIO, DiscreteBIO, SingleLayer, Projection, single_layer_potential
        get_diameter, get_diam_long, h_dist
include("fractals.jl")
include("quadrature.jl")
include("green_kernels.jl")
include("BIOs.jl")
include("diam_approx.jl")

# routine below was copied from:
# https://discourse.julialang.org/t/converting-a-matrix-into-an-array-of-arrays/17038

function slicematrix(A::AbstractMatrix{T}) where T
        m, n = size(A)
        B = Vector{T}[Vector{T}(undef, n) for _ in 1:m]
        for i in 1:m
            B[i] .= A[i, :]
        end
        return B
    end

end
