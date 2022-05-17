module IFSintegrals

using Base: Float64
using ProgressMeter
import Plots: scatter!, scatter
import Base: -, \, getindex
import Roots: find_zero, Bisection
import LinearAlgebra: norm, I, UniformScaling
import StaticArrays: SVector, SMatrix
export Similarity, Attractor, SubAttractor, Fractal, sketch_attractor, CantorSet,
        Sierpinski, SquareFlake, KochFlake, CantorN, sim_map, full_map,
        barycentre_rule, subdivide_indices, eval_green_double_integral,
        CantorDust, eval_green_single_integral_fixed_point, fixed_point,
        HelhmoltzGreen2D, dist, singular_elliptic_double_integral,
        BIO, DiscreteBIO, SingleLayer,
        single_layer_potential, far_field_pattern,
        chaos_quad, barycentre_uniform,
        slice, box, draw, get_H_minus_half_norm_function, get_H_minus_half_norm_function_from_matrix
include("similarities.jl")
include("fractals.jl")
include("partitioning.jl")
include("quadrature.jl")
include("green_kernels.jl")
include("BIOs.jl")
include("projections.jl")
include("diam_approx.jl")
include("presets.jl")
include("scattering.jl")
include("plotting.jl")

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

function constant_vector(v::Vector)
    minimum(v)==maximum(v) ? x= true : x= false
    return x
end

end
