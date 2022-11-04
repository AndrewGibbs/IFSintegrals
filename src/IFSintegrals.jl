module IFSintegrals

using Base: Float64
using ProgressMeter
import Plots: scatter!, scatter
import Base: -, \, *, getindex, isapprox
import Roots: find_zero, Bisection
import LinearAlgebra: norm, I, UniformScaling, eigvals, eigvecs
import StaticArrays: SVector, SMatrix
using Polyhedra
export Similarity, Attractor, SubAttractor, SelfSimilarFractal, sketch_attractor, CantorSet,
        Sierpinski, SquareFlake, KochFlake, CantorN, Carpet,
        sim_map, full_map,
        barycentre_rule, subdivide_indices, eval_green_double_integral,
        CantorDust, eval_green_single_integral_fixed_point, fixed_point,
        HelhmoltzGreen2D, dist, singular_elliptic_double_integral,
        SIO, DiscreteSIO, SingleLayer, Projection,
        single_layer_potential, far_field_pattern,
        chaos_quad, barycentre_uniform, long_bary, gauss_quad,
        slice, box, draw, draw!,
        get_H_minus_half_norm_function, get_H_minus_half_norm_function_from_matrix,
        s_energy, DihedralGroup
include("similarities.jl")
include("fractals.jl")
include("partitioning.jl")
include("quadrature.jl")
include("green_kernels.jl")
include("singular_homogenous_integrals.jl")
include("SIOs.jl")
include("projections.jl")
include("diam_approx.jl")
include("presets.jl")
include("screen_scattering.jl")
include("plotting.jl")
include("Jacobi_matrices.jl")
include("symmetry_groups.jl")
include("nondisjoint_singularities.jl")

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
