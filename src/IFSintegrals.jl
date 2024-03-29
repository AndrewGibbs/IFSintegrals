module IFSintegrals
#@doc read(joinpath(dirname(@__DIR__), "README.md"), String) IFSintegrals

using Base: Float64
using ProgressMeter
#import LinearOperators: AbstractLinearOperator, LinearOperator
import LinearMaps: LinearMap
import IterativeSolvers: gmres
import Plots: scatter!, scatter, plot, plot!, cgrad, Shape, current
import Base: -, \, *, +, getindex, isapprox
import Roots: find_zero, Bisection
import LinearAlgebra: norm, I, UniformScaling, eigvals, eigvecs, dot, det, cross
import StaticArrays: SVector, SMatrix
import LazySets: VPolygon, Singleton, LineSegment
import PolygonInbounds: inpoly2
import SpecialFunctions: besselh
export FractalMeasure, InvariantMeasure, SubInvariantMeasure, SelfSimilarFractal, UnionInvariantMeasure,
        CantorSet, Sierpinski, SquareFlake, KochFlake, CantorN, Carpet, Vicsek, Dragon, KochCurve,
        Similarity, sim_map, full_map,
        barycentre_rule, subdivide_indices, eval_green_double_integral,
        CantorDust, eval_green_single_integral_fixed_point, fixed_point,
        HelhmoltzGreen2D, dist, singular_elliptic_double_integral,
        SIO, DiscreteSIO, SingleLayerOperatorLaplace, SingleLayerOperatorHelmholtz,
        Projection, get_layer_potential, SingleLayerPotentialHelmholtz,
        single_layer_potential, far_field_pattern,
        chaos_quad, barycentre_uniform, long_bary, gauss_quad,
        slice, box, draw, draw!,
        get_H_minus_half_norm_function, get_H_minus_half_norm_function_from_matrix,
        s_energy, DihedralGroup,
        VolumePotential, Id
include("similarities.jl")
include("symmetry_groups.jl")
include("fractals.jl")
include("fractal_manipulations.jl")
include("partitioning.jl")
include("quadrature.jl")
include("singular_homogenous_integrals.jl")
include("nondisjoint_singularities.jl")
include("SIOs.jl")
include("Galerkin_quadrature.jl")
include("Galerkin_discretisation.jl")
include("projections.jl")
include("diam_approx.jl")
include("fractal_presets.jl")
include("Greens_kernels.jl")
include("operator_presets.jl")
include("scattering_functionals.jl")
include("plotting.jl")
include("Jacobi_matrices.jl")
include("VIE_operators.jl")

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
