[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/AndrewGibbs/IFSintegrals/HEAD)

# IFSintegrals

A toolbox for solving problems defined on attractors of iterated function systems (IFSs). This is joint work with Andrea Moiola, David Hewett, Simon Chandler-Wilde and António Caetano.

See `Quadrature example.ipynb`, `BEM Cantor set example.ipynb` and `BEM 3D example.ipynb` for examples.
[Click here to load these interactive examples in your browswer](https://mybinder.org/v2/gh/AndrewGibbs/IFSintegrals/HEAD) (without the need to install any software). I would strongly reccomend playing with these notebook files to understand what the code can do. `Quadrature example.ipynb` contains a sufficient introduction to the toolbox, and the introductory section of this should be understood first, even if you are not interested in understanding how the quadrature works.

## Installation:
To install, type the following into Julia:

`using Pkg`

`Pkg.add("https://github.com/AndrewGibbs/IFSintegrals.git")`

Some users have had issues with the GR backend during installation. If this occurs, install manually with

`Pkg.add("gr")`

and try again.

## Quadrature:
Weights and nodes for the evaluation of integrals with respect to Hausdorff (or equivalent) measure can be obtained using `barycentre_rule`, which is a generalisation of the midpoint rule to IFS attractors.

For IFS attractors which are subsets of $\mathbb{R}$, Gaussian quadrature is available using `gauss_quad`.

Certain classes of singular integrals can be evaluated using `eval_green_double_integral` and `eval_green_single_integral_fixed_point`.


## BEM:
Boundary Integral operators can be defined and discretised on attractors, using the types `BIO` and `DiscreteBIO`.
There are examples of these problems being solved in the notebook files, where the boundary integral equation for the Helmholtz equation is solved.

## Bibliography
* **Numerical Quadrature for Singular Integrals on Fractals**, A. Gibbs, D. P. Hewett, A. Moiola, [published article](https://link.springer.com/article/10.1007/s11075-022-01378-9), [arxiv preprint](http://arxiv.org/abs/2112.11793).
* **A Hausdorff-measure boundary element method for acoustic scattering by fractal screens**, António M. Caetano, Simon N. Chandler-Wilde, Andrew Gibbs, David P. Hewett, Andrea Moiola, [arxiv preprint](https://arxiv.org/abs/2212.06594)
* **A Stable Stieltjes Technique for Computing Orthogonal Polynomials and Jacobi Matrices Associated with a
Class of Singular Measures**, G. Mantica, [published article](https://link.springer.com/article/10.1007/BF02437506).

