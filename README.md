[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/AndrewGibbs/IFSintegrals/HEAD)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://andrewgibbs.github.io/IFSintegrals/dev/)
[![CI](https://github.com/AndrewGibbs/IFSintegrals/actions/workflows/CI.yml/badge.svg)](https://github.com/AndrewGibbs/IFSintegrals/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/AndrewGibbs/IFSintegrals/branch/master/graph/badge.svg?token=H4ZHZU4TMH)](https://codecov.io/gh/AndrewGibbs/IFSintegrals)
[![DOI](https://zenodo.org/badge/397262754.svg)](https://zenodo.org/doi/10.5281/zenodo.13642679)

# IFSintegrals

A toolbox for solving problems defined on attractors of iterated function systems (IFSs). This project is an implementation of multiple research projects, with Andrea Moiola, David Hewett, Simon Chandler-Wilde, António Caetano, Joshua Bannister, Botond Major and Jeevon Greewal. In particular, the following project students have contributed to the development of the code: Joshua Bannister, Jeevon Greewal.

See `Quadrature example.ipynb`, `BEM Cantor set example.ipynb` and `BEM 3D example.ipynb` for examples.
[Click here to load these interactive examples in your browser](https://mybinder.org/v2/gh/AndrewGibbs/IFSintegrals/HEAD) (without the need to install any software). I would strongly recommend playing with these notebook files to understand what the code can do. `Quadrature example.ipynb` contains a sufficient introduction to the toolbox, and the introductory section of this should be understood first, even if you are not interested in understanding how the quadrature works.

## Installation:
To install, type the following into Julia:

`using Pkg`

`Pkg.add(url="https://github.com/AndrewGibbs/IFSintegrals.git")`

## Quadrature:
Weights and nodes for the evaluation of integrals with respect to Hausdorff (or equivalent) measure can be obtained using `barycentre_rule`, which is a generalisation of the midpoint rule to IFS attractors.

For IFS attractors which are subsets of $\mathbb{R}$, Gaussian quadrature is available using `gauss_quad`.

Certain classes of singular integrals can be evaluated using `eval_green_double_integral` and `eval_green_single_integral_fixed_point`.


## BEM:
Boundary Integral operators can be defined and discretised on attractors, using the types `BIO` and `DiscreteBIO`.
There are examples of these problems being solved in the notebook files, where the boundary integral equation for the Helmholtz equation is solved.

## Bibliography
* **Numerical Quadrature for Singular Integrals on Fractals**, A. Gibbs, D. P. Hewett, A. Moiola, [published article](https://link.springer.com/article/10.1007/s11075-022-01378-9), [arxiv preprint](http://arxiv.org/abs/2112.11793).
* **Numerical evaluation of singular integrals on non-disjoint self-similar fractal sets**, Andrew Gibbs, David P. Hewett, Botond Major, [published article](https://link.springer.com/content/pdf/10.1007/s11075-023-01705-8.pdf), [arxiv preprint](https://arxiv.org/abs/2303.13141)
* **Integral equation methods for acoustic scattering by fractals**, A. M. Caetano, S. N. Chandler-Wilde, X. Claeys, A. Gibbs, D. P. Hewett, A. Moiola, [arxiv preprint](https://arxiv.org/abs/2309.02184)
* **A Hausdorff-measure boundary element method for acoustic scattering by fractal screens**, António M. Caetano, Simon N. Chandler-Wilde, Andrew Gibbs, David P. Hewett, Andrea Moiola, [arxiv preprint](https://arxiv.org/abs/2212.06594)
* **A Hausdorff-measure boundary element method for acoustic scattering by fractal screens**, António M. Caetano, Simon N. Chandler-Wilde, Andrew Gibbs, David P. Hewett, Andrea Moiola, [arxiv preprint](https://arxiv.org/abs/2212.06594)
* **A Stable Stieltjes Technique for Computing Orthogonal Polynomials and Jacobi Matrices Associated with a
Class of Singular Measures**, G. Mantica, [published article](https://link.springer.com/article/10.1007/BF02437506).

