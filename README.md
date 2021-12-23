[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/AndrewGibbs/IFSintegrals/HEAD)

# IFSintegrals

A toolbox for solving problems defined on attractors of iterated function systems (IFSs).

See 'Quadrature example.ipynb', `BEM Cantor set example.ipynb` and `BEM 3D example.ipynb` for examples.
[Click here to load these examples in your browswer](https://mybinder.org/v2/gh/AndrewGibbs/IFSintegrals/HEAD) (without the need to install any software).

I would strongly reccomend playing with the notebook files to understand what the code can do. But here is a brief explanation of the main routines. 'Quadrature example.ipynb' contains a sufficient introduction to the toolbox, even if you are not interested in understanding how the quadrature works.

## Quadrature:
Weights and nodes for the evaluation of integrals with respect to Hausdorff (or equivalent) measure can be obtained using `barycentre_rule`. Certain classes of singular integrals can be evaluated using `eval_green_double_integral` and `eval_green_single_integral_fixed_point`.

## BEM:
Boundary Integral operators can be defined and discretised on attractors, using the types `BIO` and `DiscreteBIO`.
There are examples of these problems being solved in the notebook files, where the boundary integral equation for the Helmholtz equation is solved.

This is joint work with Andrea Moiola, David Hewett, Simon Chandler-Wilde and António Caetano.

## Bibliography
* [**Numerical Quadrature for Singular Integrals on Fractals**](http://arxiv.org/abs/2112.11793)<br>
