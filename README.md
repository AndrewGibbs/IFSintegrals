[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/AndrewGibbs/IFSintegrals/HEAD)

## IFSintegrals

A toolbox for solving problems defined on attractors of iterated function systems (IFSs).
See 'Quadrature example.ipynb', `BEM Cantor set example.ipynb` and `BEM 3D example.ipynb` for examples.
[Click here to load these examples in your browswer](https://mybinder.org/v2/gh/AndrewGibbs/IFSintegrals/9a9377940a437d0d952bfcf9632c053d3e7573c3) (without the need to install any software).

I would strongly reccomend playing with the notebook files to understand what the code can do. But here is a brief explanation of the main routines.
The type `Similarity(r,d,A)` can be used to define a contracting similarity, with contraction r, translation d and rotation A.
An IFS `S` can be represented as an array of these types, which can be used to construct an `Attractor(S)` of the IFS.
The following attractors are included as presets: `CantorSet()`, `CantorDust()`, `Sierpinski()`. The attractor can be approximated using `sketch_attractor`.

Weights and nodes for the evaluation of integrals with respect to Hausdorff (or equivalent) measure can be obtained using `barycentre_rule`.
Certain classes of singular integrals can be evaluated using `eval_green_double_integral` and `eval_green_single_integral_fixed_point`.

Boundary Integral operators can be defined and discretised on attractors, using the types `BIO` and `DiscreteBIO`.
There are examples of these problems being solved in the notebook files, where the boundary integral equation for the Helmholtz equation is solved.

This is joint work with Andrea Moiola,  David Hewett, Simon Chandler-Wilde and Antonio Caetano.
