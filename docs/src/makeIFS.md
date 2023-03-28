An iterated function system is a set of $M$ affine contraction maps, or _similarities_.

Each similarity is of the form:

$s_m(x)=r_mA_mx + \delta_m,\quad x\in\mathbb{R}^n,$

where $r_m\in(0,1)$ is the contraction factor, $A_m\in\R^{n\times n}$ is a rotation matrix, and $\delta\in\R^{n}$ is a translation vector.

```@docs
Similarity
```

The iterated function system may be interpreted as a set of these maps, $\{s_m\}_{m=1}^M$. The _attractor_ $\Gamma$ is the unique non-empty compact set which satisfies

$\Gamma = S(\Gamma):=\bigcup_{m=1}^M s_m(\Gamma).$

We can construct the map $S$ defined above as follows.
```
s₂ = Similarity(1/3,2/3)
S = [s₁,s₂]
x = rand()
S(x) # applies the map S to the point x
```
This is the IFS for the Cantor Set.


The attractor $\Gamma$ will be the support of our measure $\mu$. There are other properties, such as diameter and Hausdorff dimension, which describe this type.

```@docs
SelfSimilarMeasure
```

```@docs
SelfSimilarSubMeasure
```