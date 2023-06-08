```@setup tutorial
using IFSintegrals
using Plots
```

# Constructing Fractals and Fractal measures
An iterated function system is a set of $M$ affine contraction maps, or _similarities_.

### Similarities
Each similarity is of the form:

$s_m(x)=r_mA_mx + \delta_m,\quad x\in\mathbb{R}^n,$

where $r_m\in(0,1)$ is the contraction factor, $A_m\in\R^{n\times n}$ is a rotation matrix, and $\delta\in\R^{n}$ is a translation vector.

```@docs
Similarity
```
### Attractors of IFSs, as fractal measures
The iterated function system may be interpreted as a set of these maps, $\{s_m\}_{m=1}^M$. The _attractor_ $\Gamma$ is the unique non-empty compact set which satisfies

$$\Gamma = S(\Gamma):=\bigcup_{m=1}^M s_m(\Gamma).$$

We can construct the map $S$ defined above as follows.
```@REPL tutorial
using IFSintegrals
s₁ = Similarity(1/3,2/3)
s₂ = Similarity(1/3,2/3)
S = [s₁,s₂]
x = rand()
S(x) # applies the map S to the point x
```
This is the IFS for the Cantor Set, and this can be converted into an ```InvariantMeasure``` with the command
```@REPL tutorial
 Γ = InvariantMeasure(S)
```

 The outer constructor for ```InvariantMeasure``` constructs other properties, such as diameter and Hausdorff dimension, which describe this fractal measure.
 
 Now let's try a more complicated example. The second (translation) argument of `Similarity` can also be a vector. For example:
```@example tutorial
ρ = 0.41
IFS = [
    Similarity(ρ,[0,0])
    Similarity(ρ,[1-ρ,0])
    Similarity(ρ,[(1-ρ)/2,sqrt(3)*(1-ρ)/2])
    Similarity(ρ,[(1-ρ)/2,(1-ρ)/(2*sqrt(3))])
    ]
 Γ = InvariantMeasure(IFS)
 nothing #hide
```
We can plot an approximation attractors in $\mathbb{R}$ or $\mathbb{R}^2$, as follows:
```@example tutorial
using Plots
plot(Γ,color = "black", markersize=0.75,label="My fractal")
```

## Preset fractals
Several preset fractals are available. These may be called with:
```@REPL tutorial
CantorSet()
CantorDust()
Sierpinski() # triangle/gasket
KochCurve() # Koch curve
KochFlake() # Koch snowflake
Vicsek() # Vicsek fractal
SquareFlake() # Minkowski snowflake
Carpet() # Sierpinski carpet
Dragon() # Heighway Dragon
```
Some of these presets have optional arguments, for example `CantorSet` and `CantorDust` have the default `contraction=1/3`, but this can be set to any value in $(0,1/2]$.

## Manipulating fractals

Fractals can be translated, stretched and rotated very easily, using simple arithmetic syntax. For example,
```@example tutorial
Γ = Sierpinski()
plot(Γ, markersize=0.75,
label="Sierpinski Triangle (default)")

Γ_shift = 1.5*Γ + [-2,0.5]
plot!(Γ_shift, markersize=0.75, 
label="Sierpinski Triangle (translated and stretched)",aspect_ratio=1.0)

```

The following type also describes a measure whose support is an IFS attractor. It is used when meshing a fractal, with self-similar copies of its self.
```@docs
SelfSimilarSubMeasure
```