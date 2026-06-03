# Hyperbolic Polyplet Enumerator

Goal: Extend OEIS A390200, the number of free `n`-celled polyplets in the `{4,5}` tessellation of the hyperbolic plane.  Polyplets are connected by edge or vertex adjacency; each square has 4 edge-neighbors and 8 vertex-neighbors.

Known prefix, as checked on 2026-06-03:

```text
n:    1  2   3   4    5     6      7
a(n): 1, 2, 10, 72, 710, 8026, 98353
```

## Implementation

`polyplets_exact.py` is the current exact Coxeter-group enumerator.

The exact backend models cells as cosets of the face stabilizer in the Coxeter group `[4,5]`, using matrices over `Q(sqrt(2), sqrt(5))`.  Coefficients are stored exactly as integer numerators over powers of two:

```text
(a + b*sqrt(2) + c*sqrt(5) + d*sqrt(10)) / 2^e
```

The local model checks:

- square stabilizer has order 8;
- vertex stabilizer has order 10;
- the polyplet neighbor set has size 12.

The enumerator canonicalizes each finite connected set by translating every cell in turn to the base cell and minimizing over the square stabilizer.

## Status

The program reproduces the full known prefix and gives the following candidate extension:

```text
n,count,known,status,candidates,seconds
1,1,1,ok,1,0.000
2,2,2,ok,12,0.009
3,10,10,ok,35,0.029
4,72,72,ok,238,0.253
5,710,710,ok,2150,1.951
6,8026,8026,ok,25540,24.498
7,98353,98353,ok,337234,347.390
8,1261889,,new,4725668,6753.712
```

Candidate extension:

```text
a(8) = 1261889
```

## References

- OEIS A390200.
- arXiv:2109.05331, *Extremal {p,q}-Animals*.
- arXiv:2206.14910, *Isoperimetric Formulas for Hyperbolic Animals*.
