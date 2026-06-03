# Hyperbolic Polyplet Enumerator

**Goal:** Extend [OEIS A390200](https://oeis.org/A390200), the number of free $n$-celled polyplets in the $\\{4,5\\}$ tessellation of the hyperbolic plane. Polyplets are connected by edge or vertex adjacency; each square has 4 edge-neighbors and 8 vertex-neighbors.

**Known Prefix (A390200), as checked on 2026-06-03:**  

n:    1, 2, 3, 4, 5, 6, 7  
a(n): 1, 2, 10, 72, 710, 8026, 98353

## File Structure

* **`Polyplets_Exact.py`** - The primary direct canonical expansion enumerator. Generates polyplets dynamically and canonicalizes their coordinate states against the face stabilizer.
* **`Redelmeier-and-Burnside's-Lemma-Verifier.py`** - An independent enumeration engine utilizing a memory-safe graph spanning tree to generate fixed roots, mathematically computing the free count via the Orbit-Stabilizer theorem.
* **`A119611-Verifier.py`** - A variant of the pipeline configured for strict edge-only connectivity to verify the algebraic core against an established sequence.
* **`Coxeter-Verifier.py`** - A symbolic verification script using SymPy to test that the base matrix generators strictly satisfy the $[4,5]$ Coxeter group relations.
* **`Output.json`** - Structured output data containing candidate counts, computational statistics, execution runtimes, and the verified extension data.

## Implementation & Methodology

The exact backend models cells as cosets of the face stabilizer in the Coxeter group $[4,5]$, using matrices over the algebraic number field $\mathbb{Q}(\sqrt{2}, \sqrt{5})$. Coefficients are stored exactly as integer numerators over powers of two to completely eliminate floating-point artifacts. 

The local model rigorously checks:
* The square stabilizer has order 8 ($D_4$).
* The vertex stabilizer has order 10 ($D_5$).
* The polyplet neighbor set has size 12.

To guarantee correctness, the repository uses two mathematically distinct enumeration philosophies that corroborate each other:

1. **Direct Canonicalization (`Polyplets_Exact.py`):** Expands the boundary dynamically and canonicalizes each finite connected set by translating every cell to the base cell and minimizing over the square stabilizer. 
2. **Fixed Spanning Tree & Orbit-Stabilizer (`Burnside_Enumerator.py`):** Uses a Redelmeier-style lexicographical spanning tree to generate fixed, rooted animals, bypassing memory constraints. Burnside's Lemma (the Orbit-Stabilizer theorem) is then applied analytically to the symmetries of the exact Coxeter matrices to extract the free count.

## Verification

To prove the algebraic engine, the neighbor generation logic was restricted to edge-only adjacency ($s_2$ reflections) and successfully reproduced the known prefix of **OEIS A119611** (Strict $\\{4,5\\}$ Polyominoes) perfectly up to $n=12$. The continuous geometric space is further verified symbolically via `Coxeter-Verifier.py`.

## Status & Candidate Extension

Both the direct canonicalizer and the Burnside enumerator reproduce the full known prefix and mathematically converge on the following candidate extension for $n=8$:

| $n$ | Count | Status | Candidates | Seconds |
| :--- | :--- | :--- | :--- | :--- |
| 1 | 1 | ok | 1 | 0.000 |
| 2 | 2 | ok | 12 | 0.009 |
| 3 | 10 | ok | 35 | 0.029 |
| 4 | 72 | ok | 238 | 0.253 |
| 5 | 710 | ok | 2150 | 1.951 |
| 6 | 8026 | ok | 25540 | 24.498 |
| 7 | 98353 | ok | 337234 | 347.390 |
| **8** | **1261889** | **new** | **4725668** | **6753.712** |

**Candidate extension:**
$$a(8) = 1261889$$

## References
* [OEIS A390200](https://oeis.org/A390200)
* [OEIS A119611](https://oeis.org/A119611)
* [arXiv:2109.05331](https://arxiv.org/abs/2109.05331), Extremal $\{p,q\}$-Animals
* [arXiv:2206.14910](https://arxiv.org/abs/2206.14910), Isoperimetric Formulas for Hyperbolic Animals
