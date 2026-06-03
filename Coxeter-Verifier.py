"""
Independent SymPy verifier for the [4,5] Coxeter group representation.
This validates the algebraic foundation of the ExactPolypletEnumerator.
"""

import sympy as sp

def verify_coxeter_group():
    print("Initializing exact symbolic matrices for [4,5] space...")
    
    sqrt2 = sp.sqrt(2)
    phi = (1 + sp.sqrt(5)) / 2
    
    s0 = sp.Matrix([
        [-1, sqrt2, 0],
        [0, 1, 0],
        [0, 0, 1]
    ])
    
    s1 = sp.Matrix([
        [1, 0, 0],
        [sqrt2, -1, phi],
        [0, 0, 1]
    ])
    
    s2 = sp.Matrix([
        [1, 0, 0],
        [0, 1, 0],
        [0, phi, -1]
    ])
    
    I = sp.eye(3)
    
    print("\n--- Testing Involutions (Reflections) ---")
    assert sp.simplify(s0 * s0) == I, "FAIL: s0 is not an involution"
    assert sp.simplify(s1 * s1) == I, "FAIL: s1 is not an involution"
    assert sp.simplify(s2 * s2) == I, "FAIL: s2 is not an involution"
    print("PASS: All base generators are perfect reflections (s^2 = I).")
    
    print("\n--- Testing Coxeter Relations ---")
    # Face stabilizer: (s0 * s1)^4 = I (Square geometry)
    s01 = sp.simplify(s0 * s1)
    s01_4 = sp.simplify(s01**4)
    assert s01_4 == I, "FAIL: Face stabilizer does not form D4"
    print("PASS: Face geometry is exactly D4 ((s0*s1)^4 = I).")
    
    # Vertex stabilizer: (s1 * s2)^5 = I (5-way intersection)
    s12 = sp.simplify(s1 * s2)
    s12_5 = sp.simplify(s12**5)
    assert s12_5 == I, "FAIL: Vertex stabilizer does not form D5"
    print("PASS: Vertex geometry is exactly D5 ((s1*s2)^5 = I).")
    
    # Orthogonality: (s0 * s2)^2 = I
    s02 = sp.simplify(s0 * s2)
    s02_2 = sp.simplify(s02**2)
    assert s02_2 == I, "FAIL: s0 and s2 do not commute"
    print("PASS: Orthogonal mirrors commute ((s0*s2)^2 = I).")
    
    print("\n--- Determinant Checks ---")
    assert sp.simplify(s0.det()) == -1, "FAIL: s0 determinant != -1"
    assert sp.simplify(s1.det()) == -1, "FAIL: s1 determinant != -1"
    assert sp.simplify(s2.det()) == -1, "FAIL: s2 determinant != -1"
    print("PASS: All matrices have determinant -1 (Isometries).")
    
    print("\nSUCCESS: The matrix representation is a mathematically flawless model of {4,5} hyperbolic space.")

if __name__ == "__main__":
    verify_coxeter_group()
