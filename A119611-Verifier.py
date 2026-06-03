"""
Exact Coxeter enumerator for free {4,5} hyperbolic polyominoes (Edge-only).
Validates the geometric backend against OEIS A119611.
"""

from __future__ import annotations
from dataclasses import dataclass
from functools import lru_cache
from time import perf_counter
import argparse

KNOWN_A119611 = {
    1: 1,
    2: 1,
    3: 2,
    4: 5,
    5: 16,
    6: 55,
    7: 224,
    8: 978,
    9: 4507,
    10: 21430,
    11: 104423,
    12: 517897,
    13: 2606185
    14: 13272978
    15: 68300918
    16: 354649465
    17: 1856271806
    18: 9785436128
    19: 51916467206
    20: 277045470555
}

Alg = tuple[int, int, int, int, int]
Mat = tuple[Alg, ...]
Key = Mat

ZERO: Alg = (0, 0, 0, 0, 0)
ONE: Alg = (0, 1, 0, 0, 0)
NEG_ONE: Alg = (0, -1, 0, 0, 0)
SQRT2: Alg = (0, 0, 1, 0, 0)
PHI: Alg = (1, 1, 0, 1, 0)

def alg_norm(e: int, a: int, b: int, c: int, d: int) -> Alg:
    if a == 0 and b == 0 and c == 0 and d == 0:
        return ZERO
    while e > 0 and a % 2 == 0 and b % 2 == 0 and c % 2 == 0 and d % 2 == 0:
        e -= 1
        a //= 2
        b //= 2
        c //= 2
        d //= 2
    return (e, a, b, c, d)

@lru_cache(maxsize=None)
def alg_add(x: Alg, y: Alg) -> Alg:
    ex, a, b, c, d = x
    ey, f, g, h, i = y
    e = max(ex, ey)
    sx = 1 << (e - ex)
    sy = 1 << (e - ey)
    return alg_norm(
        e,
        a * sx + f * sy,
        b * sx + g * sy,
        c * sx + h * sy,
        d * sx + i * sy,
    )

@lru_cache(maxsize=None)
def alg_mul(x: Alg, y: Alg) -> Alg:
    ex, a, b, c, d = x
    ey, f, g, h, i = y
    return alg_norm(
        ex + ey,
        a * f + 2 * b * g + 5 * c * h + 10 * d * i,
        a * g + b * f + 5 * (c * i + d * h),
        a * h + c * f + 2 * (b * i + d * g),
        a * i + b * h + c * g + d * f,
    )

@lru_cache(maxsize=None)
def alg_neg(x: Alg) -> Alg:
    return (x[0], -x[1], -x[2], -x[3], -x[4])

@lru_cache(maxsize=None)
def alg_sub(x: Alg, y: Alg) -> Alg:
    return alg_add(x, alg_neg(y))

def mat_identity() -> Mat:
    return (
        ONE, ZERO, ZERO,
        ZERO, ONE, ZERO,
        ZERO, ZERO, ONE,
    )

@lru_cache(maxsize=None)
def mat_mul(a: Mat, b: Mat) -> Mat:
    out: list[Alg] = []
    for i in range(3):
        for j in range(3):
            value = ZERO
            for k in range(3):
                value = alg_add(value, alg_mul(a[3 * i + k], b[3 * k + j]))
            out.append(value)
    return tuple(out)

def mat_key(matrix: Mat) -> Key:
    return matrix

@dataclass(frozen=True)
class Element:
    matrix: Mat
    inverse: Mat

    def __matmul__(self, other: Element) -> Element:
        return Element(
            mat_mul(self.matrix, other.matrix),
            mat_mul(other.inverse, self.inverse),
        )

def element_identity() -> Element:
    ident = mat_identity()
    return Element(ident, ident)

def coxeter_generators() -> list[Element]:
    s0 = (NEG_ONE, SQRT2, ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ONE)
    s1 = (ONE, ZERO, ZERO, SQRT2, NEG_ONE, PHI, ZERO, ZERO, ONE)
    s2 = (ONE, ZERO, ZERO, ZERO, ONE, ZERO, ZERO, PHI, NEG_ONE)
    return [Element(s, s) for s in (s0, s1, s2)]

def unique_closure(generators: list[Element]) -> list[Element]:
    identity = element_identity()
    seen: dict[Key, Element] = {mat_key(identity.matrix): identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for gen in generators:
            nxt = current @ gen
            key = mat_key(nxt.matrix)
            if key not in seen:
                seen[key] = nxt
                frontier.append(nxt)
    return list(seen.values())

class ExactPolyominoBurnsideEnumerator:
    def __init__(self) -> None:
        self.generators = coxeter_generators()
        self.face_group = unique_closure([self.generators[0], self.generators[1]])
        self.key_to_id: dict[Key, int] = {}
        self.id_to_key: list[Key] = []
        self.coset_reps: list[Element] = []
        
        self.base = self.coset_id(element_identity())
        self.neighbor_moves = self._neighbor_moves()
        self.neighbor_cache: dict[int, tuple[int, ...]] = {}
        self.inverse_cache: dict[int, tuple[Element, ...]] = {}

    def coset_id(self, element: Element) -> int:
        images = [mat_key(mat_mul(element.matrix, h.matrix)) for h in self.face_group]
        key = min(images)
        cell_id = self.key_to_id.get(key)
        if cell_id is None:
            cell_id = len(self.id_to_key)
            self.key_to_id[key] = cell_id
            self.id_to_key.append(key)
            self.coset_reps.append(element)
        return cell_id

    def _neighbor_moves(self) -> list[Element]:
        """
        MODIFIED FOR POLYOMINOES:
        Only generate the 4 strictly edge-adjacent neighbors using s2 (the edge reflection).
        """
        moves: dict[int, Element] = {}
        s2 = self.generators[2]
        
        for h in self.face_group:
            move = h @ s2
            cell_id = self.coset_id(move)
            if cell_id != self.base:
                moves.setdefault(cell_id, move)
                
        if len(moves) != 4:
            raise ValueError(f"CRITICAL ERROR: A regular {4,5} square must have exactly 4 edges. Found {len(moves)}.")
            
        return list(moves.values())

    def neighbors(self, cell_id: int) -> list[int]:
        cached = self.neighbor_cache.get(cell_id)
        if cached is not None:
            return list(cached)
        rep = self.coset_reps[cell_id]
        result = tuple(self.coset_id(rep @ move) for move in self.neighbor_moves)
        self.neighbor_cache[cell_id] = result
        return list(result)

    def anchor_transforms(self, target: int) -> tuple[Element, ...]:
        cached = self.inverse_cache.get(target)
        if cached is not None:
            return cached
        target_rep = self.coset_reps[target]
        transforms = tuple(
            Element(
                mat_mul(h.inverse, target_rep.inverse),
                mat_mul(target_rep.matrix, h.matrix),
            )
            for h in self.face_group
        )
        self.inverse_cache[target] = transforms
        return transforms

    def generate_rooted_animals(self, max_n: int) -> dict[int, list[frozenset[int]]]:
        animals = {i: [] for i in range(1, max_n + 1)}
        
        def explore(shape_list: list[int], untried_list: list[int], forbidden_set: set[int]):
            n = len(shape_list)
            animals[n].append(frozenset(shape_list))
            if n == max_n:
                return
            
            for i, v in enumerate(untried_list):
                new_shape = shape_list + [v]
                
                new_neighbors = []
                for neighbor in self.neighbors(v):
                    if neighbor > self.base and neighbor not in forbidden_set:
                        new_neighbors.append(neighbor)
                        forbidden_set.add(neighbor)
                        
                next_untried = untried_list[i + 1:] + new_neighbors
                explore(new_shape, next_untried, forbidden_set)
                
                for neighbor in new_neighbors:
                    forbidden_set.remove(neighbor)
                    
        print(f"Generating polyomino spanning tree up to n={max_n}...")
        initial_untried = [v for v in self.neighbors(self.base) if v > self.base]
        initial_forbidden = set([self.base]) | set(initial_untried)
        
        explore([self.base], initial_untried, initial_forbidden)
        return animals

    def shape_stabilizer_size(self, shape: frozenset[int]) -> int:
        sym_count = 0
        for target in shape:
            for transform in self.anchor_transforms(target):
                is_sym = True
                for cell in shape:
                    mapped_cell = self.coset_id(transform @ self.coset_reps[cell])
                    if mapped_cell not in shape:
                        is_sym = False
                        break
                if is_sym:
                    sym_count += 1
        return sym_count

    def calculate_free_counts(self, all_rooted: dict[int, list[frozenset[int]]]) -> dict[int, int]:
        free_counts = {}
        for n, rooted_shapes in all_rooted.items():
            if n == 0: continue
            
            total_sym = sum(self.shape_stabilizer_size(shape) for shape in rooted_shapes)
            divisor = 8 * n
            assert total_sym % divisor == 0, f"Math error: {total_sym} not divisible by {divisor}"
            free_counts[n] = total_sym // divisor
            
        return free_counts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=20)
    args = parser.parse_args()

    enumerator = ExactPolyominoBurnsideEnumerator()
    
    start_time = perf_counter()
    rooted_shapes = enumerator.generate_rooted_animals(args.max_n)
    generation_time = perf_counter() - start_time
    
    start_eval = perf_counter()
    free_counts = enumerator.calculate_free_counts(rooted_shapes)
    eval_time = perf_counter() - start_eval

    print("\n--- Edge-Only Polyomino Sequence ---")
    print("n,count,known,status")
    for n in range(1, args.max_n + 1):
        count = free_counts[n]
        known = KNOWN_A119611.get(n)
        status = "new" if known is None else ("ok" if known == count else "MISMATCH")
        print(f"{n},{count},{known if known is not None else ''},{status}")

    print(f"\nTime stats: Generation = {generation_time:.3f}s | Evaluation = {eval_time:.3f}s")

if __name__ == "__main__":
    main()
