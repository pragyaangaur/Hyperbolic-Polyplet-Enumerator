"""Exact Coxeter enumerator for free {4,5} hyperbolic polyplets.

This program represents the geometric reflection matrices for [4,5] in the
number field Q(sqrt(2), sqrt(5)).  Matrix equality and coset keys are exact.

The implementation is still a direct expand/canonicalize enumerator, so the
first milestone is correctness through the known A390200 prefix n <= 7; the
second milestone is improving enumeration speed enough to push n = 8.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from time import perf_counter
import argparse
import json
import os
import pickle


KNOWN_A390200 = {
    1: 1,
    2: 2,
    3: 10,
    4: 72,
    5: 710,
    6: 8026,
    7: 98353,
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


def mat_inverse_unit(matrix: Mat) -> Mat:
    """Inverse of a 3x3 group matrix with determinant +/- 1."""
    a, b, c, d, e, f, g, h, i = matrix
    cof = (
        alg_sub(alg_mul(e, i), alg_mul(f, h)),
        alg_neg(alg_sub(alg_mul(b, i), alg_mul(c, h))),
        alg_sub(alg_mul(b, f), alg_mul(c, e)),
        alg_neg(alg_sub(alg_mul(d, i), alg_mul(f, g))),
        alg_sub(alg_mul(a, i), alg_mul(c, g)),
        alg_neg(alg_sub(alg_mul(a, f), alg_mul(c, d))),
        alg_sub(alg_mul(d, h), alg_mul(e, g)),
        alg_neg(alg_sub(alg_mul(a, h), alg_mul(b, g))),
        alg_sub(alg_mul(a, e), alg_mul(b, d)),
    )
    det = alg_add(alg_add(alg_mul(a, cof[0]), alg_mul(b, cof[3])), alg_mul(c, cof[6]))
    if det == ONE:
        sign = 1
    elif det == NEG_ONE:
        sign = -1
    else:
        raise ValueError(f"Expected determinant +/-1, got {det}")
    adj = (
        cof[0], cof[3], cof[6],
        cof[1], cof[4], cof[7],
        cof[2], cof[5], cof[8],
    )
    if sign == 1:
        return adj
    return tuple(alg_neg(x) for x in adj)


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
    # Reflection matrices from the geometric representation:
    # s0 row 0 = [-1, sqrt(2), 0]
    # s1 row 1 = [sqrt(2), -1, phi]
    # s2 row 2 = [0, phi, -1]
    s0 = (
        NEG_ONE, SQRT2, ZERO,
        ZERO, ONE, ZERO,
        ZERO, ZERO, ONE,
    )
    s1 = (
        ONE, ZERO, ZERO,
        SQRT2, NEG_ONE, PHI,
        ZERO, ZERO, ONE,
    )
    s2 = (
        ONE, ZERO, ZERO,
        ZERO, ONE, ZERO,
        ZERO, PHI, NEG_ONE,
    )
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


@dataclass(frozen=True)
class EnumerationStats:
    n: int
    count: int
    candidates: int
    seconds: float


class ExactPolyplet45Enumerator:
    def __init__(self, canonical_cache_enabled: bool = True) -> None:
        self.generators = coxeter_generators()
        self.face_group = unique_closure([self.generators[0], self.generators[1]])
        self.vertex_group = unique_closure([self.generators[1], self.generators[2]])
        self.key_to_id: dict[Key, int] = {}
        self.id_to_key: list[Key] = []
        self.coset_reps: list[Element] = []
        self.inverse_cache: dict[int, tuple[Element, ...]] = {}
        self.neighbor_cache: dict[int, tuple[int, ...]] = {}
        self.canonical_cache: dict[frozenset[int], tuple[int, ...]] = {}
        self.canonical_cache_enabled = canonical_cache_enabled
        self.base = self.coset_id(element_identity())
        self.neighbor_moves = self._neighbor_moves()
        self._check_local_model()

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
        moves: dict[int, Element] = {}
        for h in self.face_group:
            for k in self.vertex_group:
                move = h @ k
                cell_id = self.coset_id(move)
                if cell_id != self.base:
                    moves.setdefault(cell_id, move)
        return list(moves.values())

    def _check_local_model(self) -> None:
        if len(self.face_group) != 8:
            raise ValueError(f"Expected square stabilizer D4 of size 8, got {len(self.face_group)}")
        if len(self.vertex_group) != 10:
            raise ValueError(f"Expected vertex stabilizer D5 of size 10, got {len(self.vertex_group)}")
        if len(self.neighbor_moves) != 12:
            raise ValueError(f"Expected 12 polyplet neighbors, got {len(self.neighbor_moves)}")

    def neighbors(self, cell_id: int) -> list[int]:
        cached = self.neighbor_cache.get(cell_id)
        if cached is not None:
            return list(cached)
        rep = self.coset_reps[cell_id]
        result = tuple(self.coset_id(rep @ move) for move in self.neighbor_moves)
        self.neighbor_cache[cell_id] = result
        return list(result)

    def anchor_transforms(self, anchor: int) -> tuple[Element, ...]:
        cached = self.inverse_cache.get(anchor)
        if cached is not None:
            return cached
        anchor_rep = self.coset_reps[anchor]
        transforms = tuple(
            Element(
                mat_mul(h.inverse, anchor_rep.inverse),
                mat_mul(anchor_rep.matrix, h.matrix),
            )
            for h in self.face_group
        )
        self.inverse_cache[anchor] = transforms
        return transforms

    def canonical_shape(self, shape: frozenset[int]) -> tuple[int, ...]:
        cached = self.canonical_cache.get(shape)
        if self.canonical_cache_enabled and cached is not None:
            return cached
        best: tuple[int, ...] | None = None
        for anchor in shape:
            for transform in self.anchor_transforms(anchor):
                normalized = tuple(
                    sorted(self.coset_id(transform @ self.coset_reps[cell]) for cell in shape)
                )
                if best is None or normalized < best:
                    best = normalized
        assert best is not None
        if self.canonical_cache_enabled:
            self.canonical_cache[shape] = best
        return best

    def checkpoint(self, checkpoint_dir: str, n: int, shapes: set[tuple[int, ...]], counts: dict[int, int], stats: list[EnumerationStats]) -> None:
        os.makedirs(checkpoint_dir, exist_ok=True)
        path = os.path.join(checkpoint_dir, f"frontier_n{n}.pkl")
        with open(path, "wb") as handle:
            pickle.dump(
                {
                    "enumerator": self,
                    "n": n,
                    "shapes": shapes,
                    "counts": counts,
                    "stats": stats,
                },
                handle,
                protocol=pickle.HIGHEST_PROTOCOL,
            )

    def enumerate(
        self,
        max_n: int,
        progress: bool = False,
        checkpoint_dir: str = "",
        start_n: int = 1,
        shapes: set[tuple[int, ...]] | None = None,
        counts: dict[int, int] | None = None,
        stats: list[EnumerationStats] | None = None,
    ) -> tuple[dict[int, int], list[EnumerationStats]]:
        if shapes is None:
            shapes = {(self.base,)}
        if counts is None:
            counts = {1: 1}
        if stats is None:
            stats = [EnumerationStats(1, 1, 1, 0.0)]
        if checkpoint_dir and start_n == 1:
            self.checkpoint(checkpoint_dir, 1, shapes, counts, stats)

        for n in range(start_n + 1, max_n + 1):
            start = perf_counter()
            next_shapes: set[tuple[int, ...]] = set()
            candidates = 0
            for index, shape_tuple in enumerate(shapes, 1):
                if progress and index % 10000 == 0:
                    print(f"progress n={n} shapes={index}/{len(shapes)} next={len(next_shapes)}", flush=True)
                shape = frozenset(shape_tuple)
                boundary = set()
                for cell in shape:
                    boundary.update(self.neighbors(cell))
                boundary.difference_update(shape)
                for new_cell in boundary:
                    candidates += 1
                    expanded = frozenset((*shape, new_cell))
                    next_shapes.add(self.canonical_shape(expanded))
            elapsed = perf_counter() - start
            shapes = next_shapes
            counts[n] = len(shapes)
            stats.append(EnumerationStats(n, len(shapes), candidates, elapsed))
            if checkpoint_dir:
                self.checkpoint(checkpoint_dir, n, shapes, counts, stats)
        return counts, stats


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=8)
    parser.add_argument("--json", type=str, default="")
    parser.add_argument("--allow-mismatch", action="store_true")
    parser.add_argument("--progress", action="store_true")
    parser.add_argument("--checkpoint-dir", type=str, default="")
    parser.add_argument("--resume", type=str, default="")
    parser.add_argument("--no-canonical-cache", action="store_true")
    args = parser.parse_args()

    if args.resume:
        with open(args.resume, "rb") as handle:
            payload = pickle.load(handle)
        enumerator = payload["enumerator"]
        counts, stats = enumerator.enumerate(
            args.max_n,
            progress=args.progress,
            checkpoint_dir=args.checkpoint_dir,
            start_n=payload["n"],
            shapes=payload["shapes"],
            counts=payload["counts"],
            stats=payload["stats"],
        )
    else:
        enumerator = ExactPolyplet45Enumerator(
            canonical_cache_enabled=not args.no_canonical_cache
        )
        counts, stats = enumerator.enumerate(
            args.max_n,
            progress=args.progress,
            checkpoint_dir=args.checkpoint_dir,
        )

    print("n,count,known,status,candidates,seconds")
    for row in stats:
        known = KNOWN_A390200.get(row.n)
        status = "new" if known is None else ("ok" if known == row.count else "MISMATCH")
        print(f"{row.n},{row.count},{known if known is not None else ''},{status},{row.candidates},{row.seconds:.3f}")

    mismatches = {
        n: (counts[n], known)
        for n, known in KNOWN_A390200.items()
        if n <= args.max_n and counts.get(n) != known
    }
    if mismatches and not args.allow_mismatch:
        raise SystemExit(f"Known-prefix mismatch: {mismatches}")

    if args.json:
        payload = {
            "counts": counts,
            "known_prefix": KNOWN_A390200,
            "stats": [row.__dict__ for row in stats],
            "backend": "exact Q(sqrt(2), sqrt(5)) matrices",
        }
        with open(args.json, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)


if __name__ == "__main__":
    main()
