// Independent check of A390200 and A119611.
// This program shares no code or arithmetic with the Python enumerators in this repository.
// The tiling is built from real hyperbolic isometries in the hyperboloid model (double precision).
// Floating point is used only to build an integer table of cells, neighbor slots and frame changes.
// Every identification is checked against the gap between distinct cells, and the program stops if a check fails.
// After that, all counting is exact integer work:
//   grow:  level-by-level set of free shapes, with a canonical form equal to the minimum BFS mask code over all 8n frames.
//   redel: Redelmeier enumeration of rooted shapes. A shape is counted when its root frame gives the minimum code,
//          and the stabilizer sum divided by 8n is reported as a second count.
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>
#include <algorithm>

using namespace std;
typedef array<double, 9> M3;

static M3 mmul(const M3 &a, const M3 &b) {
    M3 c{};
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            double s = 0;
            for (int k = 0; k < 3; k++) s += a[3 * i + k] * b[3 * k + j];
            c[3 * i + j] = s;
        }
    return c;
}
static M3 minv(const M3 &a) { // Lorentz inverse J a^T J, J=diag(1,1,-1)
    M3 t{};
    const double J[3] = {1, 1, -1};
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) t[3 * i + j] = J[i] * a[3 * j + i] * J[j];
    return t;
}
static M3 rot(double th) { return M3{cos(th), -sin(th), 0, sin(th), cos(th), 0, 0, 0, 1}; }
static M3 boost(double d) { return M3{cosh(d), 0, sinh(d), 0, 1, 0, sinh(d), 0, cosh(d)}; }
static M3 ident() { return M3{1, 0, 0, 0, 1, 0, 0, 0, 1}; }
static double maxdiff(const M3 &a, const M3 &b) {
    double m = 0;
    for (int i = 0; i < 9; i++) m = max(m, fabs(a[i] - b[i]));
    return m;
}
// center point = column 2
static inline double px(const M3 &a) { return a[2]; }
static inline double py(const M3 &a) { return a[5]; }
static inline double pz(const M3 &a) { return a[8]; }
static double hdist_cosh(const M3 &a, const M3 &b) { // -<p,q>
    return pz(a) * pz(b) - px(a) * px(b) - py(a) * py(b);
}

// ---------------- global tables ----------------
static const int P = 4, Q = 5, G = 2 * P;
static int NS;                 // neighbor slots (12 or 4)
static vector<M3> D;           // D4 elements
static vector<M3> moves;       // neighbor moves m_s
static int gmul[G][G];
static uint8_t permT[G][16], kkT[G][16];
static vector<int32_t> nb;     // ncell*NS
static vector<uint8_t> tt;     // ncell*NS
static vector<uint8_t> layer;
static int ncell;
static double worstIdErr = 0, worstMatch = 0, bestNonMatch = 1e300;

static int identifyD(const M3 &m, double &err) {
    int best = -1;
    double bd = 1e300;
    for (int g = 0; g < G; g++) {
        double d = maxdiff(m, D[g]);
        if (d < bd) bd = d, best = g;
    }
    err = bd;
    return best;
}

static void buildLocal(bool plet) {
    double coshR = 1.0 / (tan(M_PI / P) * tan(M_PI / Q));
    double coshr = cos(M_PI / Q) / sin(M_PI / P);
    double R = acosh(coshR), r = acosh(coshr);
    M3 refl{1, 0, 0, 0, -1, 0, 0, 0, 1};
    D.clear();
    for (int e = 0; e < 2; e++)
        for (int i = 0; i < P; i++) D.push_back(e ? mmul(rot(2 * M_PI * i / P), refl) : rot(2 * M_PI * i / P));
    double err;
    for (int a = 0; a < G; a++)
        for (int b = 0; b < G; b++) {
            gmul[a][b] = identifyD(mmul(D[a], D[b]), err);
            if (err > 1e-9) { fprintf(stderr, "D4 not closed\n"); exit(1); }
        }
    vector<M3> cand;
    if (plet) {
        M3 A = mmul(mmul(mmul(mmul(rot(M_PI / P), boost(R)), rot(2 * M_PI / Q)), boost(-R)), rot(-M_PI / P));
        M3 Ak = ident();
        for (int k = 1; k <= Q; k++) {
            Ak = mmul(Ak, A);
            if (k == Q) {
                if (maxdiff(Ak, ident()) > 1e-9) { fprintf(stderr, "A^5 != I\n"); exit(1); }
                break;
            }
            for (int i = 0; i < P; i++) cand.push_back(mmul(rot(2 * M_PI * i / P), Ak));
        }
    } else {
        M3 H = mmul(mmul(boost(r), rot(M_PI)), boost(-r));
        for (int i = 0; i < P; i++) cand.push_back(mmul(rot(2 * M_PI * i / P), H));
    }
    moves.clear();
    for (auto &m : cand) {
        if (hdist_cosh(m, ident()) < 1.5) {
            if (hdist_cosh(m, ident()) > 1 + 1e-9 && hdist_cosh(m, ident()) < 1.3) { fprintf(stderr, "ambiguous\n"); exit(1); }
        }
        if (hdist_cosh(m, ident()) < 1 + 1e-9) continue; // F0 itself
        bool dup = false;
        for (auto &q : moves) if (hdist_cosh(m, q) < 1 + 1e-6) dup = true;
        if (!dup) moves.push_back(m);
    }
    sort(moves.begin(), moves.end(), [](const M3 &a, const M3 &b) { return atan2(py(a), px(a)) < atan2(py(b), px(b)); });
    NS = moves.size();
    printf("local model: %d neighbor slots; cosh distances:", NS);
    for (auto &m : moves) printf(" %.4f", hdist_cosh(m, ident()));
    printf("\n");
    for (int g = 0; g < G; g++)
        for (int s = 0; s < NS; s++) {
            M3 hm = mmul(D[g], moves[s]);
            int jp = -1;
            for (int j = 0; j < NS; j++) if (hdist_cosh(hm, moves[j]) < 1 + 1e-6) jp = j;
            if (jp < 0) { fprintf(stderr, "D4 does not permute neighbors\n"); exit(1); }
            permT[g][s] = jp;
            kkT[g][s] = identifyD(mmul(minv(moves[jp]), hm), err);
            if (err > 1e-9) { fprintf(stderr, "k not in D4\n"); exit(1); }
        }
}

struct KeyHash { size_t operator()(uint64_t k) const { return k * 0x9E3779B97F4A7C15ULL; } };

static void buildBall(int RMAX) {
    vector<M3> frame;
    unordered_multimap<uint64_t, int, KeyHash> grid;
    auto key = [](long long gx, long long gy) { return (uint64_t)(gx * 1000003LL) ^ (uint64_t)(gy + (1LL << 40)); };
    auto lookup = [&](const M3 &m) -> int {
        long long gx = llround(px(m) * 2), gy = llround(py(m) * 2);
        int found = -1;
        for (long long dx = -1; dx <= 1; dx++)
            for (long long dy = -1; dy <= 1; dy++) {
                auto rng = grid.equal_range(key(gx + dx, gy + dy));
                for (auto it = rng.first; it != rng.second; ++it) {
                    double c = hdist_cosh(m, frame[it->second]);
                    // same cell: c==1 ; distinct cells: c >= cosh(2r) ~ 1.618
                    if (c < 1.2) {
                        worstMatch = max(worstMatch, fabs(c - 1));
                        if (found >= 0 && found != it->second) { fprintf(stderr, "double match\n"); exit(1); }
                        found = it->second;
                    } else bestNonMatch = min(bestNonMatch, c);
                }
            }
        return found;
    };
    auto add = [&](const M3 &m, int L) {
        int id = frame.size();
        frame.push_back(m);
        layer.push_back(L);
        grid.emplace(key(llround(px(m) * 2), llround(py(m) * 2)), id);
        return id;
    };
    add(ident(), 0);
    size_t start = 0;
    for (int L = 0; L < RMAX; L++) {
        size_t end = frame.size();
        for (size_t x = start; x < end; x++)
            for (int s = 0; s < NS; s++) {
                M3 m = mmul(frame[x], moves[s]);
                if (lookup(m) < 0) add(m, L + 1);
            }
        printf("layer %d: %zu cells\n", L + 1, frame.size() - end);
        start = end;
    }
    ncell = frame.size();
    nb.assign((size_t)ncell * NS, -1);
    tt.assign((size_t)ncell * NS, 0);
    for (int x = 0; x < ncell; x++)
        for (int s = 0; s < NS; s++) {
            M3 m = mmul(frame[x], moves[s]);
            int y = lookup(m);
            if (y < 0) {
                if (layer[x] < RMAX) { fprintf(stderr, "missing interior neighbor\n"); exit(1); }
                continue;
            }
            nb[(size_t)x * NS + s] = y;
            double err;
            tt[(size_t)x * NS + s] = identifyD(mmul(minv(frame[y]), m), err);
            worstIdErr = max(worstIdErr, err);
        }
    // structural checks
    long bad = 0;
    for (int x = 0; x < ncell; x++) {
        if (layer[x] == RMAX) continue;
        for (int s = 0; s < NS; s++) {
            int y = nb[(size_t)x * NS + s];
            for (int u = 0; u < s; u++) if (nb[(size_t)x * NS + u] == y) bad++;
            if (y == x) bad++;
            bool back = false;
            for (int u = 0; u < NS; u++) if (nb[(size_t)y * NS + u] == x) back = true;
            if (!back) bad++;
        }
    }
    printf("ball radius %d: %d cells; structural defects=%ld; worst D4-identification err=%.3g; worst same-cell |cosh-1|=%.3g; min distinct-cell cosh=%.4f\n",
           RMAX, ncell, bad, worstIdErr, worstMatch, bestNonMatch);
    if (bad || worstIdErr > 1e-3 || worstMatch > 1e-3) { fprintf(stderr, "numerical check failed\n"); exit(1); }
}

// ---------------- canonical codes ----------------
// code(S, frame (a,g)) written to out[0..n-1]; compares against ref on the fly.
// returns: -1 if code < ref, 0 if equal, +1 if greater (early exit on first difference
// when stopAtDiff). If ref==nullptr, just computes.
static inline int frameCode(const int *cells, int n, int a, int g, const uint16_t *ref, uint16_t *out, bool stopAtDiff) {
    int lc[32];
    uint8_t lg[32];
    int cnt = 1;
    lc[0] = a;
    lg[0] = g;
    int cmp = 0;
    for (int i = 0; i < n; i++) {
        if (i >= cnt) { fprintf(stderr, "disconnected shape\n"); exit(1); }
        int x = lc[i], gg = lg[i];
        uint16_t mask = 0;
        for (int s = 0; s < NS; s++) {
            int sp = permT[gg][s];
            int y = nb[(size_t)x * NS + sp];
            if (y < 0) continue;
            bool in = false;
            for (int u = 0; u < n; u++) if (cells[u] == y) { in = true; break; }
            if (!in) continue;
            mask |= (uint16_t)(1u << s);
            bool lab = false;
            for (int u = 0; u < cnt; u++) if (lc[u] == y) { lab = true; break; }
            if (!lab) {
                lc[cnt] = y;
                lg[cnt] = gmul[tt[(size_t)x * NS + sp]][kkT[gg][s]];
                cnt++;
            }
        }
        if (out) out[i] = mask;
        if (ref && cmp == 0 && mask != ref[i]) {
            cmp = mask < ref[i] ? -1 : 1;
            if (stopAtDiff) return cmp;
        }
    }
    return cmp;
}

// reconstruct cells from code with root frame (cell 0, identity)
static int reconstruct(const uint16_t *code, int n, int *cells) {
    uint8_t lg[32];
    int cnt = 1;
    cells[0] = 0;
    lg[0] = 0;
    for (int i = 0; i < n; i++) {
        int x = cells[i], gg = lg[i];
        for (int s = 0; s < NS; s++) {
            if (!(code[i] >> s & 1)) continue;
            int sp = permT[gg][s];
            int y = nb[(size_t)x * NS + sp];
            if (y < 0) { fprintf(stderr, "reconstruct out of ball\n"); exit(1); }
            bool lab = false;
            for (int u = 0; u < cnt; u++) if (cells[u] == y) { lab = true; break; }
            if (!lab) {
                cells[cnt] = y;
                lg[cnt] = gmul[tt[(size_t)x * NS + sp]][kkT[gg][s]];
                cnt++;
            }
        }
    }
    if (cnt != n) { fprintf(stderr, "reconstruct size mismatch\n"); exit(1); }
    return cnt;
}

// canonical code: min over all frames; returns stabilizer size (#frames attaining min)
static int canonical(const int *cells, int n, uint16_t *best) {
    uint16_t tmp[32];
    int stab = 0;
    bool have = false;
    for (int ai = 0; ai < n; ai++)
        for (int g = 0; g < G; g++) {
            if (!have) {
                frameCode(cells, n, cells[ai], g, nullptr, best, false);
                have = true;
                stab = 1;
                continue;
            }
            int c = frameCode(cells, n, cells[ai], g, best, tmp, false);
            if (c < 0) { memcpy(best, tmp, n * sizeof(uint16_t)); stab = 1; }
            else if (c == 0) stab++;
        }
    return stab;
}

// ---------------- method 1: grow ----------------
static void runGrow(int nmax) {
    vector<string> cur;
    {
        int cells[1] = {0};
        uint16_t b[1];
        canonical(cells, 1, b);
        cur.push_back(string((char *)b, 2));
    }
    printf("grow n=1 free=1 rootedSum=%d\n", G);
    for (int n = 2; n <= nmax; n++) {
        unordered_map<string, int> next;
        next.reserve(cur.size() * 16);
        long long candidates = 0;
        for (auto &cs : cur) {
            int m = n - 1;
            int cells[32];
            reconstruct((const uint16_t *)cs.data(), m, cells);
            vector<int> bnd;
            for (int i = 0; i < m; i++)
                for (int s = 0; s < NS; s++) {
                    int y = nb[(size_t)cells[i] * NS + s];
                    if (y < 0) { fprintf(stderr, "grow out of ball\n"); exit(1); }
                    bool in = false;
                    for (int u = 0; u < m; u++) if (cells[u] == y) in = true;
                    if (!in) bnd.push_back(y);
                }
            sort(bnd.begin(), bnd.end());
            bnd.erase(unique(bnd.begin(), bnd.end()), bnd.end());
            for (int y : bnd) {
                candidates++;
                cells[m] = y;
                uint16_t best[32];
                int st = canonical(cells, n, best);
                next.emplace(string((char *)best, 2 * n), st);
            }
        }
        long long rootedSum = 0, stabSum = 0;
        for (auto &kv : next) {
            if ((G * n) % kv.second) { fprintf(stderr, "stab does not divide 8n\n"); exit(1); }
            rootedSum += (G * n) / kv.second;
            stabSum += kv.second;
        }
        cur.clear();
        cur.reserve(next.size());
        for (auto &kv : next) cur.push_back(kv.first);
        printf("grow n=%d free=%zu candidates=%lld rooted(=sum 8n/|Stab|)=%lld\n", n, cur.size(), candidates, rootedSum);
        fflush(stdout);
    }
}

// ---------------- method 2: Redelmeier rooted + minimal-root test ----------------
struct Counters { long long rooted[32]{}, freeC[32]{}, stabSum[32]{}; };
static int NMAX;

static inline void processRooted(const int *shape, int n, Counters &C) {
    uint16_t ref[32], tmp[32];
    frameCode(shape, n, shape[0], 0, nullptr, ref, false);
    bool minimal = true;
    int stab = 0;
    for (int ai = 0; ai < n; ai++)
        for (int g = 0; g < G; g++) {
            int c = frameCode(shape, n, shape[ai], g, ref, tmp, true);
            if (c < 0) minimal = false;
            else if (c == 0) stab++;
        }
    C.rooted[n]++;
    C.stabSum[n] += stab;
    if (minimal) C.freeC[n]++;
}

static void redel(int *shape, int depth, const int *untried, int ulen, vector<uint8_t> &forb, Counters &C) {
    processRooted(shape, depth, C);
    if (depth == NMAX) return;
    int buf[512];
    for (int i = 0; i < ulen; i++) {
        int v = untried[i];
        shape[depth] = v;
        int k = 0;
        for (int j = i + 1; j < ulen; j++) buf[k++] = untried[j];
        int added[16], na = 0;
        for (int s = 0; s < NS; s++) {
            int y = nb[(size_t)v * NS + s];
            if (y < 0 || forb[y]) continue;
            forb[y] = 1;
            added[na++] = y;
            buf[k++] = y;
        }
        redel(shape, depth + 1, buf, k, forb, C);
        for (int u = 0; u < na; u++) forb[added[u]] = 0;
    }
}

struct Task { vector<int> shape, untried, forbList; };

static void collect(vector<int> &shape, const vector<int> &untried, vector<int> &forbList, vector<uint8_t> &forb, int splitDepth, vector<Task> &tasks, Counters &C) {
    int depth = shape.size();
    if (depth == splitDepth || depth == NMAX) {
        tasks.push_back({shape, untried, forbList});
        return;
    }
    processRooted(shape.data(), depth, C);
    for (size_t i = 0; i < untried.size(); i++) {
        int v = untried[i];
        vector<int> nu(untried.begin() + i + 1, untried.end());
        size_t fl = forbList.size();
        for (int s = 0; s < NS; s++) {
            int y = nb[(size_t)v * NS + s];
            if (y < 0 || forb[y]) continue;
            forb[y] = 1;
            forbList.push_back(y);
            nu.push_back(y);
        }
        shape.push_back(v);
        collect(shape, nu, forbList, forb, splitDepth, tasks, C);
        shape.pop_back();
        for (size_t u = fl; u < forbList.size(); u++) forb[forbList[u]] = 0;
        forbList.resize(fl);
    }
}

static void runRedel(int nmax, int threads) {
    NMAX = nmax;
    vector<uint8_t> forb(ncell, 0);
    Counters C0;
    vector<int> shape{0}, untried, forbList{0};
    forb[0] = 1;
    for (int s = 0; s < NS; s++) {
        int y = nb[s];
        if (!forb[y]) { forb[y] = 1; untried.push_back(y); forbList.push_back(y); }
    }
    vector<Task> tasks;
    collect(shape, untried, forbList, forb, min(nmax, 4), tasks, C0);
    printf("redel: %zu tasks, %d threads\n", tasks.size(), threads);
    fflush(stdout);
    atomic<size_t> nextT{0};
    vector<Counters> TC(threads);
    vector<thread> th;
    for (int t = 0; t < threads; t++)
        th.emplace_back([&, t]() {
            vector<uint8_t> f(ncell, 0);
            int sh[64];
            while (true) {
                size_t i = nextT++;
                if (i >= tasks.size()) break;
                Task &T = tasks[i];
                for (int c : T.forbList) f[c] = 1;
                for (size_t u = 0; u < T.shape.size(); u++) sh[u] = T.shape[u];
                redel(sh, T.shape.size(), T.untried.data(), T.untried.size(), f, TC[t]);
                for (int c : T.forbList) f[c] = 0;
            }
        });
    for (auto &x : th) x.join();
    for (auto &c : TC)
        for (int n = 0; n < 32; n++) C0.rooted[n] += c.rooted[n], C0.freeC[n] += c.freeC[n], C0.stabSum[n] += c.stabSum[n];
    for (int n = 1; n <= nmax; n++) {
        long long S = C0.stabSum[n];
        printf("redel n=%d rooted=%lld free(minimal-root)=%lld sumStab=%lld sumStab/(8n)=%lld rem=%lld\n", n, C0.rooted[n], C0.freeC[n], S, S / (G * n), S % (G * n));
    }
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: enum45 plet|omino grow|redel|both nmax [threads]\n"); return 1; }
    bool plet = string(argv[1]) == "plet";
    string method = argv[2];
    int nmax = atoi(argv[3]);
    int threads = argc > 4 ? atoi(argv[4]) : 8;
    buildLocal(plet);
    buildBall(max(1, nmax - 1));
    if (method == "grow" || method == "both") runGrow(nmax);
    if (method == "redel" || method == "both") runRedel(nmax, threads);
    return 0;
}
