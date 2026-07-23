# MIHNP 3,4,5-sample experiment for SageMath 10.8+.
# Small-root solving + Assumption 1 verification.
from sage.all import *
import itertools
import time

# ============================================================
# User configuration
# ============================================================
#
# Table 2 parameter presets:
# (NSAMPLES, PBITS, UBITS, M, TH)
#
# 3 samples:
#   (3, 512, 195, 2, QQ(3))
#   (3, 256, 102, 3, QQ(3))
#   (3, 256, 105, 4, QQ(3))
#
# 4 samples:
#   (4, 256, 101, 1, QQ(31521)/QQ(10000))
#   (4, 512, 230, 2, QQ(31521)/QQ(10000))
#
# 5 samples:
#   (5, 512, 225, 1, QQ(34028)/QQ(10000))
#
# To reproduce a particular row, manually replace NSAMPLES,
# PBITS, UBITS, M, and TH below.
#

N_RUNS = 1        # Number of independent experiments.
NSAMPLES = 3  # Number of MIHNP samples.
PBITS = 512 # Modulus bits.
UBITS = 195 # Unknown low bits of every y_i.
M = 2         # Polytope scale.
TH = QQ(3)       # sum(beta_i) <= M*TH.
DELTA = 0.99      # LLL parameter.

Q0 = 101          # First Groebner prime.
QBAD = 100        # Failed-prime limit.
QMAX = 500        # Tested-prime limit.
SEED = None       # Integer for a repeatable run sequence.

# Assumption 1 controls.
# Exact verification can be expensive, especially for five samples.
VERIFY_ASSUMPTION1 = True
USE_MSOLVE_ASSUMPTION1 = False


# Ring order: x_n > ... > x_1.  All exponent tuples used below are
# nevertheless stored in canonical order (e_1,...,e_n).
R = PolynomialRing(
    ZZ,
    names=tuple(
        "x%d" % i
        for i in range(NSAMPLES, 0, -1)
    ),
    order="lex"
)

RGENS = R.gens()                 # (x_n,...,x_1)
VARS = tuple(reversed(RGENS))    # (x_1,...,x_n)


def cmod(a, modulus):
    """Centered residue modulo modulus."""
    a = ZZ(a) % modulus
    return a - modulus if a > modulus//2 else a


def rnd_prime(n):
    """Random n-bit probable prime."""
    return ZZ(
        random_prime(
            2**n - 1,
            lbound=2**(n - 1),
            proof=False
        )
    )


def sample(p):
    """Generate one fresh NSAMPLES-sample MIHNP instance."""
    F = GF(p)
    alpha = F.random_element()
    t = []
    y = []

    while len(t) < NSAMPLES:
        ti = F.random_element()

        if ti in t or alpha + ti == 0:
            continue

        t.append(ti)
        y.append(1/(alpha + ti))

    alpha = ZZ(alpha)
    t = [ZZ(z) for z in t]
    y = [ZZ(z) for z in y]

    for ti, yi in zip(t, y):
        if yi*(alpha + ti) % p != 1:
            raise RuntimeError(
                "invalid generated MIHNP sample"
            )

    return alpha, t, y


def ex(e):
    """Convert ring exponents (e_n,...,e_1) to (e_1,...,e_n)."""
    return tuple(
        reversed(
            tuple(map(int, e))
        )
    )


def mon(e):
    """Build a monomial from canonical exponents."""
    g = R(1)

    for x, a in zip(VARS, e):
        g *= x**int(a)

    return g


def ev(f, u):
    """Evaluate at canonical root tuple u=(u_1,...,u_n)."""
    return f(
        *tuple(reversed(tuple(u)))
    )


def center_poly(f, modulus, force=None):
    """
    Center coefficients modulo modulus.

    If force is supplied, that canonical monomial coefficient must be
    1 modulo modulus and is lifted exactly as 1.
    """
    g = R(0)
    seen = False

    for ering, a in f.dict().items():
        e = ex(ering)

        if force is not None and e == force:
            if (ZZ(a) - 1) % modulus:
                raise ValueError(
                    "forced coefficient is not 1 modulo the modulus"
                )

            a = ZZ(1)
            seen = True

        else:
            a = cmod(a, modulus)

        if a:
            g += a*mon(e)

    if force is not None and not seen:
        raise ValueError(
            "forced leading exponent is absent"
        )

    return g


def supp(f):
    """Return the canonical exponent support."""
    return {
        ex(e)
        for e, a in f.dict().items()
        if a
    }


def lead(f):
    """Return the leading exponent in canonical order."""
    if not f:
        raise RuntimeError(
            "zero polynomial has no leading monomial"
        )

    return ex(
        next(
            iter(
                f.lm().dict()
            )
        )
    )


def prim(f):
    """Remove integer content and fix the leading sign."""
    if not f:
        raise RuntimeError(
            "cannot primitive-normalize zero"
        )

    d = abs(
        ZZ(
            gcd(
                [
                    ZZ(a)
                    for a in f.dict().values()
                ]
            )
        )
    )

    if not d:
        raise RuntimeError(
            "zero polynomial content"
        )

    g = R(0)

    for e, a in f.dict().items():
        g += (ZZ(a)//d)*mon(ex(e))

    if ZZ(g.monomial_coefficient(g.lm())) < 0:
        g = -g

    return g


def idx_exp(indices):
    """Squarefree canonical exponent vector on the given indices."""
    e = [0]*NSAMPLES

    for i in indices:
        e[i] = 1

    return tuple(e)


def pair_poly(i, j, known, t, p):
    """Construct monic f_ij(x_i,x_j) = 0 mod p."""
    Yi = known[i] + VARS[i]
    Yj = known[j] + VARS[j]
    d = ZZ(t[i] - t[j]) % p

    if not d:
        raise ValueError(
            "duplicate t_i modulo p"
        )

    f = (
        Yi*Yj
        + ZZ(inverse_mod(d, p))*(Yi - Yj)
    )

    return center_poly(f, p)


def extra_poly(pair, subset, p):
    """
    Construct one normalized p^2 polynomial for a four-index subset.

    For subset (a,b,c,d), product differences cancel the degree-four
    term.  The coefficient of x_b*x_c*x_d is normalized to 1.
    """
    a, b, c, d = subset
    modulus = p**2
    target = idx_exp((b, c, d))

    f_ab = pair[(a, b)]
    f_ac = pair[(a, c)]
    f_ad = pair[(a, d)]
    f_bc = pair[(b, c)]
    f_bd = pair[(b, d)]
    f_cd = pair[(c, d)]

    candidates = [
        f_ab*f_cd - f_ac*f_bd,
        f_ab*f_cd - f_ad*f_bc,
        f_ac*f_bd - f_ad*f_bc,
    ]
    candidates += [-g for g in candidates]

    for g in candidates:
        lc = ZZ(
            g.monomial_coefficient(
                mon(target)
            )
        )

        if not lc or gcd(lc, p) != 1:
            continue

        h = ZZ(
            inverse_mod(
                lc % modulus,
                modulus
            )
        )*g

        h = center_poly(
            h,
            modulus,
            force=target
        )

        if lead(h) != target:
            raise RuntimeError(
                "unexpected extra-polynomial leading monomial"
            )

        return h

    raise ValueError(
        "could not normalize extra polynomial for subset %s"
        % (subset,)
    )


def system(known, t, p):
    """Build pair polynomials and all four-subset p^2 polynomials."""
    pair = {}
    names = []
    P = []
    A = []
    E = []

    for i, j in itertools.combinations(
        range(NSAMPLES),
        2
    ):
        f = pair_poly(
            i,
            j,
            known,
            t,
            p
        )

        pair[(i, j)] = f
        names.append(
            "f%d%d" % (i + 1, j + 1)
        )
        P.append(f)
        A.append(idx_exp((i, j)))
        E.append(1)

    for subset in itertools.combinations(
        range(NSAMPLES),
        4
    ):
        f = extra_poly(
            pair,
            subset,
            p
        )

        names.append(
            "f" + "".join(
                str(i + 1)
                for i in subset
            )
        )
        P.append(f)

        _, b, c, d = subset
        A.append(idx_exp((b, c, d)))
        E.append(2)

    return names, P, A, E


def inside(b):
    """Test membership in the scaled symmetric MIHNP polytope."""
    return (
        all(
            0 <= z <= M
            for z in b
        )
        and QQ(sum(b)) <= M*TH
    )


def mons():
    """Enumerate lattice monomials."""
    B = [
        tuple(map(int, b))
        for b in itertools.product(
            range(M + 1),
            repeat=NSAMPLES
        )
        if inside(b)
    ]

    return sorted(
        B,
        key=lambda b: tuple(reversed(b))
    )


def ell_vectors(b, A):
    """Enumerate ell with sum_j ell_j*A_j <= b using DFS pruning."""
    current = [0]*len(A)
    used = [0]*NSAMPLES

    def rec(j):
        if j == len(A):
            yield tuple(current)
            return

        alpha = A[j]
        bounds = []

        for i, z in enumerate(alpha):
            if z:
                remain = b[i] - used[i]

                if remain < 0:
                    return

                bounds.append(
                    remain//z
                )

        ub = min(bounds) if bounds else 0

        for ell in range(ub + 1):
            current[j] = ell

            for i in range(NSAMPLES):
                used[i] += ell*alpha[i]

            yield from rec(j + 1)

            for i in range(NSAMPLES):
                used[i] -= ell*alpha[i]

        current[j] = 0

    yield from rec(0)


def pick(b, S, P, A, E):
    """Choose the admissible shift with largest weighted p-power."""
    best = None

    for ell in ell_vectors(b, A):
        used = tuple(
            sum(
                ell[j]*A[j][i]
                for j in range(len(A))
            )
            for i in range(NSAMPLES)
        )

        a = tuple(
            b[i] - used[i]
            for i in range(NSAMPLES)
        )

        g = mon(a)

        for z, f in zip(ell, P):
            if z:
                g *= f**z

        if (
            supp(g).issubset(S)
            and lead(g) == b
        ):
            r = sum(
                E[j]*ell[j]
                for j in range(len(E))
            )
            key = (r,) + tuple(ell)

            if best is None or key > best[0]:
                best = (key, g)

    if best is None:
        raise RuntimeError(
            "no admissible shift for monomial %s"
            % (b,)
        )

    return best[0][0], best[1]


def wt(e, X):
    """Common-bound monomial scaling weight."""
    return X**sum(e)


def make_lat(B, G, p, t, X):
    """Build the scaled coefficient lattice."""
    pos = {
        b: i
        for i, b in enumerate(B)
    }
    rows = []

    for b in B:
        r, g = G[b]
        g *= p**(t - r)
        row = [ZZ(0)]*len(B)

        for e, a in g.dict().items():
            e = ex(e)

            if e not in pos:
                raise RuntimeError(
                    "support outside monomial set"
                )

            row[pos[e]] += ZZ(a)*wt(e, X)

        rows.append(row)

    return Matrix(ZZ, rows)


def to_poly(v, B, X):
    """Recover one primitive polynomial from an LLL row."""
    f = R(0)

    for a, b in zip(v, B):
        q, r = ZZ(a).quo_rem(
            wt(b, X)
        )

        if r:
            raise RuntimeError(
                "invalid lattice scaling"
            )

        if q:
            f += q*mon(b)

    if not f:
        raise RuntimeError(
            "zero LLL polynomial"
        )

    return prim(f)


def lll_polys(L, B, X):
    """Keep the LLL row order H1,H2,... ."""
    return [
        to_poly(
            L.row(i),
            B,
            X
        )
        for i in range(L.nrows())
    ]


def _gb_is_zero_dim(G, n):
    """Decide zero-dimensionality from Groebner leading monomials."""
    pure_power_seen = [False]*n

    for g in G:
        if not g:
            continue

        e = tuple(
            next(
                iter(
                    g.lm().dict()
                )
            )
        )
        nonzero = [
            i
            for i, a in enumerate(e)
            if a
        ]

        if len(nonzero) == 1:
            pure_power_seen[
                nonzero[0]
            ] = True

    return all(pure_power_seen)


def assumption1(H):
    """Check whether the first NSAMPLES auxiliary polynomials are zero-dimensional."""
    if not VERIFY_ASSUMPTION1:
        return None

    if len(H) < NSAMPLES:
        return False

    Q = PolynomialRing(
        QQ,
        names=R.variable_names(),
        order="degrevlex"
    )
    I = Q.ideal(
        [
            Q(f)
            for f in H[:NSAMPLES]
        ]
    )

    if USE_MSOLVE_ASSUMPTION1:
        try:
            G = I.groebner_basis(
                algorithm="msolve",
                proof=False
            )

        except Exception as exc:
            raise RuntimeError(
                "msolve is unavailable or failed; install it for this "
                "SageMath environment, or set "
                "USE_MSOLVE_ASSUMPTION1 = False"
            ) from exc

    else:
        G = I.groebner_basis(
            algorithm="libsingular:slimgb"
        )

    return _gb_is_zero_dim(
        G,
        NSAMPLES
    )


def zero_polys(H, u):
    """Require H1,...,H_n to vanish; return every vanishing H_i."""
    if len(H) < NSAMPLES:
        return False, []

    is_zero = [
        ZZ(ev(f, u)) == 0
        for f in H
    ]
    ok = all(
        is_zero[:NSAMPLES]
    )
    Z = [
        f
        for f, flag in zip(H, is_zero)
        if flag
    ]

    return ok, Z if ok else []


def lift(a, q):
    """Centered CRT lift."""
    modulus = prod(q)
    z = ZZ(crt(a, q))
    return z - modulus if z > modulus//2 else z


def valid_mihnp_root(u, known, t, p):
    """Check that all reconstructed y_i give the same alpha."""
    alpha = []

    for ui, yi0, ti in zip(u, known, t):
        yi = ZZ(yi0 + ui) % p

        if not yi:
            return False

        alpha.append(
            ZZ(inverse_mod(yi, p) - ti) % p
        )

    return all(
        a == alpha[0]
        for a in alpha[1:]
    )


def solve(H, X, P, E, p, known, t):
    """Solve the selected auxiliary system modulo small primes and CRT."""
    if len(H) < NSAMPLES:
        return None, 0.0

    residues = [
        []
        for _ in range(NSAMPLES)
    ]  # Ring order x_n,...,x_1.
    moduli = []
    modulus_product = ZZ(1)
    q = ZZ(Q0 - 1)
    bad = 0
    tried = 0
    gb_time = 0.0

    while modulus_product <= 2*X:
        if bad >= QBAD or tried >= QMAX:
            return None, gb_time

        q = ZZ(next_prime(q))
        tried += 1
        F = R.change_ring(GF(q))
        Pq = []

        for f in H:
            fq = F(f)

            if fq:
                Pq.append(fq)

        if len(Pq) < NSAMPLES:
            bad += 1
            continue

        try:
            I = F.ideal(Pq)
            tic = time.perf_counter()

            try:
                I.groebner_basis()

            finally:
                gb_time += (
                    time.perf_counter()
                    - tic
                )

            V = I.variety()

        except Exception:
            bad += 1
            continue

        if len(V) != 1:
            bad += 1
            continue

        try:
            a = [
                ZZ(V[0][v])
                for v in F.gens()
            ]

        except KeyError:
            bad += 1
            continue

        moduli.append(q)

        for i in range(NSAMPLES):
            residues[i].append(a[i])

        modulus_product *= q

    ring_root = [
        lift(
            residues[i],
            moduli
        )
        for i in range(NSAMPLES)
    ]
    u = tuple(reversed(ring_root))

    if any(
        z < 0 or z >= X
        for z in u
    ):
        return None, gb_time

    if any(
        ZZ(ev(f, u)) != 0
        for f in H
    ):
        return None, gb_time

    if any(
        ZZ(ev(f, u)) % (p**e)
        for f, e in zip(P, E)
    ):
        return None, gb_time

    if not valid_mihnp_root(
        u,
        known,
        t,
        p
    ):
        return None, gb_time

    return u, gb_time


def run_once(run_id):
    """Generate and run one fresh experiment."""
    for _attempt in range(100):
        p = rnd_prime(PBITS)
        _, t, y = sample(p)
        X = ZZ(2)**UBITS
        u = tuple(
            ZZ(yi) % X
            for yi in y
        )
        known = tuple(
            ZZ(yi) - ui
            for yi, ui in zip(y, u)
        )

        try:
            names, P, A, E = system(
                known,
                t,
                p
            )

        except ValueError:
            continue

        break

    else:
        raise RuntimeError(
            "could not generate a normalizable MIHNP instance"
        )

    if not valid_mihnp_root(
        u,
        known,
        t,
        p
    ):
        raise RuntimeError(
            "invalid generated MIHNP root"
        )

    for name, f, alpha, e in zip(
        names,
        P,
        A,
        E
    ):
        if lead(f) != alpha:
            raise RuntimeError(
                "unexpected leading monomial for %s"
                % name
            )

        if ZZ(ev(f, u)) % (p**e):
            raise RuntimeError(
                "invalid congruence for %s"
                % name
            )

    B = mons()
    S = set(B)
    G = {
        b: pick(
            b,
            S,
            P,
            A,
            E
        )
        for b in B
    }
    t_power = max(
        r
        for r, _ in G.values()
    )
    L = make_lat(
        B,
        G,
        p,
        t_power,
        X
    )

    if L.rank() != L.nrows():
        raise RuntimeError(
            "rank-deficient lattice"
        )

    tic = time.perf_counter()

    try:
        LR = L.LLL(delta=DELTA)

    except TypeError:
        LR = L.LLL()

    lll_time = time.perf_counter() - tic
    H = lll_polys(LR, B, X)

    if VERIFY_ASSUMPTION1:
        tic = time.perf_counter()
        a1 = assumption1(H)
        a1_time = time.perf_counter() - tic

    else:
        a1 = None
        a1_time = None

    ok, Z = zero_polys(H, u)
    rec = None
    gb_time = 0.0

    if ok:
        rec, gb_time = solve(
            Z,
            X,
            P,
            E,
            p,
            known,
            t
        )

    success = rec == u
    output = [
        "modulus bits=%d" % p.nbits(),
        "unknown bits=%d" % UBITS,
        "lattice dimension=%d" % L.nrows(),
        "LLL time seconds=%.6f" % lll_time,
        "Groebner time seconds=%.6f" % gb_time,
        "root recovered=%s" % bool(success),
    ]

    if VERIFY_ASSUMPTION1:
        output.extend(
            [
                "Assumption 1 check time seconds=%.6f"
                % a1_time,
                "Assumption 1 holds=%s"
                % bool(a1),
            ]
        )

    print(" | ".join(output))
    return bool(success)


def main():
    """Run N_RUNS fresh experiments."""
    if NSAMPLES not in (3, 4, 5):
        raise ValueError(
            "require NSAMPLES in {3,4,5}"
        )

    if not 0 < UBITS < PBITS:
        raise ValueError(
            "require 0 < UBITS < PBITS"
        )

    if ZZ(M) != M or M < 1:
        raise ValueError(
            "require M to be a positive integer"
        )

    if ZZ(N_RUNS) != N_RUNS or N_RUNS < 1:
        raise ValueError(
            "require N_RUNS to be a positive integer"
        )

    if SEED is not None:
        set_random_seed(SEED)

    all_success = True

    for run_id in range(
        1,
        int(N_RUNS) + 1
    ):
        all_success = (
            run_once(run_id)
            and all_success
        )

    return all_success


if __name__ == "__main__":
    if not main():
        raise SystemExit(1)
