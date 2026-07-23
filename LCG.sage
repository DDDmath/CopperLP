# LCG unknown-multiplier experiment for SageMath 10.8+.
# Small-root solving + Assumption 1 verification.
from sage.all import *
import time

# ============================================================
# User configuration
# ============================================================
#
# Table 2 parameter presets:
# (NSAMPLES, PBITS, UBITS, M)
#
#   (6, 256, 72, 2)
#   (8, 256, 75, 3)
#
# NVAR is automatically set equal to NSAMPLES.
#
# To reproduce a particular row, manually replace NSAMPLES,
# PBITS, UBITS, and M below.
#

N_RUNS = 1        # Number of independent experiments.
NSAMPLES = 6      # Number of consecutive multiplicative-LCG outputs.
NVAR = NSAMPLES   # One unknown low part for every output.
PBITS = 256       # Prime modulus bits.
UBITS = 72        # Unknown low bits of every output.
M = 2             # Standard-simplex polytope scale.
DELTA = 0.99      # LLL parameter.

Q0 = 101          # First Groebner prime.
QBAD = 100        # Failed-prime limit.
QMAX = 500        # Tested-prime limit.
SEED = None       # Integer for a repeatable run sequence.

# Assumption 1 controls.
# Exact verification can be expensive for six variables.  Testing only
# small M is recommended; msolve can be faster but must be installed
# separately for the same SageMath installation.
VERIFY_ASSUMPTION1 = True
USE_MSOLVE_ASSUMPTION1 = False


# The LCG relations are
#
#     f_i = (u_{i+1}+x_{i+1})^2
#           - (u_i+x_i)(u_{i+2}+x_{i+2}) == 0 mod p,
#
# for i=1,...,NSAMPLES-2.  We need LM(f_i)=x_{i+1}^2 for all i.
# No ordinary global lex variable order can make every middle variable
# larger than both of its neighbours simultaneously.  Instead use a
# weighted-degree lexicographic order with a strictly concave positive
# weight sequence.  The weights below satisfy
#
#     2*w_{i+1} > w_i + w_{i+2}
#
# for every relation, so x_{i+1}^2 has strictly larger weighted degree
# than x_i*x_{i+2}.
LCG_WEIGHTS = tuple(
    (NSAMPLES + 1)**2 - i**2
    for i in range(NSAMPLES)
)
LCG_ORDER = TermOrder(
    "wdeglex",
    LCG_WEIGHTS
)
VAR_NAMES = tuple(
    "x%d" % (i + 1)
    for i in range(NVAR)
)
R = PolynomialRing(
    ZZ,
    names=VAR_NAMES,
    order=LCG_ORDER
)
VARS = R.gens()   # (x1,x2,...,x_NSAMPLES)


# ============================================================
# Basic polynomial utilities
# ============================================================

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


def split_low(z, X):
    """Write z = known + low with 0 <= low < X."""
    z = ZZ(z)
    low = z % X
    return z - low, low


def mon(e):
    """Build x1^e1*...*x_n^e_n from a canonical exponent tuple."""
    g = R(1)

    for x, a in zip(VARS, e):
        g *= x**int(a)

    return g


def ev(f, u):
    """Evaluate f at u=(u1,...,u_n)."""
    return f(*tuple(u))


def center_poly(f, modulus):
    """Center every coefficient of f modulo modulus."""
    g = R(0)

    for e, a in f.dict().items():
        a = cmod(a, modulus)

        if a:
            g += a*mon(tuple(map(int, e)))

    return g


def supp(f):
    """Return the exponent support in (e1,...,e_n) order."""
    return {
        tuple(map(int, e))
        for e, a in f.dict().items()
        if a
    }


def lead(f):
    """Return the leading exponent in (e1,...,e_n) order."""
    if not f:
        raise RuntimeError(
            "zero polynomial has no leading monomial"
        )

    return tuple(
        map(
            int,
            next(
                iter(
                    f.lm().dict()
                )
            )
        )
    )


def prim(f):
    """Remove integer content and fix the leading sign."""
    if not f:
        raise RuntimeError(
            "cannot primitive-normalize zero"
        )

    d = ZZ(0)

    for a in f.dict().values():
        d = gcd(
            d,
            abs(ZZ(a))
        )

    if not d:
        raise RuntimeError(
            "zero polynomial content"
        )

    g = R(0)

    for e, a in f.dict().items():
        g += (ZZ(a)//d)*mon(
            tuple(map(int, e))
        )

    if ZZ(
        g.monomial_coefficient(
            g.lm()
        )
    ) < 0:
        g = -g

    return g


def transport(f, S):
    """
    Copy f into another polynomial ring S with the same variable order.

    This explicit map is used because the lattice ring has a weighted
    order, while Assumption 1 uses degrevlex and modular solving uses lex.
    """
    if S.ngens() != NVAR:
        raise ValueError(
            "target ring has the wrong number of variables"
        )

    g = S(0)
    y = S.gens()
    K = S.base_ring()

    for e, a in f.dict().items():
        term = K(a)

        for yi, z in zip(y, e):
            if z:
                term *= yi**int(z)

        g += term

    return g


# ============================================================
# Multiplicative LCG instance generation
# ============================================================

def sample(p):
    """
    Generate a nondegenerate multiplicative LCG instance

        s_{i+1} = multiplier*s_i mod p.

    The multiplier and seed are nonzero, multiplier != 1, and the sampled
    outputs are required to be distinct.  The latter only removes small-
    order degenerate test instances; it does not change the equations.
    """
    for _ in range(1000):
        multiplier = ZZ.random_element(
            2,
            p
        )
        seed = ZZ.random_element(
            1,
            p
        )
        seq = []
        cur = seed

        for _ in range(NSAMPLES):
            seq.append(ZZ(cur))
            cur = ZZ(multiplier*cur % p)

        if len(set(seq)) != NSAMPLES:
            continue

        if any(
            (seq[i + 1]**2 - seq[i]*seq[i + 2]) % p
            for i in range(NSAMPLES - 2)
        ):
            raise RuntimeError(
                "invalid generated LCG sequence"
            )

        return (
            ZZ(multiplier),
            ZZ(seed),
            tuple(seq),
        )

    raise RuntimeError(
        "could not generate a nondegenerate LCG instance"
    )


# ============================================================
# LCG polynomial system
# ============================================================

def relation_lead(i):
    """Leading exponent of relation f_{i+1}, with 0-based i."""
    if not 0 <= i < NSAMPLES - 2:
        raise ValueError(
            "invalid LCG relation index"
        )

    e = [0]*NVAR
    e[i + 1] = 2
    return tuple(e)


def system(known, p):
    """
    Build f_1,...,f_{NSAMPLES-2} and their leading data.

    Here s_i=known_i+x_i and

        f_i = s_{i+1}^2 - s_i*s_{i+2} == 0 mod p.
    """
    if len(known) != NSAMPLES:
        raise ValueError(
            "wrong number of known LCG output parts"
        )

    Y = [
        ZZ(known[i]) + VARS[i]
        for i in range(NSAMPLES)
    ]
    names = []
    polys = []
    leading = []
    powers = []

    for i in range(NSAMPLES - 2):
        f = center_poly(
            Y[i + 1]**2
            - Y[i]*Y[i + 2],
            p
        )
        names.append(
            "f%d" % (i + 1)
        )
        polys.append(f)
        leading.append(
            relation_lead(i)
        )
        powers.append(1)

    return (
        names,
        polys,
        leading,
        powers,
    )


def valid_lcg_root(u, X, known, p):
    """
    Validate a candidate using the public LCG structure and root bounds.

    The reconstructed outputs must be nonzero residues below p and all
    consecutive ratios must give one common multiplier modulo p.
    """
    if len(u) != NSAMPLES:
        return False

    if any(
        z < 0 or z >= X
        for z in u
    ):
        return False

    seq = tuple(
        ZZ(known[i] + u[i])
        for i in range(NSAMPLES)
    )

    if any(
        s <= 0 or s >= p
        for s in seq
    ):
        return False

    multiplier = ZZ(
        seq[1]*inverse_mod(
            seq[0],
            p
        )
    ) % p

    return all(
        seq[i + 1]
        == ZZ(multiplier*seq[i] % p)
        for i in range(NSAMPLES - 1)
    )


# ============================================================
# Standard-simplex polytope and admissible shifts
# ============================================================

def inside(b):
    """Test b in M*{z_i>=0, sum(z_i)<=1}."""
    return (
        len(b) == NVAR
        and min(b) >= 0
        and sum(b) <= M
    )


def monomial_key(b):
    """Ascending key compatible with the weighted-degree lex order."""
    return (
        sum(
            LCG_WEIGHTS[i]*b[i]
            for i in range(NVAR)
        ),
        tuple(b),
    )


def mons():
    """Enumerate every integer monomial exponent in the scaled simplex."""
    B = []
    prefix = []

    def rec(i, remaining):
        if i == NVAR:
            B.append(
                tuple(prefix)
            )
            return

        for z in range(remaining + 1):
            prefix.append(z)
            rec(
                i + 1,
                remaining - z
            )
            prefix.pop()

    rec(
        0,
        int(M)
    )

    return sorted(
        B,
        key=monomial_key
    )


def ell_vectors(b, leading):
    """Enumerate ell with sum_j ell_j*leading_j <= b using DFS pruning."""
    current = [0]*len(leading)
    used = [0]*NVAR

    def rec(j):
        if j == len(leading):
            yield tuple(current)
            return

        alpha = leading[j]
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

            for i in range(NVAR):
                used[i] += ell*alpha[i]

            yield from rec(j + 1)

            for i in range(NVAR):
                used[i] -= ell*alpha[i]

        current[j] = 0

    yield from rec(0)


def pick(b, S, polys, leading, powers):
    """Choose the admissible shift with largest weighted p-power."""
    best = None

    for ell in ell_vectors(
        b,
        leading
    ):
        used = tuple(
            sum(
                ell[j]*leading[j][i]
                for j in range(len(leading))
            )
            for i in range(NVAR)
        )
        a = tuple(
            b[i] - used[i]
            for i in range(NVAR)
        )
        g = mon(a)

        for z, f in zip(ell, polys):
            if z:
                g *= f**z

        if (
            supp(g).issubset(S)
            and lead(g) == b
        ):
            r = sum(
                powers[j]*ell[j]
                for j in range(len(powers))
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


# ============================================================
# Lattice construction and LLL polynomials
# ============================================================

def wt(e, X):
    """Equal-bound monomial scaling weight."""
    return ZZ(X)**sum(e)


def make_lat(B, G, p, t, X):
    """Build the scaled coefficient lattice."""
    pos = {
        b: i
        for i, b in enumerate(B)
    }
    rows = []

    for b in B:
        r, g0 = G[b]

        if t < r:
            raise RuntimeError(
                "target p-power is smaller than a shift p-power"
            )

        g = p**(t - r)*g0
        row = [ZZ(0)]*len(B)

        for e, a in g.dict().items():
            e = tuple(map(int, e))

            if e not in pos:
                raise RuntimeError(
                    "support outside monomial set"
                )

            row[pos[e]] += ZZ(a)*wt(
                e,
                X
            )

        rows.append(row)

    return Matrix(
        ZZ,
        rows
    )


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


# ============================================================
# Assumption 1 verification
# ============================================================

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
    """
    Check whether H1,...,H_NVAR generate a zero-dimensional ideal over QQ.
    """
    if not VERIFY_ASSUMPTION1:
        return None

    if len(H) < NVAR:
        return False

    Q = PolynomialRing(
        QQ,
        names=VAR_NAMES,
        order="degrevlex"
    )
    I = Q.ideal(
        [
            transport(f, Q)
            for f in H[:NVAR]
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
        NVAR
    )


def zero_polys(H, u):
    """
    Require H1,...,H_NVAR to vanish; return every vanishing H_i.

    This is the same experimental oracle-assisted selection used in the
    preceding CIHNP/MIHNP/ECHNP scripts: the true root is used only to
    identify which LLL rows are genuine integer-zero polynomials.
    """
    if len(H) < NVAR:
        return False, []

    is_zero = [
        ZZ(ev(f, u)) == 0
        for f in H
    ]
    ok = all(
        is_zero[:NVAR]
    )
    Z = [
        f
        for f, flag in zip(H, is_zero)
        if flag
    ]

    return ok, Z if ok else []


# ============================================================
# Modular Groebner solving and CRT lifting
# ============================================================

def lift(a, q):
    """Centered CRT lift."""
    modulus = prod(q)
    z = ZZ(crt(a, q))
    return z - modulus if z > modulus//2 else z


def solve(H, X, polys, powers, p, known):
    """
    Solve the selected auxiliary system modulo small primes and lift by CRT.

    The lattice ring uses weighted degree order, but each finite-field solve
    is performed in a fresh lexicographic ring x1>...>x_n.  A prime is used
    only when the full vanishing system has exactly one GF(q)-rational root.
    """
    if len(H) < NVAR:
        return None, 0.0

    residues = [
        []
        for _ in range(NVAR)
    ]  # x1,...,x_n.
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
        F = PolynomialRing(
            GF(q),
            names=VAR_NAMES,
            order="lex"
        )
        Pq = []

        for f in H:
            fq = transport(f, F)

            if fq:
                Pq.append(fq)

        if len(Pq) < NVAR:
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

        for i in range(NVAR):
            residues[i].append(a[i])

        modulus_product *= q

    u = tuple(
        lift(
            residues[i],
            moduli
        )
        for i in range(NVAR)
    )

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
        for f, e in zip(polys, powers)
    ):
        return None, gb_time

    if not valid_lcg_root(
        u,
        X,
        known,
        p
    ):
        return None, gb_time

    return u, gb_time


# ============================================================
# One experiment and repeated-run driver
# ============================================================

def run_once(run_id):
    """Generate and run one fresh unknown-multiplier LCG experiment."""
    p = rnd_prime(PBITS)
    _, _, seq = sample(p)
    X = ZZ(2)**UBITS
    u = tuple(
        ZZ(s) % X
        for s in seq
    )
    known = tuple(
        ZZ(s) - u[i]
        for i, s in enumerate(seq)
    )
    names, polys, leading, powers = system(
        known,
        p
    )

    if not valid_lcg_root(
        u,
        X,
        known,
        p
    ):
        raise RuntimeError(
            "invalid generated LCG root"
        )

    for name, f, alpha, e in zip(
        names,
        polys,
        leading,
        powers
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
    expected_dimension = binomial(
        NVAR + M,
        M
    )

    if len(B) != expected_dimension:
        raise RuntimeError(
            "unexpected simplex lattice dimension"
        )

    S = set(B)
    G = {
        b: pick(
            b,
            S,
            polys,
            leading,
            powers
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
        LR = L.LLL(
            delta=DELTA
        )

    except TypeError:
        LR = L.LLL()

    lll_time = time.perf_counter() - tic
    H = lll_polys(
        LR,
        B,
        X
    )

    if VERIFY_ASSUMPTION1:
        tic = time.perf_counter()
        a1 = assumption1(H)
        a1_time = time.perf_counter() - tic

    else:
        a1 = None
        a1_time = None

    ok, Z = zero_polys(
        H,
        u
    )
    rec = None
    gb_time = 0.0

    if ok:
        rec, gb_time = solve(
            Z,
            X,
            polys,
            powers,
            p,
            known
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

    print(
        " | ".join(output)
    )
    return bool(success)


def main():
    """Run N_RUNS fresh experiments."""
    if NSAMPLES < 3 or NVAR != NSAMPLES:
        raise ValueError(
            "require NSAMPLES=NVAR>=3"
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

    if len(LCG_WEIGHTS) != NVAR or any(
        ZZ(w) != w or w <= 0
        for w in LCG_WEIGHTS
    ):
        raise ValueError(
            "require positive integral LCG term-order weights"
        )

    if any(
        2*LCG_WEIGHTS[i + 1]
        <= LCG_WEIGHTS[i] + LCG_WEIGHTS[i + 2]
        for i in range(NSAMPLES - 2)
    ):
        raise ValueError(
            "LCG weights must be strictly concave"
        )

    if Q0 < 2 or QBAD < 1 or QMAX < 1:
        raise ValueError(
            "invalid modular-solving limits"
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
