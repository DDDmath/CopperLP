# LIPH unknown-degree POKE experiment for SageMath 10.8+.
# Small-root solving + Assumption 1 verification.
from sage.all import *
import time


N_RUNS = 1         # Number of independent experiments.
NVAR = 2           # Unknowns are x1 and x2.
PBITS = 1024       # Prime modulus bits.
X1_BITS = 149      # Bound X1=2^X1_BITS for x1.
X2_BITS = 360      # Bound X2=2^X2_BITS for x2.
M = 4              # LIPH polytope scale.
TH = QQ(301)/QQ(500)  # Rational approximation to theta ~= 0.60199.
DELTA = 0.99       # LLL parameter.

Q0 = 101           # First Groebner prime.
QBAD = 100         # Failed-prime limit.
QMAX = 500         # Tested-prime limit.
SEED = None        # Integer for a repeatable run sequence.

# Instance controls.
# Synthetic mode samples a fresh algebraic LIPH/POKE congruence each run.
# Manual mode expects
#     MANUAL_INSTANCE = (p, a, b, c, x1_root, x2_root)
# where the two roots are supplied only for experiment validation.
USE_SYNTHETIC_INSTANCE = True
MANUAL_INSTANCE = None

# Assumption 1 controls.
# Exact verification is usually inexpensive in two variables, but it can
# still be disabled.  msolve must be installed separately for the same
# SageMath installation before USE_MSOLVE_ASSUMPTION1 can be enabled.
VERIFY_ASSUMPTION1 = True
USE_MSOLVE_ASSUMPTION1 = False


# Algebraic relation:
#
#     f(x1,x2) = x1^2*x2 + a*x1*x2 + b*x2 + c == 0 mod p.
#
# Use degree-lexicographic order with x1 > x2 (equivalently x2 < x1).
# The unique degree-three term is therefore the required leading monomial
#
#     LM(f) = x1^2*x2,
#
# with canonical exponent alpha=(2,1).
VAR_NAMES = ("x1", "x2")
R = PolynomialRing(
    ZZ,
    names=VAR_NAMES,
    order="deglex"
)
x1, x2 = R.gens()
VARS = R.gens()
LEADING_EXPONENT = (2, 1)


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


def mon(e):
    """Build x1^e1*x2^e2 from a canonical exponent tuple."""
    if len(e) != NVAR:
        raise ValueError(
            "wrong monomial exponent length"
        )

    g = R(1)

    for x, z in zip(VARS, e):
        g *= x**int(z)

    return g


def ev(f, u):
    """Evaluate f at u=(u1,u2)."""
    if len(u) != NVAR:
        raise ValueError(
            "wrong root tuple length"
        )

    return f(*tuple(u))


def center_poly(f, modulus):
    """Center every coefficient of f modulo modulus."""
    g = R(0)

    for e, a in f.dict().items():
        a = cmod(a, modulus)

        if a:
            g += a*mon(
                tuple(map(int, e))
            )

    return g


def supp(f):
    """Return the exponent support in canonical (e1,e2) order."""
    return {
        tuple(map(int, e))
        for e, a in f.dict().items()
        if a
    }


def lead(f):
    """Return the leading exponent in canonical (e1,e2) order."""
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
    Copy f into another two-variable polynomial ring S.

    The lattice ring uses deglex, Assumption 1 uses degrevlex over QQ,
    and modular solving uses lex over GF(q).  Explicit coefficient and
    monomial transport avoids relying on an implicit coercion between rings
    with different term orders.
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
# LIPH / unknown-degree POKE algebraic instance
# ============================================================

def system(a, b, c, p):
    """Build the single public LIPH polynomial and its leading data."""
    f = center_poly(
        x1**2*x2
        + ZZ(a)*x1*x2
        + ZZ(b)*x2
        + ZZ(c),
        p
    )

    return (
        ("f",),
        (f,),
        (LEADING_EXPONENT,),
        (1,),
    )


def sample(p, bounds):
    """
    Generate a nondegenerate synthetic algebraic LIPH/POKE instance.

    Choose 0 < r1 < X1 and 0 < r2 < X2, sample public coefficients a,b,
    and define

        c = -r2*(r1^2 + a*r1 + b) mod p.

    This models the algebraic congruence layer after the LIPH/POKE reduction;
    it does not simulate an entire protocol transcript.
    """
    X1, X2 = map(ZZ, bounds)
    F = GF(p)

    for _ in range(1000):
        r1 = ZZ.random_element(
            1,
            X1
        )
        r2 = ZZ.random_element(
            1,
            X2
        )
        a = cmod(
            ZZ(F.random_element()),
            p
        )
        b = cmod(
            ZZ(F.random_element()),
            p
        )

        # Excluding zero coefficients and a zero denominator only removes
        # degenerate synthetic tests; the public equation itself is unchanged.
        if not a or not b:
            continue

        denominator = (
            r1**2
            + a*r1
            + b
        ) % p

        if not denominator:
            continue

        c = cmod(
            -r2*denominator,
            p
        )

        if not c:
            continue

        names, polys, leading, powers = system(
            a,
            b,
            c,
            p
        )
        f = polys[0]
        u = (
            ZZ(r1),
            ZZ(r2),
        )

        if lead(f) != LEADING_EXPONENT:
            raise RuntimeError(
                "unexpected leading monomial in generated instance"
            )

        if ZZ(ev(f, u)) % p:
            raise RuntimeError(
                "invalid generated LIPH congruence"
            )

        return (
            ZZ(a),
            ZZ(b),
            ZZ(c),
            u,
            names,
            polys,
            leading,
            powers,
        )

    raise RuntimeError(
        "could not generate a nondegenerate LIPH instance"
    )


def load_manual_instance(bounds):
    """Load and validate MANUAL_INSTANCE."""
    if MANUAL_INSTANCE is None:
        raise ValueError(
            "manual mode requires MANUAL_INSTANCE="
            "(p,a,b,c,x1_root,x2_root)"
        )

    if len(MANUAL_INSTANCE) != 6:
        raise ValueError(
            "MANUAL_INSTANCE must contain six integers"
        )

    p, a, b, c, r1, r2 = map(
        ZZ,
        MANUAL_INSTANCE
    )

    if p <= 2 or not p.is_prime(proof=False):
        raise ValueError(
            "manual modulus must be an odd prime"
        )

    if any(
        ZZ(Xi) >= p
        for Xi in bounds
    ):
        raise ValueError(
            "manual root bounds must be smaller than the modulus"
        )

    a = cmod(a, p)
    b = cmod(b, p)
    c = cmod(c, p)
    u = (
        ZZ(r1),
        ZZ(r2),
    )
    names, polys, leading, powers = system(
        a,
        b,
        c,
        p
    )

    if not valid_liph_root(
        u,
        bounds,
        a,
        b,
        c,
        p
    ):
        raise ValueError(
            "manual root is outside the bounds or does not satisfy f=0 mod p"
        )

    return (
        p,
        a,
        b,
        c,
        u,
        names,
        polys,
        leading,
        powers,
    )


def valid_liph_root(u, bounds, a, b, c, p):
    """Validate a candidate using the public equation and both root bounds."""
    if len(u) != NVAR or len(bounds) != NVAR:
        return False

    if any(
        z < 0 or z >= bounds[i]
        for i, z in enumerate(u)
    ):
        return False

    r1, r2 = map(ZZ, u)
    value = (
        r1**2*r2
        + ZZ(a)*r1*r2
        + ZZ(b)*r2
        + ZZ(c)
    )

    return value % p == 0


# ============================================================
# LIPH polytope and admissible shifts
# ============================================================

def inside(b):
    """
    Test membership in M*P(TH), where

        P(theta) = {(u,v)>=0 : v<=1, u-2v<=theta}.
    """
    if len(b) != NVAR:
        return False

    b1, b2 = map(ZZ, b)

    return (
        b1 >= 0
        and b2 >= 0
        and b2 <= M
        and QQ(b1 - 2*b2) <= M*TH
    )


def expected_dimension():
    """Closed-form lattice dimension for the scaled LIPH polytope."""
    return ZZ(M + 1)*(
        ZZ(M)
        + ZZ(floor(M*TH))
        + 1
    )


def monomial_key(b):
    """Ascending key compatible with degree-lex order x1>x2."""
    b1, b2 = b
    return (
        b1 + b2,
        b1,
        b2,
    )


def mons():
    """Enumerate every lattice monomial exponent in M*P(TH)."""
    B = []
    max_b1 = int(
        floor(
            M*(QQ(2) + TH)
        )
    )

    for b1 in range(max_b1 + 1):
        for b2 in range(int(M) + 1):
            b = (
                int(b1),
                int(b2),
            )

            if inside(b):
                B.append(b)

    return sorted(
        B,
        key=monomial_key
    )


def pick(b, S, f):
    """
    Choose the admissible shift with the largest p-adic power.

    For alpha=(2,1), candidates are

        x^(b-ell*alpha) * f^ell,

    with ell>=0 and ell*alpha<=b.  Since f=0 mod p, the shift has
    p-adic weight r=ell.
    """
    b1, b2 = b
    max_ell = min(
        b1//LEADING_EXPONENT[0],
        b2//LEADING_EXPONENT[1]
    )
    best = None

    for ell in range(max_ell + 1):
        a = (
            b1 - ell*LEADING_EXPONENT[0],
            b2 - ell*LEADING_EXPONENT[1],
        )
        g = mon(a)*f**ell

        if (
            supp(g).issubset(S)
            and lead(g) == b
        ):
            key = (
                ell,
                ell,
            )

            if best is None or key > best[0]:
                best = (
                    key,
                    g,
                )

    if best is None:
        raise RuntimeError(
            "no admissible shift for monomial %s"
            % (b,)
        )

    return best[0][0], best[1]


# ============================================================
# Lattice construction and LLL polynomials
# ============================================================

def wt(e, bounds):
    """Heterogeneous monomial scaling X1^e1*X2^e2."""
    if len(e) != NVAR or len(bounds) != NVAR:
        raise ValueError(
            "wrong scaling-vector length"
        )

    return prod(
        ZZ(Xi)**int(z)
        for Xi, z in zip(bounds, e)
    )


def make_lat(B, G, p, t, bounds):
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
                bounds
            )

        rows.append(row)

    return Matrix(
        ZZ,
        rows
    )


def to_poly(v, B, bounds):
    """Recover one primitive polynomial from an LLL row."""
    f = R(0)

    for a, b in zip(v, B):
        q, r = ZZ(a).quo_rem(
            wt(b, bounds)
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


def lll_polys(L, B, bounds):
    """Keep the LLL row order H1,H2,... ."""
    return [
        to_poly(
            L.row(i),
            B,
            bounds
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
            for i, z in enumerate(e)
            if z
        ]

        if len(nonzero) == 1:
            pure_power_seen[
                nonzero[0]
            ] = True

    return all(pure_power_seen)


def assumption1(H):
    """Check whether H1,H2 generate a zero-dimensional ideal over QQ."""
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
    Require H1,H2 to vanish and return every vanishing H_i.

    As in the preceding CIHNP-style experiment scripts, the simulated true
    root is used only to identify LLL rows that vanish over ZZ.  This is an
    oracle-assisted experiment-validation step, not an oracle-free attack API.
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


def solve(H, bounds, polys, powers, p, a, b, c):
    """
    Solve the selected auxiliary system modulo small primes and lift by CRT.

    A finite-field prime is accepted only when the complete selected system
    has exactly one GF(q)-rational root.  Because every selected H_i vanishes
    at the true integer root, that unique modular root must be its reduction.
    """
    if len(H) < NVAR:
        return None, 0.0

    residues = [
        []
        for _ in range(NVAR)
    ]  # x1,x2 in lex-ring generator order.
    moduli = []
    modulus_product = ZZ(1)
    target_product = 2*max(
        ZZ(Xi)
        for Xi in bounds
    )
    q = ZZ(Q0 - 1)
    bad = 0
    tried = 0
    gb_time = 0.0

    while modulus_product <= target_product:
        if bad >= QBAD or tried >= QMAX:
            return None, gb_time

        q = ZZ(
            next_prime(q)
        )
        tried += 1
        F = PolynomialRing(
            GF(q),
            names=VAR_NAMES,
            order="lex"
        )
        Pq = []

        for f in H:
            fq = transport(
                f,
                F
            )

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
            root_mod_q = [
                ZZ(V[0][v])
                for v in F.gens()
            ]

        except KeyError:
            bad += 1
            continue

        moduli.append(q)

        for i in range(NVAR):
            residues[i].append(
                root_mod_q[i]
            )

        modulus_product *= q

    u = tuple(
        lift(
            residues[i],
            moduli
        )
        for i in range(NVAR)
    )

    if any(
        z < 0 or z >= bounds[i]
        for i, z in enumerate(u)
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

    if not valid_liph_root(
        u,
        bounds,
        a,
        b,
        c,
        p
    ):
        return None, gb_time

    return u, gb_time


# ============================================================
# One experiment and repeated-run driver
# ============================================================

def run_once(run_id):
    """Generate/load and run one LIPH small-root experiment."""
    bounds = (
        ZZ(2)**X1_BITS,
        ZZ(2)**X2_BITS,
    )

    if USE_SYNTHETIC_INSTANCE:
        p = rnd_prime(PBITS)
        (
            a,
            b,
            c,
            u,
            names,
            polys,
            leading,
            powers,
        ) = sample(
            p,
            bounds
        )

    else:
        (
            p,
            a,
            b,
            c,
            u,
            names,
            polys,
            leading,
            powers,
        ) = load_manual_instance(
            bounds
        )

    if not valid_liph_root(
        u,
        bounds,
        a,
        b,
        c,
        p
    ):
        raise RuntimeError(
            "invalid generated or loaded LIPH root"
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

    if len(B) != expected_dimension():
        raise RuntimeError(
            "unexpected LIPH lattice dimension"
        )

    if len(B) != len(set(B)) or any(
        not inside(z)
        for z in B
    ):
        raise RuntimeError(
            "invalid LIPH monomial enumeration"
        )

    S = set(B)
    f = polys[0]
    G = {
        beta: pick(
            beta,
            S,
            f
        )
        for beta in B
    }
    t_power = max(
        r
        for r, _ in G.values()
    )

    # For this polytope, beta=(2M,M) admits f^M and every beta has b2<=M,
    # so the generalized target power must be exactly t_M=M.
    if t_power != M:
        raise RuntimeError(
            "unexpected LIPH target p-power"
        )

    L = make_lat(
        B,
        G,
        p,
        t_power,
        bounds
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
        bounds
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
            bounds,
            polys,
            powers,
            p,
            a,
            b,
            c
        )

    success = rec == u
    output = [
        "modulus bits=%d" % p.nbits(),
        "unknown bits=(%d,%d)" % (
            X1_BITS,
            X2_BITS
        ),
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
    """Run N_RUNS experiments."""
    if ZZ(PBITS) != PBITS or PBITS < 3:
        raise ValueError(
            "require PBITS to be an integer at least 3"
        )

    if NVAR != 2:
        raise ValueError(
            "this script requires NVAR=2"
        )

    if ZZ(X1_BITS) != X1_BITS or not 0 < X1_BITS < PBITS:
        raise ValueError(
            "require integral 0 < X1_BITS < PBITS"
        )

    if ZZ(X2_BITS) != X2_BITS or not 0 < X2_BITS < PBITS:
        raise ValueError(
            "require integral 0 < X2_BITS < PBITS"
        )

    if ZZ(M) != M or M < 1:
        raise ValueError(
            "require M to be a positive integer"
        )

    if QQ(TH) <= 0:
        raise ValueError(
            "require TH to be positive"
        )

    if not QQ(1)/QQ(4) < QQ(DELTA) < 1:
        raise ValueError(
            "require 1/4 < DELTA < 1"
        )

    if ZZ(N_RUNS) != N_RUNS or N_RUNS < 1:
        raise ValueError(
            "require N_RUNS to be a positive integer"
        )

    if any(
        ZZ(z) != z
        for z in (Q0, QBAD, QMAX)
    ) or Q0 < 2 or QBAD < 1 or QMAX < 1:
        raise ValueError(
            "invalid modular-solving limits"
        )

    if not USE_SYNTHETIC_INSTANCE and MANUAL_INSTANCE is None:
        raise ValueError(
            "manual mode requires MANUAL_INSTANCE"
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
