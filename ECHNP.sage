# ECHNP 3-sample experiment for SageMath 10.8+.
# Small-root solving + Assumption 1 verification.
from sage.all import *
import itertools
import time

# ============================================================
# User configuration
# ============================================================
#
# Table 2 parameter presets: (PBITS, UBITS, M)
#
#   (512, 96, 1)
#   (256, 69, 2)
#
# Both experiments use:
#
#   NSAMPLES = 3
#   NVAR = 4
#
# The TH parameters remain unchanged. To reproduce a
# particular row, manually replace PBITS, UBITS, and M below.
#

N_RUNS = 1        # Number of independent experiments.
NSAMPLES = 3      # Number of ECHNP oracle sample pairs.
NVAR = 4          # x1 for x(P), then x2,x3,x4 for the three sample sums.
PBITS = 256       # Modulus bits.
UBITS = 69        # Unknown low bits of each individual x-coordinate.
M = 2             # Polytope scale.
TH = (
    QQ(140)/QQ(100),
    QQ(232)/QQ(100),
    QQ(236)/QQ(100),
)                  # ECHNP polytope parameters (theta1,theta2,theta3).
DELTA = 0.99       # LLL parameter.

Q0 = 101           # First Groebner prime.
QBAD = 100         # Failed-prime limit.
QMAX = 500         # Tested-prime limit.
SEED = None        # Integer for a repeatable run sequence.

# Assumption 1 controls.
# Exact verification can be expensive.  Testing only small instances is
# recommended; msolve can be faster but must be installed separately.
VERIFY_ASSUMPTION1 = True
USE_MSOLVE_ASSUMPTION1 = False


# Ring order: x4 > x3 > x2 > x1.  All exponent tuples used below are
# nevertheless stored in canonical order (e1,e2,e3,e4).
#
# With this order the required leading monomials are
#     LM(f2)  = x1^2*x2,  LM(f3)  = x1^2*x3,
#     LM(f4)  = x1^2*x4,  LM(f23) = x2*x3,
#     LM(f24) = x2*x4,    LM(f34) = x3*x4.
R = PolynomialRing(
    ZZ,
    names=tuple(
        "x%d" % i
        for i in range(NVAR, 0, -1)
    ),
    order="lex"
)

RGENS = R.gens()                 # (x4,x3,x2,x1)
VARS = tuple(reversed(RGENS))    # (x1,x2,x3,x4)


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


def ex(e):
    """Convert ring exponents (e4,e3,e2,e1) to (e1,e2,e3,e4)."""
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
    """Evaluate at canonical root tuple u=(u1,u2,u3,u4)."""
    return f(
        *tuple(reversed(tuple(u)))
    )


def center_poly(f, modulus):
    """Center all polynomial coefficients modulo modulus."""
    g = R(0)

    for ering, a in f.dict().items():
        a = cmod(a, modulus)

        if a:
            g += a*mon(ex(ering))

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


# ============================================================
# Elliptic-curve ECHNP instance generation
# ============================================================

def sample(p):
    """
    Generate one fresh three-sample ECHNP instance over GF(p).

    The curve is y^2 = x^3 + curve_a*x + curve_b.  The secret point is
    P, the public base point is R0, and Q_i=[t_i]R0.  Each Q_i is chosen
    with a distinct x-coordinate and nonzero y-coordinate.  These
    conditions also make the pair-relation linear system nonsingular.
    """
    F = GF(p)

    while True:
        curve_a = F.random_element()
        curve_b = F.random_element()

        if 4*curve_a**3 + 27*curve_b**2 == 0:
            continue

        E = EllipticCurve(
            F,
            [curve_a, curve_b]
        )
        O = E(0)
        secret = E.random_point()
        base = E.random_point()

        if (
            secret == O
            or base == O
            or base[1] == 0
        ):
            continue

        x_secret = ZZ(secret[0])
        t_values = []
        points = []

        for _ in range(500):
            if len(points) == NSAMPLES:
                break

            ti = ZZ.random_element(1, p)

            if ti in t_values:
                continue

            Q = ti*base

            if Q == O or Q[1] == 0:
                continue

            plus = secret + Q
            minus = secret - Q

            if plus == O or minus == O:
                continue

            qx = ZZ(Q[0])

            if qx == x_secret:
                continue

            if any(
                ZZ(Q0[0]) == qx
                for Q0 in points
            ):
                continue

            t_values.append(ti)
            points.append(Q)

        if len(points) == NSAMPLES:
            return (
                E,
                ZZ(curve_a),
                ZZ(curve_b),
                secret,
                base,
                tuple(t_values),
                tuple(points),
            )


def leak_sample(secret, Q, X):
    """
    Combine the two oracle leaks for Q.

    If
        S_plus  = x(P+Q) = known_plus  + error_plus,
        S_minus = x(P-Q) = known_minus + error_minus,
    then the polynomial variable is
        error_plus + error_minus < 2*X,
    while the known constant is known_plus + known_minus.
    """
    s_plus = ZZ((secret + Q)[0])
    s_minus = ZZ((secret - Q)[0])
    known_plus, error_plus = split_low(
        s_plus,
        X
    )
    known_minus, error_minus = split_low(
        s_minus,
        X
    )
    known_sum = known_plus + known_minus
    error_sum = error_plus + error_minus

    if not 0 <= error_sum < 2*X:
        raise RuntimeError(
            "invalid sample-sum bound"
        )

    return (
        ZZ(Q[0]),
        ZZ(known_sum),
        ZZ(error_sum),
    )


# ============================================================
# ECHNP polynomial system
# ============================================================

def original_poly(
    qx,
    known_x,
    known_sum,
    curve_a,
    curve_b,
    p,
    sample_index
):
    """
    Construct the original sample polynomial modulo p.

    With X=known_x+x1 and Y=known_sum+x_{sample_index+2},

        Y*(X-q)^2
        - 2*(q*X^2 + (curve_a+q^2)*X + curve_a*q + 2*curve_b)
        == 0 mod p.
    """
    q = ZZ(qx)
    X = known_x + VARS[0]
    Y = known_sum + VARS[sample_index + 1]
    f = (
        Y*(X - q)**2
        - 2*(
            q*X**2
            + (curve_a + q**2)*X
            + curve_a*q
            + 2*curve_b
        )
    )

    return center_poly(
        f,
        p
    )


def relation_coefficients(qi, qj, curve_a, curve_b, p):
    """
    Compute the six coefficients of the pair relation over GF(p).

    For the full-coordinate rational functions Y_i(X),Y_j(X), solve

        Y_i*Y_j
        + (r1*X+r0)*Y_j
        + (s1*X+s0)*Y_i
        + (t1*X+t0) = 0.
    """
    F = GF(p)
    U = PolynomialRing(
        F,
        names=("T",)
    )
    T = U.gen()
    qi = F(qi)
    qj = F(qj)
    curve_a = F(curve_a)
    curve_b = F(curve_b)

    Ni = 2*(
        qi*T**2
        + (curve_a + qi**2)*T
        + curve_a*qi
        + 2*curve_b
    )
    Nj = 2*(
        qj*T**2
        + (curve_a + qj**2)*T
        + curve_a*qj
        + 2*curve_b
    )
    Di = (T - qi)**2
    Dj = (T - qj)**2

    terms = [
        T*Nj*Di,       # r1
        Nj*Di,         # r0
        T*Ni*Dj,       # s1
        Ni*Dj,         # s0
        T*Di*Dj,       # t1
        Di*Dj,         # t0
    ]
    target = -Ni*Nj
    mat = Matrix(
        F,
        [
            [
                term[d]
                for term in terms
            ]
            for d in range(6)
        ]
    )
    vec = vector(
        F,
        [
            target[d]
            for d in range(6)
        ]
    )

    try:
        sol = mat.solve_right(vec)

    except Exception as exc:
        raise ValueError(
            "singular ECHNP pair-relation system"
        ) from exc

    r1, r0, s1, s0, t1, t0 = sol
    identity = (
        Ni*Nj
        + (r1*T + r0)*Nj*Di
        + (s1*T + s0)*Ni*Dj
        + (t1*T + t0)*Di*Dj
    )

    if identity:
        raise RuntimeError(
            "invalid ECHNP pair relation"
        )

    return tuple(
        ZZ(z)
        for z in sol
    )


def pair_poly(
    i,
    j,
    known_x,
    known_sums,
    qx,
    curve_a,
    curve_b,
    p
):
    """Construct f_{i+2,j+2}=0 mod p in shifted variables."""
    if not 0 <= i < j < NSAMPLES:
        raise ValueError(
            "invalid sample-pair indices"
        )

    r1, r0, s1, s0, t1, t0 = relation_coefficients(
        qx[i],
        qx[j],
        curve_a,
        curve_b,
        p
    )
    X = known_x + VARS[0]
    Yi = known_sums[i] + VARS[i + 1]
    Yj = known_sums[j] + VARS[j + 1]
    f = (
        Yi*Yj
        + (r1*X + r0)*Yj
        + (s1*X + s0)*Yi
        + (t1*X + t0)
    )

    return center_poly(
        f,
        p
    )


def original_lead(i):
    """Leading exponent of the i-th original sample polynomial."""
    e = [0]*NVAR
    e[0] = 2
    e[i + 1] = 1
    return tuple(e)


def pair_lead(i, j):
    """Leading exponent of the pair polynomial on samples i,j."""
    e = [0]*NVAR
    e[i + 1] = 1
    e[j + 1] = 1
    return tuple(e)


def system(
    known_x,
    known_sums,
    qx,
    curve_a,
    curve_b,
    p
):
    """Build f2,f3,f4,f23,f24,f34 and their leading data."""
    names = []
    polys = []
    leading = []
    powers = []

    for i in range(NSAMPLES):
        f = original_poly(
            qx[i],
            known_x,
            known_sums[i],
            curve_a,
            curve_b,
            p,
            i
        )
        names.append(
            "f%d" % (i + 2)
        )
        polys.append(f)
        leading.append(
            original_lead(i)
        )
        powers.append(1)

    for i, j in itertools.combinations(
        range(NSAMPLES),
        2
    ):
        f = pair_poly(
            i,
            j,
            known_x,
            known_sums,
            qx,
            curve_a,
            curve_b,
            p
        )
        names.append(
            "f%d%d" % (i + 2, j + 2)
        )
        polys.append(f)
        leading.append(
            pair_lead(i, j)
        )
        powers.append(1)

    return names, polys, leading, powers


# ============================================================
# ECHNP polytope and admissible shifts
# ============================================================

def inside(b):
    """Test membership in the scaled ECHNP polytope."""
    b1, b2, b3, b4 = b
    theta1, theta2, theta3 = TH

    return (
        min(b) >= 0
        and QQ(b1) <= M*theta1
        and b2 <= M
        and b3 <= M
        and b4 <= M
        and QQ(b1 + b2 + b3) <= M*theta2
        and QQ(b1 + b2 + b4) <= M*theta2
        and QQ(b1 + b3 + b4) <= M*theta2
        and QQ(b1 + b2 + b3 + b4) <= M*theta3
    )


def mons():
    """Enumerate lattice monomials in ascending ring-monomial order."""
    theta1, _, _ = TH
    B = []

    for b1 in range(
        int(floor(M*theta1)) + 1
    ):
        for tail in itertools.product(
            range(M + 1),
            repeat=NSAMPLES
        ):
            b = (
                int(b1),
                *tuple(map(int, tail))
            )

            if inside(b):
                B.append(b)

    return sorted(
        B,
        key=lambda b: tuple(reversed(b))
    )


def ell_vectors(b, leading):
    """Enumerate ell with sum_j ell_j*leading_j <= b by DFS pruning."""
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

def wt(e, bounds):
    """Heterogeneous monomial scaling weight."""
    return prod(
        ZZ(Xi)**int(a)
        for Xi, a in zip(bounds, e)
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
        g = p**(t - r)*g0
        row = [ZZ(0)]*len(B)

        for e, a in g.dict().items():
            e = ex(e)

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
            for i, a in enumerate(e)
            if a
        ]

        if len(nonzero) == 1:
            pure_power_seen[
                nonzero[0]
            ] = True

    return all(pure_power_seen)


def assumption1(H):
    """Check whether the first four auxiliary polynomials are zero-dimensional."""
    if not VERIFY_ASSUMPTION1:
        return None

    if len(H) < NVAR:
        return False

    Q = PolynomialRing(
        QQ,
        names=R.variable_names(),
        order="degrevlex"
    )
    I = Q.ideal(
        [
            Q(f)
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
    """Require H1,...,H4 to vanish; return every vanishing H_i."""
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


def valid_echnp_root(
    u,
    bounds,
    known_x,
    known_sums,
    qx,
    curve_a,
    curve_b,
    p
):
    """Validate a candidate using only the public ECHNP equations and bounds."""
    if len(u) != NVAR:
        return False

    if any(
        z < 0 or z >= bounds[i]
        for i, z in enumerate(u)
    ):
        return False

    Xcoord = ZZ(known_x + u[0])

    if not 0 <= Xcoord < p:
        return False

    F = GF(p)

    if not F(
        Xcoord**3
        + curve_a*Xcoord
        + curve_b
    ).is_square():
        return False

    for i in range(NSAMPLES):
        Ysum = ZZ(
            known_sums[i]
            + u[i + 1]
        )
        q = ZZ(qx[i])

        if not 0 <= Ysum <= 2*(p - 1):
            return False

        if (Xcoord - q) % p == 0:
            return False

        value = (
            Ysum*(Xcoord - q)**2
            - 2*(
                q*Xcoord**2
                + (curve_a + q**2)*Xcoord
                + curve_a*q
                + 2*curve_b
            )
        )

        if value % p:
            return False

    return True


def solve(
    H,
    bounds,
    polys,
    powers,
    p,
    known_x,
    known_sums,
    qx,
    curve_a,
    curve_b
):
    """Solve the selected auxiliary system modulo small primes and CRT."""
    if len(H) < NVAR:
        return None, 0.0

    residues = [
        []
        for _ in range(NVAR)
    ]  # Ring order x4,x3,x2,x1.
    moduli = []
    modulus_product = ZZ(1)
    target_product = 2*max(
        ZZ(z)
        for z in bounds
    )
    q = ZZ(Q0 - 1)
    bad = 0
    tried = 0
    gb_time = 0.0

    while modulus_product <= target_product:
        if bad >= QBAD or tried >= QMAX:
            return None, gb_time

        q = ZZ(next_prime(q))
        tried += 1
        F = R.change_ring(
            GF(q)
        )
        Pq = []

        for f in H:
            fq = F(f)

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

    ring_root = [
        lift(
            residues[i],
            moduli
        )
        for i in range(NVAR)
    ]
    u = tuple(
        reversed(ring_root)
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

    if not valid_echnp_root(
        u,
        bounds,
        known_x,
        known_sums,
        qx,
        curve_a,
        curve_b,
        p
    ):
        return None, gb_time

    return u, gb_time


# ============================================================
# One experiment and repeated-run driver
# ============================================================

def run_once(run_id):
    """Generate and run one fresh ECHNP experiment."""
    p = rnd_prime(PBITS)
    (
        _,
        curve_a,
        curve_b,
        secret,
        _,
        _,
        sample_points,
    ) = sample(p)

    X = ZZ(2)**UBITS
    bounds = (
        X,
        2*X,
        2*X,
        2*X,
    )
    known_x, root_x = split_low(
        ZZ(secret[0]),
        X
    )
    qx = []
    known_sums = []
    root_sums = []

    for Q in sample_points:
        qi, known_i, root_i = leak_sample(
            secret,
            Q,
            X
        )
        qx.append(qi)
        known_sums.append(known_i)
        root_sums.append(root_i)

    qx = tuple(qx)
    known_sums = tuple(known_sums)
    u = (
        ZZ(root_x),
        *tuple(map(ZZ, root_sums))
    )
    names, polys, leading, powers = system(
        known_x,
        known_sums,
        qx,
        curve_a,
        curve_b,
        p
    )

    if not valid_echnp_root(
        u,
        bounds,
        known_x,
        known_sums,
        qx,
        curve_a,
        curve_b,
        p
    ):
        raise RuntimeError(
            "invalid generated ECHNP root"
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
            known_x,
            known_sums,
            qx,
            curve_a,
            curve_b
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
    if NSAMPLES != 3 or NVAR != 4:
        raise ValueError(
            "this script requires NSAMPLES=3 and NVAR=4"
        )

    if not 0 < UBITS < PBITS:
        raise ValueError(
            "require 0 < UBITS < PBITS"
        )

    if ZZ(M) != M or M < 1:
        raise ValueError(
            "require M to be a positive integer"
        )

    if len(TH) != 3 or any(
        QQ(z) <= 0
        for z in TH
    ):
        raise ValueError(
            "require three positive polytope parameters"
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
