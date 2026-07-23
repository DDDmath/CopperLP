# ============================================================
#  Algorithm 1: Computation of the asymptotic coefficients
#
#  Input:
#      polys = [f_1, ..., f_n]
#      P_ineqs = [((a_1, ..., a_k), b), ...]
#          representing the rational polytope
#          P = {x in R^k : <a_t, x> <= b_t, t in T}
#      leading_exponents = [alpha_1, ..., alpha_n]
#      e = [e_1, ..., e_n]
#
#  Output:
#      Delta, Y, Vert(Y), affine restrictions varphi_i,
#      full-dimensional cells Q_i, and the asymptotic coefficients
#          sigma_1, ..., sigma_k, sigma_0.
#
#      The asymptotic small-root condition is
#          prod_{j=1}^k X_j^(sigma_j) < M^(sigma_0 - epsilon).
#
#      If X_1 = ... = X_k = X, then
#          X < M^(theta - epsilon'),
#          theta = sigma_0 / (sigma_1 + ... + sigma_k).
# ============================================================

from sage.all import *

try:
    from sage.geometry.triangulation.point_configuration import PointConfiguration
except ImportError:
    PointConfiguration = None


# ============================================================
#  1. Basic utilities
# ============================================================

def decimal6(x):
    """
    Return a decimal string with six digits after the decimal point.
    """
    return "{:.6f}".format(float(x))


def matrix_and_rhs_from_face_inequalities(P_ineqs):
    """
    Construct the matrix A and vector b from the face inequalities
        <a_t, x> <= b_t, t in T.

    Output:
        A: |T| x k matrix with rows a_t^T.
        b: |T|-dimensional vector with entries b_t.
    """
    if len(P_ineqs) == 0:
        raise ValueError("P_ineqs must contain at least one face inequality.")

    A = matrix(QQ, [list(vector(QQ, a)) for a, _ in P_ineqs])
    b = vector(QQ, [QQ(bt) for _, bt in P_ineqs])
    return A, b


def polyhedron_from_leq(P_ineqs, ambient_dim=None):
    """
    Construct a Sage Polyhedron from inequalities of the form
        <a_t, x> <= b_t.

    Sage uses inequalities in the form c + A*x >= 0.  Therefore
        <a_t, x> <= b_t
    is encoded as
        b_t - <a_t, x> >= 0.
    """
    ieqs = []
    for a, b in P_ineqs:
        a = vector(QQ, a)
        ieqs.append([QQ(b)] + list(-a))
    return Polyhedron(ieqs=ieqs, ambient_dim=ambient_dim)


def check_polytope_P(P, k, require_origin=True, require_nonnegative=True):
    """
    Check the assumptions on P used in Algorithm 1:
        1. P is nonempty;
        2. P has ambient dimension k;
        3. P is full-dimensional;
        4. P is bounded;
        5. P is rational;
        6. optionally, P is contained in R_{>=0}^k;
        7. optionally, 0 is contained in P.
    """
    if P.is_empty():
        raise ValueError("P must be nonempty.")

    if P.ambient_dim() != k:
        raise ValueError(f"The ambient dimension of P must be {k}.")

    if P.dim() != k:
        raise ValueError("P must be full-dimensional.")

    if not P.is_compact():
        raise ValueError("P must be bounded; unbounded regions are not covered by Algorithm 1.")

    try:
        vertices = [vector(QQ, v) for v in P.vertices_list()]
    except Exception:
        raise ValueError("P must be a rational polytope.")

    if require_nonnegative:
        for v in vertices:
            if any(x < 0 for x in v):
                raise ValueError("P should be contained in R_{>=0}^k.")

    if require_origin:
        zero = vector(QQ, [0] * k)
        if not P.contains(zero):
            print("Warning: 0 is not contained in P.")
            print("For the usual monomial family mP cap L, one normally assumes 0 in P.")

    return True


# Backward-compatible name.
check_region_P = check_polytope_P


def support_exponents(f):
    """
    Return the exponent set A(f).
    """
    return [vector(QQ, exponent) for exponent in f.exponents()]


def leading_exponent(f, alpha=None):
    """
    Return the leading exponent alpha_j.

    If alpha is given, it is used directly. Otherwise, the leading monomial
    is computed from the current Sage monomial order.
    """
    if alpha is not None:
        return vector(QQ, alpha)
    return vector(QQ, f.lm().degrees())


def lattice_determinant(k, lattice_basis=None):
    """
    Compute det(L).

    If lattice_basis is None, then L = Z^k. Otherwise, lattice_basis is a
    k x k integer matrix whose rows are a basis of L.
    """
    if lattice_basis is None:
        B = identity_matrix(ZZ, k)
    else:
        B = matrix(ZZ, lattice_basis)
        if B.nrows() != k or B.ncols() != k:
            raise ValueError("lattice_basis must be a k x k integer matrix.")

    det_L = abs(B.det())
    if det_L == 0:
        raise ValueError("lattice_basis must be invertible.")
    return QQ(det_L), B


def format_affine(varphi, varnames=None):
    """
    Format an affine function
        varphi(x) = c0 + c^T x.
    """
    c0, c = varphi
    k = len(c)
    if varnames is None:
        varnames = [f"x{i+1}" for i in range(k)]

    terms = []
    if c0 != 0:
        terms.append(str(c0))
    for i in range(k):
        if c[i] != 0:
            terms.append(f"({c[i]})*{varnames[i]}")

    if not terms:
        return "0"
    return " + ".join(terms).replace("+ (-", "- (")


# ============================================================
#  2. Delta, dual feasible polyhedron Y, affine restrictions,
#     and cells Q_i
# ============================================================

def delta_matrix_from_face_inequalities(P_ineqs, support_sets, alphas):
    """
    Compute the matrix Delta = (Delta_t(f_j))_{t in T, 1 <= j <= n}, where
        Delta_t(f_j) = max_{u in A(f_j)} <a_t, u - alpha_j>.
    """
    A, _ = matrix_and_rhs_from_face_inequalities(P_ineqs)
    r = A.nrows()
    n = len(support_sets)
    Delta = matrix(QQ, r, n)

    for t in range(r):
        a_t = vector(QQ, A.row(t))
        for j in range(n):
            Delta[t, j] = max(a_t.dot_product(u - alphas[j]) for u in support_sets[j])

    return Delta


# Backward-compatible name.
delta_matrix_from_ineqs = delta_matrix_from_face_inequalities


def dual_feasible_polyhedron_Y(Delta, e):
    """
    Construct the dual feasible polyhedron
        Y = {z in R^{|T|}_{>=0} : Delta^T z >= e}.
    """
    r = Delta.nrows()
    n = Delta.ncols()
    e = vector(QQ, e)

    ieqs = []

    # z_t >= 0 for every t in T.
    for t in range(r):
        row = [QQ(0)] + [QQ(0)] * r
        row[1 + t] = QQ(1)
        ieqs.append(row)

    # Delta^T z >= e.
    for j in range(n):
        col = list(Delta.column(j))
        ieqs.append([-e[j]] + col)

    Y = Polyhedron(ieqs=ieqs)
    if Y.is_empty():
        raise ValueError("The dual feasible polyhedron Y is empty.")
    return Y


# Backward-compatible name.
dual_polyhedron_Y = dual_feasible_polyhedron_Y


def affine_restriction_from_dual_vertex(z, A, b):
    """
    Given a vertex z^(i) of Y, construct the affine restriction
        L_i(x) = (b - A*x)^T z^(i)
               = b^T z^(i) - (A^T z^(i))^T x.

    Output:
        (c0, c), representing L_i(x) = c0 + c^T x.
    """
    z = vector(QQ, z)
    c0 = b.dot_product(z)
    c = -(A.transpose() * z)
    return QQ(c0), vector(QQ, c)


def affine_from_dual_vertex(z, P_ineqs):
    """
    Backward-compatible wrapper for affine_restriction_from_dual_vertex.
    """
    A, b = matrix_and_rhs_from_face_inequalities(P_ineqs)
    return affine_restriction_from_dual_vertex(z, A, b)


def distinct_affine_restrictions_from_Y(Y, A, b):
    """
    Enumerate Vert(Y) = {z^(1), ..., z^(W)} and construct the distinct
    restrictions among
        L_1|P, ..., L_W|P.

    Since P is full-dimensional, equality on P is the same as equality as
    affine functions in the ambient space.
    """
    varphis = []
    seen = set()

    vertices_Y = [vector(QQ, v) for v in Y.vertices_list()]
    if len(vertices_Y) == 0:
        raise ValueError("Y has no vertices; Algorithm 1 cannot enumerate Vert(Y).")

    for z in vertices_Y:
        varphi = affine_restriction_from_dual_vertex(z, A, b)
        key = (QQ(varphi[0]), tuple(vector(QQ, varphi[1])))
        if key not in seen:
            seen.add(key)
            varphis.append(varphi)

    return varphis


def distinct_affines_from_Y(Y, P_ineqs):
    """
    Backward-compatible wrapper for distinct_affine_restrictions_from_Y.
    """
    A, b = matrix_and_rhs_from_face_inequalities(P_ineqs)
    return distinct_affine_restrictions_from_Y(Y, A, b)


def cell_polyhedron(P, i, varphis):
    """
    Construct the cell
        Q_i = P cap intersection_h {x in R^k : varphi_i(x) <= varphi_h(x)}.
    """
    c0_i, c_i = varphis[i]
    comparison_ieqs = []

    for h, (c0_h, c_h) in enumerate(varphis):
        if h == i:
            continue

        # c0_i + c_i^T*x <= c0_h + c_h^T*x
        # is equivalent to
        # (c_i - c_h)^T*x <= c0_h - c0_i.
        a = c_i - c_h
        b = c0_h - c0_i
        comparison_ieqs.append([QQ(b)] + list(-a))

    if comparison_ieqs:
        return P.intersection(Polyhedron(ieqs=comparison_ieqs))
    return P


# Backward-compatible name.
piece_polyhedron = cell_polyhedron


# ============================================================
#  3. Integration by triangulation and simplex integration
# ============================================================

def triangulate_polytope(Q):
    """
    Return a triangulation of Q.

    Each simplex is represented by a list of vertices [v_0, ..., v_k].
    """
    if Q.is_empty():
        return []

    if Q.dim() != Q.ambient_dim():
        return []

    vertices = [vector(QQ, v) for v in Q.vertices_list()]
    k = Q.ambient_dim()

    if len(vertices) == 0:
        return []
    if len(vertices) == k + 1:
        return [vertices]

    try:
        triangulation = Q.triangulate()
        simplices_idx = triangulation.simplices() if hasattr(triangulation, "simplices") else triangulation
        simplices = []
        for simplex_idx in simplices_idx:
            simplex_idx = list(simplex_idx)
            if len(simplex_idx) == k + 1:
                simplices.append([vertices[i] for i in simplex_idx])
        if len(simplices) > 0:
            return simplices
    except Exception:
        pass

    if PointConfiguration is None:
        raise RuntimeError("This Sage environment cannot triangulate the polytope.")

    point_configuration = PointConfiguration([tuple(v) for v in vertices])
    triangulation = point_configuration.triangulate()
    simplices_idx = triangulation.simplices() if hasattr(triangulation, "simplices") else triangulation

    simplices = []
    for simplex_idx in simplices_idx:
        simplex_idx = list(simplex_idx)
        if len(simplex_idx) == k + 1:
            simplices.append([vertices[i] for i in simplex_idx])

    return simplices


def simplex_volume(simplex):
    """
    Compute the Euclidean volume of a k-dimensional simplex.
    """
    k = len(simplex[0])
    v0 = simplex[0]
    M = matrix(QQ, [list(v - v0) for v in simplex[1:]])
    return abs(M.det()) / factorial(k)


def integrate_affine_on_simplex(simplex, varphi):
    """
    Integrate an affine function on a simplex:
        integral_S varphi = Vol(S) * average of varphi over the vertices.
    """
    c0, c = varphi
    vol = simplex_volume(simplex)
    avg = sum(c0 + c.dot_product(v) for v in simplex) / QQ(len(simplex))
    return vol * avg


def integrate_coordinate_on_simplex(simplex, j):
    """
    Integrate the coordinate function x_j on a simplex.
    Here j is zero-based.
    """
    vol = simplex_volume(simplex)
    avg = sum(v[j] for v in simplex) / QQ(len(simplex))
    return vol * avg


def integrate_affine_on_polytope(Q, varphi):
    """
    Integrate an affine function over a full-dimensional polytope Q.
    """
    return sum(integrate_affine_on_simplex(S, varphi) for S in triangulate_polytope(Q))


def integrate_coordinate_on_polytope(Q, j):
    """
    Integrate the coordinate function x_j over a full-dimensional polytope Q.
    Here j is zero-based.
    """
    return sum(integrate_coordinate_on_simplex(S, j) for S in triangulate_polytope(Q))


# ============================================================
#  4. Common bound when X_1 = ... = X_k
# ============================================================

def common_X_bound_exponent(sigma, sigma0):
    """
    If X_1 = ... = X_k = X, then Algorithm 1 gives
        X^(sigma_1 + ... + sigma_k) < M^(sigma_0 - epsilon).

    Therefore the common exponent is
        theta = sigma_0 / (sigma_1 + ... + sigma_k).
    """
    sigma_sum = sum(sigma)
    if sigma_sum == 0:
        raise ValueError("sum_j sigma_j is zero, so the common X exponent is undefined.")
    return QQ(sigma0) / QQ(sigma_sum)


# ============================================================
#  5. Main routine: Algorithm 1 with P given by face inequalities
# ============================================================

def compute_asymptotic_coefficients_from_face_inequalities(polys,
                                                           P_ineqs,
                                                           leading_exponents=None,
                                                           e=None,
                                                           lattice_basis=None,
                                                           require_origin=True,
                                                           require_nonnegative=True,
                                                           X_names=None,
                                                           M_name="M"):
    """
    Compute the asymptotic coefficients in Algorithm 1.

    Input:
        polys = [f_1, ..., f_n].
        P_ineqs = [((a_1, ..., a_k), b), ...], meaning <a_t, x> <= b_t.
        leading_exponents = [alpha_1, ..., alpha_n].
        e = [e_1, ..., e_n]; if omitted, e = (1, ..., 1)^T.
        lattice_basis is optional; if omitted, L = Z^k.

    Output:
        A dictionary containing A, b, Delta, Y, Vert(Y), varphi_i, Q_i,
        sigma_1, ..., sigma_k, sigma_0, and the final asymptotic condition.

    Note:
        Algorithm 1 is stated for L = Z^k.  The optional lattice_basis parameter
        keeps the standard det(L) rescaling, so all coefficients are divided by
        det(L) when a sublattice L is supplied.
    """
    if len(polys) == 0:
        raise ValueError("polys must contain at least one polynomial.")

    n = len(polys)
    k = polys[0].parent().ngens()

    if e is None:
        e = [1] * n
    e = vector(QQ, e)
    if len(e) != n:
        raise ValueError("The length of e must equal the number of polynomials.")

    P = polyhedron_from_leq(P_ineqs, ambient_dim=k)
    check_polytope_P(P, k,
                     require_origin=require_origin,
                     require_nonnegative=require_nonnegative)

    A, b = matrix_and_rhs_from_face_inequalities(P_ineqs)
    support_sets = [support_exponents(f) for f in polys]

    if leading_exponents is None:
        alphas = [leading_exponent(f) for f in polys]
    else:
        if len(leading_exponents) != n:
            raise ValueError("The length of leading_exponents must equal the number of polynomials.")
        alphas = [leading_exponent(polys[j], leading_exponents[j]) for j in range(n)]

    det_L, lattice_basis_matrix = lattice_determinant(k, lattice_basis=lattice_basis)

    # Algorithm 1, lines 1--3: Delta_t(f_j) and the data A, b, Delta, e.
    Delta = delta_matrix_from_face_inequalities(P_ineqs, support_sets, alphas)

    # Algorithm 1, lines 4--5: Y and Vert(Y).
    Y = dual_feasible_polyhedron_Y(Delta, e)
    vertices_Y = [vector(QQ, v) for v in Y.vertices_list()]

    # Algorithm 1, lines 6--8: L_i(x) and distinct restrictions varphi_i.
    varphis = distinct_affine_restrictions_from_Y(Y, A, b)

    # Algorithm 1, lines 9--13: cells Q_i and sigma_0.
    cells = []
    sigma0_num = QQ(0)

    for i in range(len(varphis)):
        Q_i = cell_polyhedron(P, i, varphis)
        if (not Q_i.is_empty()) and (Q_i.dim() == k):
            integral_i = integrate_affine_on_polytope(Q_i, varphis[i])
            cells.append({
                "index": i,
                "Q": Q_i,
                "varphi": varphis[i],
                "integral": integral_i
            })
            sigma0_num += integral_i

    # Algorithm 1, lines 14--15: sigma_j = integral_P x_j dx.
    sigma = [integrate_coordinate_on_polytope(P, j) / det_L for j in range(k)]
    sigma0 = sigma0_num / det_L

    if X_names is None:
        X_names = [f"X{j+1}" for j in range(k)]

    lhs = " * ".join([f"{X_names[j]}^({sigma[j]})" for j in range(k)])
    asymptotic_condition = f"{lhs} < {M_name}^({sigma0} - epsilon)"

    theta = common_X_bound_exponent(sigma, sigma0)
    theta_decimal = decimal6(theta)
    common_X_bound = f"X < {M_name}^({theta_decimal} - epsilon')"

    result = {
        "P": P,
        "P_ineqs": P_ineqs,
        "A": A,
        "b": b,
        "T_size": len(P_ineqs),
        "polys": polys,
        "support_sets": support_sets,
        "alphas": alphas,
        "e": e,
        "Delta": Delta,
        "Y": Y,
        "vertices_Y": vertices_Y,
        "varphis": varphis,
        "cells": cells,
        "sigma": sigma,
        "sigma0": sigma0,
        "detL": det_L,
        "lattice_basis_matrix": lattice_basis_matrix,
        "asymptotic_condition": asymptotic_condition,
        "common_X_exponent": theta,
        "common_X_exponent_decimal": theta_decimal,
        "common_X_bound": common_X_bound,

        # Backward-compatible keys from the original script.
        "supports": support_sets,
        "phis": varphis,
        "pieces": cells,
        "I": sigma,
        "C": sigma0,
        "bound": asymptotic_condition
    }

    return result


# Backward-compatible name used by the original script.
def asymptotic_bound_from_face_inequalities(polys,
                                            P_ineqs,
                                            leading_exponents=None,
                                            e=None,
                                            lattice_basis=None,
                                            require_origin=True,
                                            require_nonnegative=True,
                                            X_names=None,
                                            N_name="N"):
    return compute_asymptotic_coefficients_from_face_inequalities(
        polys=polys,
        P_ineqs=P_ineqs,
        leading_exponents=leading_exponents,
        e=e,
        lattice_basis=lattice_basis,
        require_origin=require_origin,
        require_nonnegative=require_nonnegative,
        X_names=X_names,
        M_name=N_name
    )


# ============================================================
#  6. Output
# ============================================================

def print_result(res, varnames=None):
    k = len(res["sigma"])
    if varnames is None:
        varnames = [f"x{j+1}" for j in range(k)]

    print("==================================================")
    print("Input polytope P = {x in R^k : <a_t, x> <= b_t, t in T}:")
    for t, (a, b) in enumerate(res["P_ineqs"], start=1):
        print(f"  t = {t}: <{vector(QQ, a)}, x> <= {QQ(b)}")
    print()

    print("==================================================")
    print("Leading exponents alpha_j:")
    for j, alpha in enumerate(res["alphas"], start=1):
        print(f"  alpha_{j} = {alpha}")
    print()

    print("==================================================")
    print("Algorithm 1 data:")
    print("  A =")
    print(res["A"])
    print("  b =", res["b"])
    print("  e =", res["e"])
    print()

    print("==================================================")
    print("Delta matrix, where Delta[t,j] = max_{u in A(f_j)} <a_t, u - alpha_j>:")
    print(res["Delta"])
    print()

    print("==================================================")
    print("Dual feasible polyhedron:")
    print("  Y = {z in R^{|T|}_{>=0} : Delta^T z >= e}")
    print("  Number of vertices in Vert(Y):", len(res["vertices_Y"]))
    print("  Number of distinct affine restrictions varphi_i:", len(res["varphis"]))
    print()

    print("==================================================")
    print("Piecewise-affine value function phi_P(x) = min_i varphi_i(x):")
    for cell in res["cells"]:
        i = cell["index"]
        varphi_str = format_affine(cell["varphi"], varnames=varnames)
        print(f"  cell Q_{i}:")
        print(f"    varphi_{i}(x) = {varphi_str}")
        print(f"    dim(Q_{i}) = {cell['Q'].dim()}")
        print(f"    integral_Q_{i} varphi_{i}(x) dx = {cell['integral']}")
    print()

    print("==================================================")
    print("Asymptotic coefficients:")
    for j in range(k):
        print(f"  sigma_{j+1} = integral_P {varnames[j]} dx = {res['sigma'][j]}")
    print(f"  sigma_0 = integral_P phi_P(x) dx = {res['sigma0']}")
    print(f"  det(L) = {res['detL']}")
    print()

    print("==================================================")
    print("Asymptotic small-root condition:")
    print(" ", res["asymptotic_condition"])
    print()

    print("==================================================")
    print("If X_1 = ... = X_k = X, then")
    print("  common exponent theta =", res["common_X_exponent_decimal"])
    print(" ", res["common_X_bound"])
    print("==================================================")


# ============================================================
#  7. Paper examples from Section 4
# ============================================================
# The coefficients of the polynomials are irrelevant for Algorithm 1;
# only the support A(f_j), the leading exponent alpha_j, and the modulus
# exponent e_j are used. Therefore all nonzero coefficients are set to 1.
#
# The paper has five application classes. Since MIHNP is treated for
# 3, 4, and 5 samples separately, this script provides seven executable cases:
#   1. CIHNP-CSURF
#   2. MIHNP with 3 samples
#   3. MIHNP with 4 samples
#   4. MIHNP with 5 samples
#   5. ECHNP with 3 samples
#   6. LCG with unknown multiplier
#   7. LIPH for unknown-degree POKE
# ============================================================

import time


def exponent_vector(k, positions=None):
    """
    Construct an exponent vector in Z_{>=0}^k, except that this helper is also
    used for inequality normals, where negative entries may occur.
    Indices are zero-based.
    """
    v = [ZZ(0)] * k
    if positions is None:
        return tuple(v)
    if isinstance(positions, dict):
        for i, a in positions.items():
            v[i] += ZZ(a)
    else:
        for i in positions:
            v[i] += ZZ(1)
    return tuple(v)


def polynomial_from_support(R, support):
    """
    Return a polynomial over R with exactly the prescribed support.
    All nonzero coefficients are set to 1.
    """
    gens = R.gens()
    f = R(0)
    for exp in support:
        if len(exp) != len(gens):
            raise ValueError("An exponent vector has the wrong length.")
        mon = R(1)
        for i, a in enumerate(exp):
            mon *= gens[i] ** ZZ(a)
        f += mon
    return f


def polynomial_ring_with_names(k):
    names = [f"x{i+1}" for i in range(k)]
    return PolynomialRing(QQ, names), names


def nonnegative_inequalities(k):
    """Return inequalities x_i >= 0, encoded as -x_i <= 0."""
    return [(exponent_vector(k, {i: -1}), QQ(0)) for i in range(k)]


def coordinate_upper_inequality(k, i, bound):
    """Return x_i <= bound."""
    return (exponent_vector(k, {i: 1}), QQ(bound))


def linear_inequality(k, coeffs, bound):
    """Return sum_i coeffs[i]*x_i <= bound."""
    return (exponent_vector(k, coeffs), QQ(bound))


def sum_inequality(k, bound):
    """Return x_1 + ... + x_k <= bound."""
    return (tuple([ZZ(1)] * k), QQ(bound))


def pairwise_bilinear_support(k, i, j):
    """Support of x_i*x_j + a*x_i + b*x_j + c in k variables."""
    return [
        exponent_vector(k, [i, j]),
        exponent_vector(k, [i]),
        exponent_vector(k, [j]),
        exponent_vector(k),
    ]


def squarefree_support_degree_at_most(k, subset, max_degree):
    """Squarefree support in variables indexed by subset, total degree <= max_degree."""
    subset = list(subset)
    support = []
    for mask in range(1 << len(subset)):
        if ZZ(mask).popcount() <= max_degree:
            pos = [subset[r] for r in range(len(subset)) if (mask >> r) & 1]
            support.append(exponent_vector(k, pos))
    return support


def print_summary(name, res, expected=None):
    print("==================================================")
    print(name)
    if expected is not None:
        print("  expected from paper:", expected)
    print("  sigma =", res["sigma"])
    print("  sigma0 =", res["sigma0"])
    print("  common exponent theta =", res["common_X_exponent_decimal"])
    print(" ", res["common_X_bound"])


def run_case(name, polys, P_ineqs, leading_exponents, e, varnames,
             expected=None, detailed=False):
    start = time.perf_counter()
    res = compute_asymptotic_coefficients_from_face_inequalities(
        polys=polys,
        P_ineqs=P_ineqs,
        leading_exponents=leading_exponents,
        e=e,
        lattice_basis=None,
        require_origin=True,
        require_nonnegative=True,
        X_names=["X"] * len(varnames),
        M_name="M"
    )
    elapsed = time.perf_counter() - start
    if detailed:
        print_result(res, varnames=varnames)
    else:
        print_summary(name, res, expected=expected)
        if name.startswith("LIPH"):
            # In the LIPH application, the paper fixes X_2 = M^(6/17) and
            # optimizes the exponent of X_1, rather than imposing X_1 = X_2.
            liph_x1_exp = (res["sigma0"] - (QQ(6) / QQ(17)) * res["sigma"][1]) / res["sigma"][0]
            print("  with X2 = M^(6/17): X1 < M^({:.6f} - epsilon')".format(float(liph_x1_exp)))
        print("  running time: {:.6f} seconds".format(elapsed))
    return res


# ------------------------------------------------------------
#  7.1 CIHNP-CSURF
# ------------------------------------------------------------

def make_CIHNP_CSURF_case():
    k = 3
    R, varnames = polynomial_ring_with_names(k)

    supp_f1 = [
        (2, 0, 0), (1, 2, 0), (1, 1, 0), (1, 0, 0),
        (0, 2, 0), (0, 1, 0), (0, 0, 0),
    ]
    supp_f2 = [
        (0, 0, 2), (2, 0, 1), (1, 0, 1), (0, 0, 1),
        (2, 0, 0), (1, 0, 0), (0, 0, 0),
    ]
    polys = [polynomial_from_support(R, supp_f1), polynomial_from_support(R, supp_f2)]

    theta1 = QQ(158) / QQ(100)   # 1.58
    theta2 = QQ(251) / QQ(100)   # 2.51
    theta3 = QQ(419) / QQ(100)   # 4.19
    P_ineqs = nonnegative_inequalities(k) + [
        coordinate_upper_inequality(k, 1, 1),
        linear_inequality(k, {1: 1, 2: 1}, theta1),
        linear_inequality(k, {0: 2, 1: 1}, theta2),
        linear_inequality(k, {0: 2, 1: 1, 2: 4}, theta3),
    ]

    return {
        "name": "CIHNP-CSURF",
        "polys": polys,
        "P_ineqs": P_ineqs,
        "leading_exponents": [(2, 0, 0), (0, 0, 2)],
        "e": [1, 1],
        "varnames": varnames,
        "expected": "theta=(1.58, 2.51, 4.19), X < M^0.2662",
    }


# ------------------------------------------------------------
#  7.2 MIHNP with 3, 4, and 5 samples
# ------------------------------------------------------------

def make_MIHNP_3_samples_case():
    k = 3
    R, varnames = polynomial_ring_with_names(k)
    pairs = [(0, 1), (0, 2), (1, 2)]

    polys = [polynomial_from_support(R, pairwise_bilinear_support(k, i, j)) for i, j in pairs]
    leading_exponents = [exponent_vector(k, [i, j]) for i, j in pairs]

    theta = QQ(3)
    P_ineqs = nonnegative_inequalities(k)
    P_ineqs += [coordinate_upper_inequality(k, i, 1) for i in range(k)]
    P_ineqs += [sum_inequality(k, theta)]

    return {
        "name": "MIHNP with 3 samples",
        "polys": polys,
        "P_ineqs": P_ineqs,
        "leading_exponents": leading_exponents,
        "e": [1, 1, 1],
        "varnames": varnames,
        "expected": "theta=3, X < M^(11/24) approx M^0.4583",
    }


def make_MIHNP_4_samples_case():
    k = 4
    R, varnames = polynomial_ring_with_names(k)
    pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]

    polys = [polynomial_from_support(R, pairwise_bilinear_support(k, i, j)) for i, j in pairs]
    leading_exponents = [exponent_vector(k, [i, j]) for i, j in pairs]
    e = [1] * len(polys)

    supp_4321 = squarefree_support_degree_at_most(k, [0, 1, 2, 3], 3)
    polys.append(polynomial_from_support(R, supp_4321))
    leading_exponents.append((0, 1, 1, 1))
    e.append(2)

    theta = (QQ(315) / QQ(100))
    P_ineqs = nonnegative_inequalities(k)
    P_ineqs += [coordinate_upper_inequality(k, i, 1) for i in range(k)]
    P_ineqs += [sum_inequality(k, theta)]

    return {
        "name": "MIHNP with 4 samples",
        "polys": polys,
        "P_ineqs": P_ineqs,
        "leading_exponents": leading_exponents,
        "e": e,
        "varnames": varnames,
        "expected": "theta=3.15, X < M^0.5336",
    }


def make_MIHNP_5_samples_case():
    k = 5
    R, varnames = polynomial_ring_with_names(k)
    pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]

    polys = [polynomial_from_support(R, pairwise_bilinear_support(k, i, j)) for i, j in pairs]
    leading_exponents = [exponent_vector(k, [i, j]) for i, j in pairs]
    e = [1] * len(polys)

    from itertools import combinations
    for subset in combinations(range(k), 4):
        subset = list(subset)
        polys.append(polynomial_from_support(R, squarefree_support_degree_at_most(k, subset, 3)))
        leading_exponents.append(exponent_vector(k, subset[1:]))
        e.append(2)

    theta = (QQ(340) / QQ(100))
    P_ineqs = nonnegative_inequalities(k)
    P_ineqs += [coordinate_upper_inequality(k, i, 1) for i in range(k)]
    P_ineqs += [sum_inequality(k, theta)]

    return {
        "name": "MIHNP with 5 samples",
        "polys": polys,
        "P_ineqs": P_ineqs,
        "leading_exponents": leading_exponents,
        "e": e,
        "varnames": varnames,
        "expected": "theta=3.40, X < M^0.5794",
    }


# ------------------------------------------------------------
#  7.3 ECHNP with 3 samples
# ------------------------------------------------------------

def make_ECHNP_3_samples_case():
    k = 4
    R, varnames = polynomial_ring_with_names(k)

    polys = []
    leading_exponents = []

    for i in [1, 2, 3]:
        support = [
            exponent_vector(k, {0: 2, i: 1}),
            exponent_vector(k, {0: 1, i: 1}),
            exponent_vector(k, {i: 1}),
            exponent_vector(k, {0: 2}),
            exponent_vector(k, {0: 1}),
            exponent_vector(k),
        ]
        polys.append(polynomial_from_support(R, support))
        leading_exponents.append(exponent_vector(k, {0: 2, i: 1}))

    for i, j in [(1, 2), (1, 3), (2, 3)]:
        support = [
            exponent_vector(k, [i, j]),
            exponent_vector(k, [0, j]),
            exponent_vector(k, [0, i]),
            exponent_vector(k, [j]),
            exponent_vector(k, [i]),
            exponent_vector(k, [0]),
            exponent_vector(k),
        ]
        polys.append(polynomial_from_support(R, support))
        leading_exponents.append(exponent_vector(k, [i, j]))

    theta1, theta2, theta3 = (QQ(140) / QQ(100)), (QQ(232) / QQ(100)), (QQ(236) / QQ(100))
    P_ineqs = nonnegative_inequalities(k)
    P_ineqs += [coordinate_upper_inequality(k, 0, theta1)]
    P_ineqs += [coordinate_upper_inequality(k, i, 1) for i in [1, 2, 3]]
    P_ineqs += [
        linear_inequality(k, {0: 1, 1: 1, 2: 1}, theta2),
        linear_inequality(k, {0: 1, 1: 1, 3: 1}, theta2),
        linear_inequality(k, {0: 1, 2: 1, 3: 1}, theta2),
        linear_inequality(k, {0: 1, 1: 1, 2: 1, 3: 1}, theta3),
    ]

    return {
        "name": "ECHNP with 3 samples",
        "polys": polys,
        "P_ineqs": P_ineqs,
        "leading_exponents": leading_exponents,
        "e": [1] * len(polys),
        "varnames": varnames,
        "expected": "theta=(1.40, 2.32, 2.36), X < M^0.4071",
    }


# ------------------------------------------------------------
#  7.4 LCG with unknown multiplier
# ------------------------------------------------------------

def make_LCG_unknown_multiplier_case(k_lcg=6):
    if k_lcg < 3:
        raise ValueError("LCG requires k_lcg >= 3.")

    k = k_lcg
    R, varnames = polynomial_ring_with_names(k)
    polys = []
    leading_exponents = []

    for i in range(k - 2):
        support = [
            exponent_vector(k, {i + 1: 2}),
            exponent_vector(k, [i, i + 2]),
            exponent_vector(k, [i]),
            exponent_vector(k, [i + 1]),
            exponent_vector(k, [i + 2]),
            exponent_vector(k),
        ]
        polys.append(polynomial_from_support(R, support))
        leading_exponents.append(exponent_vector(k, {i + 1: 2}))

    P_ineqs = nonnegative_inequalities(k) + [sum_inequality(k, 1)]

    return {
        "name": f"LCG with unknown multiplier, k={k_lcg}",
        "polys": polys,
        "P_ineqs": P_ineqs,
        "leading_exponents": leading_exponents,
        "e": [1] * len(polys),
        "varnames": varnames,
        "expected": "X < M^(1/2 - 1/k); default k=6 gives X < M^(1/3)",
    }


# ------------------------------------------------------------
#  7.5 LIPH for unknown-degree POKE
# ------------------------------------------------------------

def make_LIPH_unknown_degree_case():
    k = 2
    R, varnames = polynomial_ring_with_names(k)

    support = [(2, 1), (1, 1), (0, 1), (0, 0)]
    polys = [polynomial_from_support(R, support)]

    # Algorithm 1 is implemented for rational polytopes.
    # The paper states that an optimal rational choice is close to
    # (2*sqrt(223)-10)/33 ≈ 0.6020. We use 301/500 = 0.602.
    theta = (QQ(301) / QQ(500))
    P_ineqs = nonnegative_inequalities(k) + [
        coordinate_upper_inequality(k, 1, 1),
        linear_inequality(k, {0: 1, 1: -2}, theta),
    ]

    return {
        "name": "LIPH for unknown-degree POKE",
        "polys": polys,
        "P_ineqs": P_ineqs,
        "leading_exponents": [(2, 1)],
        "e": [1],
        "varnames": varnames,
        "expected": "theta approx (2*sqrt(223)-10)/33 approx 0.6020; X1 < M^0.20195 when X2=M^(6/17)",
    }


# ------------------------------------------------------------
#  7.6 Run all paper examples
# ------------------------------------------------------------

def paper_example_cases(include_lcg=True, include_liph=True, k_lcg=6):
    """
    Return the cases from the paper.

    The first five entries correspond to the concrete optimized cases displayed
    before LCG/LIPH: CIHNP, MIHNP-3, MIHNP-4, MIHNP-5, and ECHNP.
    If include_lcg/include_liph are True, the remaining two application classes
    are also included.
    """
    cases = [
        make_CIHNP_CSURF_case(),
        make_MIHNP_3_samples_case(),
        make_MIHNP_4_samples_case(),
        make_MIHNP_5_samples_case(),
        make_ECHNP_3_samples_case(),
    ]
    if include_lcg:
        cases.append(make_LCG_unknown_multiplier_case(k_lcg=k_lcg))
    if include_liph:
        cases.append(make_LIPH_unknown_degree_case())
    return cases


def run_paper_examples(include_lcg=True, include_liph=True, k_lcg=6, detailed=False):
    """Run Algorithm 1 on the examples from Section 4."""
    total_start = time.perf_counter()
    results = {}
    for case in paper_example_cases(include_lcg=include_lcg, include_liph=include_liph, k_lcg=k_lcg):
        results[case["name"]] = run_case(
            name=case["name"],
            polys=case["polys"],
            P_ineqs=case["P_ineqs"],
            leading_exponents=case["leading_exponents"],
            e=case["e"],
            varnames=case["varnames"],
            expected=case["expected"],
            detailed=detailed,
        )
    total_elapsed = time.perf_counter() - total_start
    print("==================================================")
    print("Total running time for all selected paper examples: {:.6f} seconds".format(total_elapsed))
    print("==================================================")
    return results


if __name__ == "__main__":
    # By default this runs the first five optimized cases plus LCG and LIPH.
    # To run only the first five concrete cases, use:
    #     run_paper_examples(include_lcg=False, include_liph=False)
    results = run_paper_examples(include_lcg=True, include_liph=True, k_lcg=6, detailed=False)
