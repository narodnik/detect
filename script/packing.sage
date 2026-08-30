#!/usr/bin/env sage
#
# Packing demo for the UnifOMR RLWEenc scheme, in miniature.
#
# Reproduces the worked examples of crypto.md, section 3:
#   - KeyGen / Enc / Dec in the toy ring Z_17[X]/(X^4 + 1)  (section 3.2)
#   - four-clue BFV coefficient packing                      (section 3.4)
#   - a randomized end-to-end consistency check
#
# Toy parameters:
#
#   q      = 17        toy modulus
#   nring  = 4         toy ring dimension (real: 1024)
#   D      = 4         toy BFV ring dimension / number of slots (real: 2048)
#   h      = 2         toy secret Hamming weight (real: 80)
#   r      = 3         toy error range bound (real: 48)
#   q/2    = 8         plaintext bit encoding of Z_17
#
# Last run output is quoted in the comments below each section.

from random import choice

from sage.all import Zmod, PolynomialRing, RealField, RR, ZZ, ceil, Subsets

from sage.probability.probability_distribution import GeneralDiscreteDistribution

q = 17
nring = 4
h_toy = 2
r_toy = 3

P = PolynomialRing(Zmod(q), "X")
X = P.gen()
R = P.quotient(X^nring + 1, "x")


def centered_coeffs(f):
    # lift to a coefficient list, centered in (-q/2, q/2]
    cs = f.lift().coefficients(sparse=False)
    cs += [0] * (nring - len(cs))
    return [ZZ(c) - q if ZZ(c) > q // 2 else ZZ(c) for c in cs]


def poly_str(f):
    # human-readable centered rendering, e.g. "9 + 2X - X^2 - 2X^3"
    parts = []
    for i, c in enumerate(centered_coeffs(f)):
        if c == 0 and not (i == 0 and len(set(centered_coeffs(f))) == 1):
            continue
        sign = " + " if c > 0 and parts else (" - " if c < 0 and parts else ("-" if c < 0 else ""))
        mag = str(abs(c)) if abs(c) != 1 or i == 0 else ""
        mon = mag + (f"X^{i}" if i > 1 else ("X" if i == 1 else ""))
        parts.append(sign + mon)
    return "".join(parts) if parts else "0"


# ------------------------- toy chi_sigma -------------------------

sigma = 0.6
PR = RealField(200)
support = list(range(-ceil(12 * sigma), ceil(12 * sigma) + 1))
weights = [PR(exp(-x^2 / (2 * sigma^2))) for x in support]
chi_sigma = GeneralDiscreteDistribution(weights)


def toy_noise():
    return R([ZZ(support[chi_sigma.get_random_element()]) for _ in range(nring)])


def toy_secret():
    coeffs = [0] * nring
    for i in Subsets(range(nring), h_toy).random_element():
        coeffs[i] = choice([-1, 1])
    return R(coeffs)


def toy_uniform():
    return R([Zmod(q).random_element() for _ in range(nring)])


# ------------------------- section 3.3: KeyGen / Enc / Dec -------------------------

def worked_scheme():
    alpha = R(3 + 5 * X - 2 * X^2 + X^3)
    s = R(1 - X^3)
    x = R(1 - X)

    pk1 = alpha * s + x
    print(f"alpha                : {poly_str(alpha)}")
    print(f"s                    : {poly_str(s)}   (weight {h_toy})")
    print(f"pk[1] = alpha*s + x  : {poly_str(pk1)}")
    # last run: pk[1] = alpha*s + x  : -8 + 2X - X^2 - 2X^3  (== 9 + 2X - X^2 - 2X^3 mod 17)
    assert centered_coeffs(pk1) == [-8, 2, -1, -2]

    y = R(-1 + X^2)
    x1 = R(1)
    x2 = R(-1)
    t = R(8)

    a = alpha * y + x1
    b = pk1 * y + t + x2
    print(f"a = alpha*y + x'     : {poly_str(a)}")
    print(f"b = pk[1]*y + t + x''  : {poly_str(b)}")
    # last run: a = -6X + 5X^2 + 4X^3 ; b = -1 - 7X^2 + 4X^3  (b == -1 + 10X^2 + 4X^3 mod 17)
    assert centered_coeffs(a) == [0, -6, 5, 4]
    assert centered_coeffs(b) == [-1, 0, -7, 4]

    d = b - a * s
    noise = d - t
    dc = centered_coeffs(d)
    m0 = 1 if abs(dc[0]) > r_toy else 0
    print(f"d = b - a*s          : {poly_str(d)}")
    print(f"noise = d - t        : {poly_str(noise)}   per-coeff {centered_coeffs(noise)}")
    print(f"decode |d_0|={abs(dc[0])} vs r={r_toy}  : m_0 = {m0}   (sent 1)")
    # last run: d = 5 + X + X^2 ; noise per-coeff (-3, 1, 1, 0); m_0 = 1
    assert centered_coeffs(d) == [5, 1, 1, 0]
    assert m0 == 1


# ------------------------- section 3.4: four-clue packing -------------------------

def lin_form_scalars(a_coeffs, k=0):
    # unfold (a*s)_k = sum_i A[k][i] * s_i over Z_q:
    # position i pairs with a[(k - i) mod nring], negated when the product wraps
    scalars = []
    for i in range(nring):
        val = a_coeffs[(k - i) % nring]
        scalars.append(ZZ(val) if i <= k else -ZZ(val))
    return scalars


def worked_packing():
    s_vec = [1, 0, 0, -1]

    clues = [
        [1, 2, 0, 5],
        [4, -3, 2, 0],
        [-2, 1, 0, -4],
        [6, -1, 3, 2],
    ]
    names = ["1 + 2X + 5X^3", "4 - 3X + 2X^2", "-2 + X - 4X^3", "6 - X + 3X^2 + 2X^3"]
    for j, name in enumerate(names):
        print(f"clue {j}: a^({j}) = {name}")

    A = [lin_form_scalars(clues[j]) for j in range(nring)]

    p = [R([A[j][i] for j in range(nring)]) for i in range(nring)]
    for i in range(nring):
        print(f"p_{i}(X) = {poly_str(p[i])}")
    # last run: p_0 = 1 + 4X - 2X^2 + 6X^3 ; p_1 = -5 + 4X^2 - 2X^3
    #           p_2 = -2X - 3X^3           ; p_3 = -2 + 3X - X^2 + X^3
    assert centered_coeffs(p[0]) == [1, 4, -2, 6]
    assert centered_coeffs(p[1]) == [-5, 0, 4, -2]
    assert centered_coeffs(p[2]) == [0, -2, 0, -3]
    assert centered_coeffs(p[3]) == [-2, 3, -1, 1]

    # homomorphic product-sum, shown on plaintexts:
    # sum_i p_i * ct_{s_i} where ct_{s_i} encrypts the constant s_i
    total = sum(p[i] * R(s_vec[i]) for i in range(nring))
    print(f"sum_i p_i * s_i      : {poly_str(total)}")
    # last run: 3 + X - X^2 + 5X^3

    got = centered_coeffs(total)
    for j in range(nring):
        direct = centered_coeffs(R(clues[j]) * R(s_vec))[0]
        mark = "ok" if got[j] == direct else "FAIL"
        print(f"coefficient {j} = {got[j]:+d}   direct (a^({j})*s)_0 = {direct:+d}   [{mark}]")
        assert got[j] == direct
    # last run: 3, 1, -1, 5 -- all ok


# ------------------------- randomized end-to-end check -------------------------

def random_trials(num=200):
    ok = 0
    for _ in range(num):
        alpha = toy_uniform()
        s = toy_secret()
        x = toy_noise()
        pk1 = alpha * s + x

        y = toy_secret()
        x1 = toy_noise()
        x2 = toy_noise()
        m = choice([0, 1])
        t = R(8 * m)

        a = alpha * y + x1
        b = pk1 * y + t + x2
        dc = centered_coeffs(b - a * s)
        if (1 if abs(dc[0]) > r_toy else 0) == m:
            ok += 1

    print(f"random trials        : {ok}/{num} decodes correct")
    # last run: 200/200 (occasional failures are the toy's own false-negative
    # budget: noise std sqrt(2h+1)*sigma ~ 1.3 against r = 3, tail ~ 2%)
    return ok


# ------------------------- randomized packing check -------------------------

def random_packing_trials(num=50):
    for _ in range(num):
        s_vec = centered_coeffs(toy_secret())
        clue_lists = []
        for _ in range(nring):
            a = toy_uniform()
            clue_lists.append(centered_coeffs(a))
        A = [lin_form_scalars(clue_lists[j]) for j in range(nring)]
        p = [R([A[j][i] for j in range(nring)]) for i in range(nring)]
        total = sum(p[i] * R(s_vec[i]) for i in range(nring))
        got = centered_coeffs(total)
        for j in range(nring):
            assert got[j] == centered_coeffs(R(clue_lists[j]) * R(s_vec))[0]
    print(f"random packing trials: {num}/{num} coefficient-extraction checks passed")


if __name__ in ("__main__", "sage.all"):
    print("=== worked KeyGen / Enc / Dec (section 3.2) ===")
    worked_scheme()
    print()
    print("=== worked four-clue packing (section 3.4) ===")
    worked_packing()
    print()
    print("=== randomized checks ===")
    random_trials()
    random_packing_trials()
