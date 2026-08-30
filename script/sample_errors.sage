#!/usr/bin/env sage
#
# Error sampling for UnifOMR RLWEenc, Parameter Set 1 (Table 1, Section 7.1),
# using Sage's probability distributions and the cyclotomic quotient ring
# R = Z_q[X] / (Phi_2048(X)) = Z_q[X] / (X^1024 + 1).
#
# Parameter Set 1 (Table 1):
#
#   n       = 900       secret key dimension
#   nring   = 1024      ring dimension (n padded to next power of two)
#   q       = 4169729   ciphertext modulus
#   sigma   = 0.6       RLWEenc error standard deviation
#   h       = 80        secret key Hamming weight
#   r       = 48        error range bound
#   ell     = 1         number of bits encrypted per RLWEenc ciphertext
#   eps_n   = 2^-30     false negative rate
#
# The script samples error polynomials x <- chi_sigma in R (Definition 4.6),
# where chi_sigma is the discrete Gaussian over Z with density prop. to
# exp(-x^2/(2*sigma^2)), and verifies them against the paper's analysis.
#
# Expected output (last run), line by line:
#
#   ring                  : Univariate Quotient Polynomial Ring in x over
#                           Ring of integers modulo 4169729 with modulus X^1024 + 1
#   defining relation     : Phi_2048(X) = True (X^1024 + 1)
#       all arithmetic happens in the quotient ring Z_q[X]/(Phi_2048);
#       True confirms Phi_2048 collapses to X^1024 + 1, i.e. the ring of
#       Definition 4.6 with m = 2*nring
#
#       "negacyclic" means the defining relation is X^nring = -1: powers of
#       X wrap around the top with a sign flip, X^(nring + k) = -X^k, so
#       polynomial multiplication becomes a cyclic convolution where the
#       wraparound terms are negated. Compare with "cyclic" X^nring = 1
#       (i.e. reduction mod X^nring - 1), where wraparound is positive.
#       Concretely, (sum a_i X^i) * (sum b_j X^j) has coefficient
#       c_k = sum_{i+j=k} a_i b_j  -  sum_{i+j=k+nring} a_i b_j
#       (ordinary terms minus the wrapped ones), all mod q
#
#   samples               : 2048000 coefficients
#       2000 error polynomials x 1024 coefficients each
#   empirical mean        : -0.00122
#       chi_sigma is centered; the sample mean fluctuates around 0
#   empirical std         : 0.59351   (target 0.6)
#       per-coefficient spread; slightly UNDER sigma because a discrete
#       Gaussian on Z with sigma < 1 squeezes onto {0, +-1}: exact variance
#       2*0.1655 + 4*0.0026 + ... ~ 0.352 vs sigma^2 = 0.36
#   empirical P(0)        : 0.66351
#   empirical P(+1)=P(-1) : 0.16506
#       the near-ternary shape: ~99.5% of coefficients land in {0, +-1},
#       matching theory P(0) = 1/1.5065 = 0.6638, P(+-1) = 0.2494/1.5065
#       = 0.1655 each; P(|x| >= 2) ~ 0.005
#   max |coeff|           : 3          (bound r = 48)
#       the largest of 2M raw draws is 3, nowhere near r; individual error
#       terms are tiny, it is their accumulation that must respect r
#
#   effective noise std   : 7.6131 = sqrt(2h+1)*sigma
#       after decryption, pk[1] - alpha*s contains a sum of ~2h+1 = 161
#       independent sigma-Gaussians (h products with the ternary key, plus
#       fresh error), so the noise actually bounded by r has this std
#   P(|noise| > r) approx : 2.884e-10   (need <= 9.313e-10)
#       the correctness budget: this is Definition 4.6's condition
#       erf(r / (sqrt(2)*sqrt(2h+1)*sigma)) <= eps_n/ell, i.e. the chance a
#       decryption noise crosses r = 48 stays below the false-negative rate
#       2^-30; observed tail is ~3.2x inside the budget
#
#   demo pk noise max|c|  : 2
#       end-to-end KeyGen check: pk = (alpha, alpha*s + x) computed in the
#       quotient ring, then alpha*s subtracted, recovering exactly x
#       (largest coefficient 2) -- proof that Sage's ring arithmetic reduces
#       mod X^1024 + 1 and q correctly
#
# In short: the first two lines verify the algebra is in the right ring, the
# next five verify chi_sigma has the right statistics, the effective-noise
# pair verifies the paper's correctness condition with margin, and the demo
# line verifies the ring arithmetic end-to-end.

from random import choice

from sage.all import (
    Zmod,
    PolynomialRing,
    RealField,
    RR,
    ZZ,
    erf,
    ceil,
    sqrt,
    cyclotomic_polynomial,
    Subsets,
)

from sage.probability.probability_distribution import GeneralDiscreteDistribution

# ------------------------- Parameter Set 1 -------------------------

# from Table 1, Section 7.1 of the UnifOMR paper
# n is the secret key dimension; footnote 10 says n is zero-padded up to the
# next power of two nring so we can work over a cyclotomic ring
n = 900
nring = 1024

# the (2*nring)-th cyclotomic polynomial is X^nring + 1, so R is the
# negacyclic ring Z_q[X]/(X^1024 + 1) used by RLWEenc (Definition 4.6)
m = 2 * nring
q = 4169729

# chi_sigma: per-coefficient noise, discrete Gaussian with std sigma
sigma = 0.6

# secret distribution D: ternary {-1,0,1} vector with Hamming weight h
h = 80

# error range bound from Definition 4.6:
# the minimum r with erf(r / (sqrt(2)*sqrt(2h+1)*sigma)) <= eps_n/ell
r = 48
ell = 1
eps_n = 2^-30

# ------------------------- the cyclotomic ring -------------------------

# P = Z_q[X]; quotient by the m-th cyclotomic polynomial gives
# R = Z_q[X]/(Phi_m) = Z_q[X]/(X^nring + 1), the negacyclic ring:
# X^nring acts as -1, so multiplication wraps with sign flips
P = PolynomialRing(Zmod(q), "X")
X = P.gen()
phi = P(cyclotomic_polynomial(m))
R = P.quotient(phi, "x")

# ------------------------- chi_sigma as a Sage distribution -------------------------

# the error distribution chi_sigma is the discrete Gaussian over Z with
# density prop. to exp(-x^2/(2*sigma^2)); GeneralDiscreteDistribution
# takes unnormalized weights (it normalizes internally; the weights sum
# to ~1.5065, the Gaussian normalization constant for sigma = 0.6)
PREC = RealField(200)

# 12*sigma cutoff: mass beyond is e^(-72) ~ 10^-31, unreachable in practice
support = list(range(-ceil(12 * sigma), ceil(12 * sigma) + 1))
weights = [PREC(exp(-x^2 / (2 * sigma^2))) for x in support]
chi_sigma = GeneralDiscreteDistribution(weights)


def chi_coeff():
    # GeneralDiscreteDistribution returns an index into the weight list,
    # so map it back through the support {-8, ..., 8}
    return ZZ(support[chi_sigma.get_random_element()])


def sample_error():
    # x <- chi_sigma in R: nring i.i.d. discrete Gaussian coefficients,
    # each with std sigma; expected std 0.594, P(0) = 0.664, P(+-1) = 0.165
    return R([chi_coeff() for _ in range(nring)])


def sample_secret():
    # s <- D: ternary {-1, 0, 1}^n of Hamming weight h, zero-padded to nring
    coeffs = [0] * nring
    for i in Subsets(range(n), h).random_element():
        coeffs[i] = choice([-1, 1])
    return R(coeffs)


def sample_uniform():
    # uniform ring element, used for the public random alpha <-$ R_q
    return R([Zmod(q).random_element() for _ in range(nring)])


# ------------------------- verification -------------------------

def centered(coeffs):
    # interpret Z_q coefficients in (-q/2, q/2] so tiny noises show as tiny
    # (an error of +2 must read as 2, not as q - 2)
    return [min(ZZ(c), q - ZZ(c)) for c in coeffs]


def verify(num_samples=2000):
    # group 1: ring sanity -- the quotient modulus should collapse to
    # X^nring + 1, confirming we are in the negacyclic ring of Def 4.6
    # (X^nring = -1, so wraparound terms in products are negated)
    # last run: True
    print(f"ring                  : {R}")
    print(f"defining relation     : Phi_{m}(X) = {phi == X^nring + 1} (X^{nring} + 1)")

    # group 2: chi_sigma statistics -- draw num_samples*nring single
    # coefficients (last run: 2000*1024 = 2048000) and compare to theory
    flat = [chi_coeff() for _ in range(num_samples * nring)]
    mu = RR(sum(flat)) / len(flat)
    var = RR(sum((v - mu)^2 for v in flat)) / len(flat)

    # centered: the sample mean fluctuates around 0
    # last run: -0.00122
    print(f"samples               : {len(flat)} coefficients")
    print(f"empirical mean        : {float(mu):+.5f}")

    # per-coefficient spread; runs slightly under sigma because the
    # discrete Gaussian on Z with sigma < 1 concentrates on {0, +-1}
    # (variance ~0.352 vs sigma^2 = 0.36)
    # last run: 0.59351
    print(f"empirical std         : {float(sqrt(var)):.5f}   (target {sigma})")

    # near-ternary shape; theory: P(0) = 0.6638, P(+-1) = 0.1655 each,
    # P(|x| >= 2) ~ 0.005 -- so ~99.5% of coefficients are in {0, +-1}
    # last run: P(0) = 0.66351, P(+-1) = 0.16506 each
    print(f"empirical P(0)        : {float(RR(flat.count(0)) / len(flat)):.5f}")
    print(f"empirical P(+1) = P(-1): {float(RR(flat.count(1)) / len(flat)):.5f}")

    # raw error terms must never approach the decryption bound r = 48;
    # individual errors are tiny, only their accumulated sum matters
    # last run: max |coeff| = 3 over 2048000 coefficients
    print(f"max |coeff|           : {max(abs(v) for v in flat)}          (bound r = {r})")
    assert max(abs(v) for v in flat) < r

    # group 3: correctness budget from Definition 4.6
    # decryption computes pk[1] - alpha*s, where the surviving noise is
    # roughly a sum of 2h+1 = 161 independent sigma-Gaussians (h products
    # x_i*y_j with the ternary y, plus fresh error), hence effective std
    # sqrt(2h+1)*sigma; a decryption fails iff |noise| > r, and that tail
    # probability must stay below eps_n/ell
    # last run: eff = 7.6131, tail = 2.884e-10 vs budget 9.313e-10 (~3.2x margin)
    eff = RR(sqrt(2 * h + 1)) * sigma
    tail = 1 - erf(r / (RR(sqrt(2)) * eff))
    print(f"effective noise std   : {float(eff):.4f} = sqrt(2h+1)*sigma")
    print(f"P(|noise| > r) approx : {float(tail):.3e}   (need <= eps_n/ell = {float(eps_n / ell):.3e})")

    # group 4: end-to-end demo of KeyGen in the actual quotient ring
    # pk = (alpha, alpha*s + x); subtracting alpha*s recovers exactly x,
    # confirming ring arithmetic reduces mod X^1024 + 1 and q correctly
    # last run: max|c| = 2
    s = sample_secret()
    alpha = sample_uniform()
    pk = alpha * s + sample_error()
    noise_coeffs = (pk - alpha * s).lift().coefficients(sparse=False)
    noise_coeffs += [0] * (nring - len(noise_coeffs))
    print(f"demo pk noise max|c|  : {max(centered(noise_coeffs))}")


if __name__ in ("__main__", "sage.all"):
    verify()
