//! The finite field $\mathbb{F}_q$ and the negacyclic polynomial ring
//! $R_q = \mathbb{F}_q[X]/(X^{N} + 1)$ over it.
//!
//! The field is generic over its (prime) modulus; the ring is generic over
//! modulus and dimension, with multiplication implementing the negacyclic
//! convolution: with $X^N = -1$, coefficient $k$ of a product is
//!
//! $$c_k = \sum_{i+j=k} a_i b_j \;-\; \sum_{i+j=k+N} a_i b_j \pmod q .$$
//!
//! The type alias [`Rq`] instantiates the ring at Parameter Set 1.

use std::ops::{Add, Mul, Neg, Sub};

use crate::param;

/// An element of the prime finite field $\mathbb{F}_p = \mathbb{Z}/p\mathbb{Z}$.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Hash)]
pub struct Fp<const P: u64>(u64);

impl<const P: u64> Fp<P> {
    /// Reduce a `u64` into the field.
    pub fn new(v: u64) -> Self {
        debug_assert!(P > 1);
        Self(v % P)
    }

    /// Reduce a (possibly negative) integer into the field, e.g. the
    /// centered noise coefficients.
    pub fn from_i64(v: i64) -> Self {
        Self(v.rem_euclid(P as i64) as u64)
    }

    pub fn zero() -> Self {
        Self(0)
    }

    pub fn one() -> Self {
        Self(1 % P)
    }

    pub fn is_zero(&self) -> bool {
        self.0 == 0
    }

    /// Scalar exponentiation by square-and-multiply.
    pub fn pow(self, mut e: u64) -> Self {
        let mut result = Self::one();
        let mut base = self;
        while e > 0 {
            if !e.is_multiple_of(2) {
                result = result * base;
            }
            base = base * base;
            e >>= 1;
        }
        result
    }

    /// Multiplicative inverse via Fermat's little theorem; only valid
    /// when $p$ is prime.
    pub fn invert(self) -> Option<Self> {
        if self.0 == 0 {
            return None;
        }
        Some(self.pow(P - 2))
    }

    /// Lift into $(-p/2, p/2]$, the natural view for small noise values.
    pub fn lift_centered(self) -> i64 {
        if self.0 <= P / 2 {
            self.0 as i64
        } else {
            self.0 as i64 - P as i64
        }
    }

    pub fn raw(self) -> u64 {
        self.0
    }
}

impl<const P: u64> Add for Fp<P> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self((self.0 + rhs.0) % P)
    }
}

impl<const P: u64> Sub for Fp<P> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self((self.0 + P - rhs.0) % P)
    }
}

impl<const P: u64> Mul for Fp<P> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Self((self.0 as u128 * rhs.0 as u128 % P as u128) as u64)
    }
}

impl<const P: u64> Neg for Fp<P> {
    type Output = Self;
    fn neg(self) -> Self {
        Self((P - self.0) % P)
    }
}

/// A polynomial of the negacyclic ring
/// $R = \mathbb{F}_p[X]/(X^N + 1)$, stored as its $N$ coefficients.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct PolyFp<const P: u64, const N: usize>([u64; N]);

impl<const P: u64, const N: usize> PolyFp<P, N> {
    /// The zero polynomial.
    pub fn zero() -> Self {
        Self([0; N])
    }

    /// Build from (possibly negative, possibly fewer than $N$) integer
    /// coefficients, reducing each into $\mathbb{F}_p$ and zero-padding.
    pub fn from_i64_coeffs(coeffs: &[i64]) -> Self {
        debug_assert!(coeffs.len() <= N);
        let mut arr = [0u64; N];
        for (a, c) in arr.iter_mut().zip(coeffs) {
            *a = Fp::<P>::from_i64(*c).raw();
        }
        Self(arr)
    }

    /// Construct directly from already-reduced $\mathbb{F}_p$
    /// coefficients.
    pub fn from_raw(coeffs: [u64; N]) -> Self {
        Self(coeffs)
    }

    pub fn coeff(&self, i: usize) -> Fp<P> {
        Fp(self.0[i])
    }

    /// Coefficient read out in $(-p/2, p/2]$, e.g. for $|d_i| \le r$
    /// decoding checks.
    pub fn coeff_centered(&self, i: usize) -> i64 {
        self.coeff(i).lift_centered()
    }

    /// Multiply every coefficient by a scalar.
    pub fn scalar_mul(&self, s: Fp<P>) -> Self {
        let arr = self.0.map(|c| (Fp::<P>::new(c) * s).raw());
        Self(arr)
    }

    /// Negacyclic convolution: multiply modulo $X^N + 1$, where wrapped
    /// terms flip sign. Schoolbook $O(N^2)$, accumulated in `i128`
    /// (requires $N \cdot p^2 < 2^{127}$, comfortably true for the
    /// Parameter Set 1 moduli).
    pub fn mul(&self, other: &Self) -> Self {
        let mut acc = [0i128; N];
        for (i, a_val) in self.0.iter().enumerate() {
            if *a_val == 0 {
                continue;
            }
            let a = *a_val as i128;
            for (j, b_val) in other.0.iter().enumerate() {
                if *b_val == 0 {
                    continue;
                }
                let prod = a * *b_val as i128;
                let k = i + j;
                if k < N {
                    acc[k] += prod;
                } else {
                    acc[k - N] -= prod;
                }
            }
        }
        let arr = acc.map(|v| v.rem_euclid(P as i128) as u64);
        Self(arr)
    }

    /// Negacyclic multiplication where `t` is ternary (every
    /// coefficient in $\{-1, 0, 1\}$): each nonzero $t_i = \pm 1$ adds
    /// or subtracts one rotated copy of `self`, skipping the $t_i = 0$
    /// positions entirely — for a weight-$h$ secret this is $h \cdot N$
    /// additions instead of $N^2$ multiplications.
    ///
    /// Accumulated in `i128` without intermediate reduction (requires
    /// $N \cdot p < 2^{127}$, so any modulus up to $2^{116}$; in
    /// particular the BFV modulus $2^{60}$).
    pub fn mul_ternary(&self, t: &Self) -> Self {
        let mut acc = [0i128; N];
        for (i, t_val) in t.0.iter().enumerate() {
            let sign: i128 = if *t_val == 0 {
                continue;
            } else if *t_val == 1 {
                1
            } else if *t_val == P - 1 {
                -1
            } else {
                panic!("mul_ternary: coefficient {i} = {t_val} is not ternary")
            };
            for (j, p_val) in self.0.iter().enumerate() {
                let prod = sign * (*p_val as i128);
                let k = i + j;
                if k < N {
                    acc[k] += prod;
                } else {
                    acc[k - N] -= prod;
                }
            }
        }
        let arr = acc.map(|v| v.rem_euclid(P as i128) as u64);
        Self(arr)
    }

    /// Negacyclic multiplication where every coefficient of `plain` is
    /// canonically small ($< 2^{44}$) while `self` is arbitrary: the
    /// plaintext-by-ciphertext product of the BFV layer, where `plain`
    /// carries RLWEenc scalars lifted from $\mathbb{Z}_q$ ($q < 2^{22}$).
    ///
    /// Unlike [`PolyFp::mul`] this needs no per-product reduction even
    /// at the BFV modulus $P = 2^{60}$: the whole convolution sums to
    /// at most $N^2 \cdot 2^{44} \cdot P < 2^{127}$ in `i128`. Zero
    /// coefficients of `plain` are skipped, so packing $D$ clues into
    /// fewer coefficients costs proportionally less.
    pub fn mul_plain(&self, plain: &Self) -> Self {
        debug_assert!(plain.0.iter().all(|c| *c < (1u64 << 44)));
        let mut acc = [0i128; N];
        for (i, p_val) in plain.0.iter().enumerate() {
            if *p_val == 0 {
                continue;
            }
            let p = *p_val as i128;
            for (j, b_val) in self.0.iter().enumerate() {
                if *b_val == 0 {
                    continue;
                }
                let prod = p * (*b_val as i128);
                let k = i + j;
                if k < N {
                    acc[k] += prod;
                } else {
                    acc[k - N] -= prod;
                }
            }
        }
        let arr = acc.map(|v| v.rem_euclid(P as i128) as u64);
        Self(arr)
    }

    /// Reinterpret the coefficients in $\mathbb{F}_{P_2}[X]/(X^N + 1)$,
    /// re-centering each one first (requires every coefficient to fit
    /// canonically, $|v| \le P_2/2$) — e.g. packing RLWEenc scalars
    /// from $R_q$ into the BFV ring, or re-embedding a ternary secret
    /// at the modulus-switched modulus $Q'$.
    pub fn lift_to<const P2: u64>(&self) -> PolyFp<P2, N> {
        let mut out = [0u64; N];
        for (o, i) in out.iter_mut().zip(0..N) {
            let centered = self.coeff(i).lift_centered();
            debug_assert!(centered.abs() <= (P2 as i64) / 2);
            *o = Fp::<P2>::from_i64(centered).raw();
        }
        PolyFp::<P2, N>::from_raw(out)
    }

    /// Largest absolute centered coefficient, e.g. for noise readouts.
    pub fn max_centered_abs(&self) -> i64 {
        self.0
            .iter()
            .map(|c| Fp::<P>::new(*c).lift_centered().abs())
            .max()
            .unwrap_or(0)
    }
}

impl<const P: u64, const N: usize> Add for &PolyFp<P, N> {
    type Output = PolyFp<P, N>;
    fn add(self, rhs: Self) -> PolyFp<P, N> {
        let arr = std::iter::zip(self.0.iter(), rhs.0.iter())
            .map(|(a, b)| Fp::<P>::new(*a + *b).raw())
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        PolyFp(arr)
    }
}

impl<const P: u64, const N: usize> Sub for &PolyFp<P, N> {
    type Output = PolyFp<P, N>;
    fn sub(self, rhs: Self) -> PolyFp<P, N> {
        let arr = std::iter::zip(self.0.iter(), rhs.0.iter())
            .map(|(a, b)| Fp::<P>::new((*a + P) - *b).raw())
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        PolyFp(arr)
    }
}

impl<const P: u64, const N: usize> Neg for &PolyFp<P, N> {
    type Output = PolyFp<P, N>;
    fn neg(self) -> PolyFp<P, N> {
        let arr = self
            .0
            .iter()
            .map(|c| Fp::<P>::new(*c).neg().raw())
            .collect::<Vec<_>>();
        PolyFp(arr.try_into().unwrap())
    }
}

impl<const P: u64, const N: usize> Mul for &PolyFp<P, N> {
    type Output = PolyFp<P, N>;
    fn mul(self, rhs: Self) -> PolyFp<P, N> {
        PolyFp::mul(self, rhs)
    }
}

/// The RLWEenc ring $R_q = \mathbb{F}_q[X]/(X^{n'} + 1)$ at Parameter
/// Set 1.
pub type Rq = PolyFp<{ param::Q }, { param::N_RING }>;

#[cfg(test)]
mod tests {
    use super::*;

    type Toy = PolyFp<17, 4>;

    #[test]
    fn field_arithmetic() {
        let a = Fp::<17>::from_i64(-1);
        assert_eq!(a.raw(), 16);
        assert_eq!(a.lift_centered(), -1);
        let b = Fp::<17>::from_i64(5);
        assert_eq!((a + b).lift_centered(), 4);
        assert_eq!((a * b).lift_centered(), -5);
        assert_eq!((-a).lift_centered(), 1);
        assert_eq!(b.invert().unwrap() * b, Fp::one());
        assert!(Fp::<17>::zero().invert().is_none());
        assert_eq!(Fp::<7>::from_i64(2).pow(3).raw(), 1);
    }

    #[test]
    fn negacyclic_wrap_flips_signs() {
        let x = Toy::from_i64_coeffs(&[0, 1]);
        let sq = &x * &x;
        assert_eq!(sq.coeff_centered(2), 1);
        let x3 = Toy::from_i64_coeffs(&[0, 0, 0, 1]);
        let wrapped = &x * &x3;
        assert_eq!(wrapped.coeff_centered(0), -1);
        let a = Toy::from_i64_coeffs(&[3, 5, -2, 1]);
        let b = Toy::from_i64_coeffs(&[1, 0, 0, -1]);
        let ab = &a * &b;
        assert_eq!(ab.coeff_centered(0), 8);
        assert_eq!(ab.coeff_centered(1), 3);
        assert_eq!(ab.coeff_centered(2), -1);
        assert_eq!(ab.coeff_centered(3), -2);
    }

    #[test]
    fn worked_keygen_enc_dec() {
        let alpha = Toy::from_i64_coeffs(&[3, 5, -2, 1]);
        let s = Toy::from_i64_coeffs(&[1, 0, 0, -1]);
        let x = Toy::from_i64_coeffs(&[1, -1]);
        let pk1 = &(&alpha * &s) + &x;
        assert_eq!(pk1.coeff_centered(0), -8);

        let y = Toy::from_i64_coeffs(&[-1, 0, 1, 0]);
        let x1 = Toy::from_i64_coeffs(&[1]);
        let x2 = Toy::from_i64_coeffs(&[-1]);
        let t = Toy::from_i64_coeffs(&[8]);
        let a = &(&alpha * &y) + &x1;
        let b = &(&(&pk1 * &y) + &t) + &x2;
        assert_eq!(a.max_centered_abs(), 6);
        assert_eq!(b.coeff_centered(0), -1);

        let d = &b - &(&a * &s);
        assert_eq!(d.coeff_centered(0), 5);
        let noise = &d - &t;
        assert_eq!(
            (0..4).map(|i| noise.coeff_centered(i)).collect::<Vec<_>>(),
            vec![-3, 1, 1, 0]
        );
        assert!(d.coeff_centered(0).abs() > 3);
    }

    #[test]
    fn poly_ops_agree_with_field_ops() {
        let a = Toy::from_i64_coeffs(&[3, 5, -2, 1]);
        let zero = Toy::zero();
        assert_eq!(&a + &zero, a);
        assert_eq!(&a - &zero, a);
        assert_eq!(&a + &(-&a), zero);
        assert_eq!(&zero * &a, zero);
        let two = Fp::<17>::from_i64(2);
        let doubled = a.scalar_mul(two);
        assert_eq!((&a + &a), doubled);
    }

    #[test]
    fn specialized_muls_agree_with_schoolbook() {
        let mut dense = [0u64; 4];
        for (i, c) in dense.iter_mut().enumerate() {
            *c = Fp::<17>::from_i64(5 * i as i64 - 3).raw();
        }
        let a = Toy::from_raw(dense);
        let t = Toy::from_i64_coeffs(&[1, 0, 0, -1]);
        let p = Toy::from_i64_coeffs(&[2, -3, 0, 4]);
        assert_eq!(a.mul_ternary(&t), a.mul(&t));
        assert_eq!(a.mul(&p), a.mul_plain(&p));
        assert_eq!(a.mul_ternary(&t), (&a * &t));
    }

    #[test]
    fn lift_to_preserves_coefficients() {
        let a = Toy::from_i64_coeffs(&[3, 5, -2, 1]);
        let lifted = a.lift_to::<35>();
        assert_eq!(
            (0..4).map(|i| lifted.coeff_centered(i)).collect::<Vec<_>>(),
            (0..4).map(|i| a.coeff_centered(i)).collect::<Vec<_>>()
        );
        let ternary = Toy::from_i64_coeffs(&[-1, 0, 0, 1]);
        let small_mod = ternary.lift_to::<5>();
        assert_eq!(small_mod.coeff_centered(0), -1);
        assert_eq!(small_mod.coeff_centered(3), 1);
        assert_eq!(small_mod.coeff_centered(1), 0);
    }

    #[test]
    fn mul_ternary_rejects_non_ternary_input() {
        let a = Toy::from_i64_coeffs(&[1, 2]);
        let bad = Toy::from_i64_coeffs(&[1, 2]);
        let result = std::panic::catch_unwind(|| a.mul_ternary(&bad));
        assert!(result.is_err());
    }
}
