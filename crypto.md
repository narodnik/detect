---
title: "Sampling RLWEenc Errors for UnifOMR"
subtitle: "A step-by-step walkthrough of `script/sample_errors.sage`, Parameter Set 1"
date: "August 2026"
lang: en
toc: true
numbersections: true
---

# Overview

This document walks through [`script/sample_errors.sage`](script/sample_errors.sage), a
SageMath script that samples the error terms of the **RLWEenc** public-key
encryption scheme used by **UnifOMR** (oblivious message detection/retrieval),
instantiated with **Parameter Set 1** of Table 1, Section 7.1 of the paper.

Why sample errors at all?

* RLWEenc's security rests on the Ring Learning With Errors (RLWE) assumption:
  without the noise $x$, the secret could be recovered from a public key by
  linear algebra over $R_q$ in polynomial time. The Gaussian noise $x$ is what
  makes the problem hard.
* Correctness is a fight *against* that same noise: decryption succeeds only
  if the accumulated noise stays inside the bound $r = 48$. The whole
  parameter set is calibrated so that this holds with probability
  $1 - \epsilon_n$, $\epsilon_n = 2^{-30}$.

So the script does two things: it **builds** the sampler $\chi_\sigma$, and it
**verifies** both the statistical shape of the samples and the paper's
correctness arithmetic.

Everything below references the UnifOMR paper as "the paper" (Definition 4.6
for the scheme, Section 7.1 / Table 1 for the parameters, Theorem 6.1 for the
error analysis).

# The RLWEenc scheme (Definition 4.6)

RLWEenc is a simple PKE over the negacyclic ring $R_q = \mathbb{Z}_q[X]/(X^{n} + 1)$.

**Size of $R_q$.** An element of $R_q$ is a polynomial with $n' = 1024$
coefficients, each an element of $\mathbb{Z}_q$ with $q = 4169729 <
2^{22}$, so each coefficient needs 22 bits. A ring element therefore takes

$$
1024 \times 22 \text{ bits} = 22528 \text{ bits} = \mathbf{2816 \text{ bytes}}
$$

packed optimally (in practice one stores coefficients as 32-bit words,
i.e. 4 KiB per ring element). This is the unit of communication for public
keys and ciphertexts, and it is why the paper works hard to keep $q$ small.

**The distributions.** Three random sources feed the scheme. With the
embedding arrows made explicit:

$$
s \leftarrow D: \quad \{-1,0,1\}^{n} \hookrightarrow R_q, \qquad
s = \sum_{i=0}^{n-1} s_i X^i, \quad \text{exactly } h \text{ of the } s_i \text{ nonzero},
$$

i.e. $D$ is *not* a distribution on all of $R_q$ — it lives on the tiny
subset of ternary, fixed-weight polynomials, embedded into $R_q$ through the
natural map (where $-1$ is represented as $q - 1$). The other two:

$$
\alpha \leftarrow\$ R_q \;\; \text{(uniform over the whole ring)}, \qquad
x, x', x'' \leftarrow \chi_\sigma \;\; \text{(i.i.d. } D_{\mathbb{Z},\sigma} \text{ coefficients, §6)}.
$$

## The moving parts

* **Message.** $\vec m = (m_0, \dots, m_{\ell-1})$ with each
  $m_i \in \{0,1\}$ a single **bit** of the message. The message length is
  $\ell \leq n$: each ciphertext carries $\ell$ bits, one per plaintext
  polynomial coefficient that gets decoded. For Parameter Set 1,
  $\ell = 1$: one RLWEenc ciphertext encrypts exactly **one bit** — enough
  for OMR, where the payload is a detection flag ("this ciphertext is
  pertinent to me: yes/no").
* **Secret key.** $s \leftarrow D$, a *ternary* vector: exactly $h$ of its
  $n$ coefficients are nonzero, each $\pm1$; the rest are 0. The number of
  nonzero entries is the **Hamming weight** — the standard measure of
  sparsity for vectors, defined as
  $\mathrm{wt}(v) = |\{i : v_i \neq 0\}|$. Here $\mathrm{wt}(s) = h = 80$.
* **Noise.** $x, x', x'' \leftarrow \chi_\sigma$: ring elements with i.i.d.
  discrete Gaussian coefficients of standard deviation $\sigma$ (§6).
* **Encryption randomness.** $y \leftarrow$ ternary means $y$ is drawn from
  the *same* kind of distribution as the secret key $s$: a ternary vector of
  Hamming weight $h$. Two different ciphertexts use independent $y$'s.
* **Public randomness.** $\alpha$ uniform in $R_q$.

## The algorithms, and what each value is doing

$$
\begin{aligned}
\textsf{KeyGen}: \quad & \mathrm{pk} = (\underbrace{\alpha}_{\mathrm{pk}[0]},\;
                             \underbrace{\alpha s + x}_{\mathrm{pk}[1]}) \\
\textsf{Enc}(\vec m): \quad & t = \sum_{i=0}^{\ell-1} \frac{q}{2} m_i X^i, \quad
   y \leftarrow \text{ternary}, \\
   & a = \alpha y + x', \quad b = \mathrm{pk}[1]\, y + t + x'' \\
\textsf{Dec}(\mathrm{ct}): \quad & d = b - a s, \quad
   m_i = \mathbb 1[\,|d_i| > r\,]
\end{aligned}
$$

Here $\mathrm{pk}[1] = \alpha s + x$ denotes the **second component** of the
public key (0-indexed), i.e. the "noisy product" part.

Where each value comes from:

| value | source | role |
|---|---|---|
| $\alpha$ | uniform in $R_q$ | public base |
| $s$ | $D$ (ternary, weight $h$) | secret key |
| $x, x', x''$ | $\chi_\sigma$ (per-coefficient $D_{\mathbb{Z},\sigma}$) | noise |
| $y$ | $D$ (ternary, weight $h$) | ephemeral encryption key |
| $m_i$ | $\{0,1\}$ | message bits |
| $\frac{q}{2} m_i$ | $\mathbb{Z}_q$ | encoded bits building $t \in R_q$ |

**The ElGamal-shaped structure.** As a reminder, ElGamal in a cyclic group
$\langle g \rangle$ of order $q'$ is:

$$
\textsf{KeyGen}: \mathrm{sk} = s \leftarrow \mathbb{Z}_{q'}, \;\;
\mathrm{pk} = g^s; \qquad
\textsf{Enc}(m): r \leftarrow \mathbb{Z}_{q'}, \;\;
\mathrm{ct} = (\,g^r,\;\; m \cdot \mathrm{pk}^r\,); \qquad
\textsf{Dec}: m = \mathrm{ct}[1] \,/\, \mathrm{ct}[0]^{s},
$$

where decryption works because
$\mathrm{ct}[0]^s = (g^r)^s = g^{rs} = \mathrm{pk}^r$, so the two DH
"masks" cancel, leaving $m$. RLWEenc is the additive, noisy cousin of this:

$$
g^s \;\leftrightarrow\; \alpha s + x, \qquad
g^r \;\leftrightarrow\; \alpha y + x', \qquad
m\cdot \mathrm{pk}^r \;\leftrightarrow\; \mathrm{pk}[1]\,y + t + x'',
$$

with exponent multiplication $g^{rs}$ replaced by ring multiplication
$\alpha s y = \alpha y s$ (commutativity of $R_q$ plays the role of the
exponent law $rs = sr$), and with a fresh Gaussian added at each step.

* $\alpha$ plays the role of the **public generator** $g$: a public random
  base element (in ElGamal a group generator, here a uniform ring element —
  invertible with overwhelming probability, so it "generates" in spirit).
* $y$ is an **ephemeral secret**, freshly drawn per encryption, exactly like
  ElGamal's per-message $r$: it blinds the ciphertext through the
  Diffie–Hellman-flavored pair of products $\alpha y$ (sent as $a$) and
  $\mathrm{pk}[1]\,y$ (whose $\alpha s y$ part cancels $\alpha y s$ at
  decryption time, leaving the message).
* $x'$ and $x''$ are fresh noise terms that **blind the two published
  values**: $a = \alpha y + x'$ and $b = \mathrm{pk}[1]y + t + x''$. Without
  them, $a, b$ would be exact ring elements and the secret could be solved
  for by linear algebra; with them, the RLWE assumption says $(a, b)$ is
  pseudorandom.
* $x$ in $\mathrm{pk}[1] = \alpha s + x$ plays the same blinding role for
  the public key itself.

## Why the $q/2$ term in $t$

The plaintext bits must live in $\mathbb{Z}_q$, so they need an embedding
into that group. The map

$$
m_i \in \{0, 1\} \;\longmapsto\; \frac{q}{2}\,m_i \in \{0,\; q/2\}
$$

places the two possible bit values at **maximum distance** from each other
on the cyclic group $\mathbb{Z}_q$: going up from 0 you reach $q/2$ after
$q/2$ steps, and going down (mod $q$) also after $q/2$ steps. No other
pair of points in $\mathbb{Z}_q$ is further apart. Decoding then tolerates
noise of up to $\approx q/4$ in either direction before a 0 could be
mistaken for a 1, and the decision itself is trivial:

* $m_i = 0$: $d_i = \text{noise} \in [-r, r]$ — small, decode 0;
* $m_i = 1$: $d_i = q/2 + \text{noise}$ — $\approx q/2 \approx 2.08\times10^6$,
  astronomically far above the threshold $r = 48$, decode 1.

With $q = 4169729$ the gap $q/2$ is about $43{,}000\times$ the error bound —
correctness is never in danger from the encoding; only the *statistical*
budget $\Pr[|\text{noise}| > r] \le \epsilon_n$ matters. (The same encoding
also gives the *wrong-key* behaviour the OMR detector exploits: decrypting
with a mismatched key makes $a s'$ uniform, so $d$ lands at a uniform point
of $\mathbb{Z}_q$, which is $> r$ with probability $1 - 2r/q \approx 1$.)

Our script implements the three samplers (§6–§7) and checks that decryption
noise behaves as the paper's analysis predicts (§8.3).

# Packing: where every value lives

This section goes through the scheme value by value, stating exactly what
each object is, which space it lives in, and where each coefficient comes
from. Everything reduces to one picture: **a ring element is a coefficient
vector wearing polynomial clothes**, and every distribution in the scheme is
just a rule for filling that vector.

## The objects, coefficient by coefficient

All of the following live in $R_q = \mathbb{Z}_q[X]/(X^{1024}+1)$, i.e. each
is a polynomial $\sum_{i=0}^{1023} c_i X^i$ with coefficients reduced
mod $q$ (equivalently: a vector in $\mathbb{Z}_q^{1024}$ whose products wrap
negacyclically):

$$
\begin{aligned}
s &= s_0 + s_1 X + s_2 X^2 + \dots + s_{1023} X^{1023},
  && s_i \leftarrow \begin{cases} \pm1 & \text{for exactly } h=80 \text{ positions} \\ 0 & \text{otherwise} \end{cases} \\[4pt]
\alpha &= \alpha_0 + \alpha_1 X + \dots + \alpha_{1023} X^{1023},
  && \alpha_i \leftarrow\$ \mathbb{Z}_q \text{ (uniform)} \\[4pt]
x = x' = x'' &= \dots \text{ three independent draws of } \dots
  && x_i \leftarrow D_{\mathbb{Z},\,0.6} \text{ (i.i.d., §6)} \\[4pt]
y &= y_0 + y_1 X + \dots + y_{1023} X^{1023},
  && y_i \leftarrow \begin{cases} \pm1 & \text{for exactly } h=80 \text{ positions} \\ 0 & \text{otherwise} \end{cases} \\[4pt]
t &= \tfrac{q}{2} m_0 + \tfrac{q}{2} m_1 X + \dots + \tfrac{q}{2} m_{\ell-1} X^{\ell-1},
  && \tfrac{q}{2} m_i \in \{0,\; q/2\} \subset \mathbb{Z}_q,\; m_i \in \{0,1\}
\end{aligned}
$$

(The message $t$ is zero-padded from degree $\ell - 1 = 0$ up to the full
ring dimension.)

## A worked KeyGen, Enc, and Dec, in miniature

To see the packing *and* the negacyclic wrap in action, take a toy ring
$R_q = \mathbb{Z}_{17}[X]/(X^4+1)$ ($n' = 4$, $q = 17$, toy $r = 3$,
$\ell = 1$ so $\frac{q}{2}m_0 = 8m_0$). Draw

$$
\alpha = 3 + 5X - 2X^2 + X^3 \leftarrow\$ R_q, \qquad
s = 1 - X^3 \leftarrow D \;(h=2), \qquad
x = 1 - X \leftarrow \chi_\sigma .
$$

**KeyGen.** Multiply out $\alpha \cdot s$ before reduction:

$$
\alpha s = 3 + 5X - 2X^2 + X^3 - 3X^3 - 5X^4 + 2X^5 - X^6 .
$$

Now reduce with $X^4 = -1$, i.e. $X^4 \mapsto -1$, $X^5 \mapsto -X$,
$X^6 \mapsto -X^2$ — this is the negacyclic wrap, where the wrapped terms
flip sign:

$$
\alpha s = 3 + 5X - 2X^2 + X^3 - 3X^3 + 5 - 2X + X^2
         = 8 + 3X - X^2 - 2X^3 .
$$

Note how $s$ having weight 2 made only two "copies" of $\alpha$ appear,
with one of them sign-flipped — this is exactly the mechanism behind the
variance count $h + h + 1$ in §8.3.1. Finally add the noise:

$$
\mathrm{pk}[1] = \alpha s + x = (8 + 3X - X^2 - 2X^3) + (1 - X)
              = 9 + 2X - X^2 - 2X^3 \;\in R_q .
$$

**Enc of the bit $m_0 = 1$.** Draw the ephemeral key
$y = -1 + X^2 \leftarrow D$ and fresh noise $x' = 1$, $x'' = -1$, and
encode the message as $t = 8m_0 = 8$. Compute the ciphertext
$(a, b) = (\alpha y + x',\; \mathrm{pk}[1]\, y + t + x'')$:

$$
\alpha y = -3 - 5X + 5X^2 + 4X^3 - 2X^4 + X^5
         \overset{X^4=-1}{=} -1 - 6X + 5X^2 + 4X^3,
$$

$$
\mathrm{pk}[1]\, y = -9 - 2X + 10X^2 + 4X^3 - X^4 - 2X^5
                   = -8 + 10X^2 + 4X^3,
$$

so

$$
a = \alpha y + x' = -6X + 5X^2 + 4X^3, \qquad
b = \mathrm{pk}[1]y + t + x'' = -1 + 10X^2 + 4X^3 .
$$

**Dec.** Compute $a\,s$ (another negacyclic product — wrap terms at
$X^4, X^5, X^6$ again flip sign):

$$
a s = -6X + 5X^2 + 4X^3 + 6X^4 - 5X^5 - 4X^6
    = -6 - X + 9X^2 + 4X^3,
$$

and subtract:

$$
d = b - a s = (-1 + 10X^2 + 4X^3) - (-6 - X + 9X^2 + 4X^3)
  = 5 + X + X^2 .
$$

**Check what happened.** All the $\alpha$-terms cancelled exactly as
designed ($\alpha s y = \alpha y s$), leaving

$$
d = t + \underbrace{(x y + x'' - x' s)}_{\text{noise}}
  = 8 + (-3 + X + X^2) = 5 + X + X^2 .
$$

The per-coefficient noise is $(-3, 1, 1, 0)$ — each entry a sum of
$\pm1$-weighted Gaussian draws, exactly the $\sqrt{2h+1}\,\sigma$
accumulation of §8.3.1 (here $2h+1 = 5$ toy terms). Decoding the single
payload coefficient: $|d_0| = 5 > r = 3$, so the recipient outputs
$m_0 = 1$ — correct. Had $m_0 = 0$, then $d_0$ would be the pure noise
$-3$, and $|-3| \le r$ would decode to $0$.

**One more view — the linear form.** The constant coefficient of $a s$ is
$a_1 s_3$ with a wrap: the product $a_1 X \cdot s_3 X^3$ produces
$a_1 s_3 X^4 = -a_1 s_3$, so

$$
(a s)_0 = -a_1\, s_3 = -(-6)(-1) = -6 .
$$

This is precisely the "public scalar $\times$ secret coefficient"
shape $A[0][i]\, s_i$ that the detector's BFV circuit evaluates
homomorphically (next subsection) — the negacyclic ring multiplication,
unfolded.

In the real scheme this arithmetic happens 1024 coefficients wide, mod
$q = 4169729$ — same structure, bigger numbers.

## Packing into BFV: coefficients, slots, and constants

The detector side of UnifOMR evaluates the *decryption* of many clues at
once under BFV homomorphic encryption, and the way the values are packed
into BFV plaintexts is the clever part of the construction. The BFV
plaintext ring is $\mathbb{Z}_t[X]/(X^{D}+1)$ with $D = 2048$ and
$t \equiv 1 \!\! \pmod{2D}$, which (by CRT) decomposes into $D$
independent "$\mathbb{Z}_t$-slots" — one SIMD lane each (§4.5 of the paper).

Step by step, for the first decoded bit ($j = 0$) of $D$ different clues:

1. **Break the ring multiply into scalars.** Decryption needs
   $b[0] - (a\cdot s)_0$, and the constant coefficient of the negacyclic
   product is a *linear form* in the secret coefficients:
   $(a s)_0 = \sum_{i=0}^{n-1} A[0][i]\, s_i$, where the public scalars
   $A[0][i] \in \mathbb{Z}_t$ collect the coefficients of $a$ together with
   their wrap-around signs. The ring structure of $R_q$ has been unfolded
   into a plain matrix–vector product $A \cdot \vec s$ over $\mathbb{Z}_t$.

2. **Encrypt the secret as constants.** Each secret coefficient is
   encrypted *by itself* as a constant polynomial
   $s_i(X) := s_i \in \mathbb{Z}_t[X]/(X^D+1)$, degree 0:
   $\mathrm{ct}_{s_i} \leftarrow \textsf{BFV.Enc}(s_i(X))$.

3. **Pack $D$ clues into coefficients.** For the $i$-th secret coefficient,
   build $p_i(X) = \sum_{j=0}^{D-1} A_j[0][i]\, X^j$, where $A_j$ is the
   matrix of the $j$-th of the $D$ packed clues — i.e. the $D$ SIMD scalars
   are placed in the $D$ *coefficients* of one BFV plaintext.

4. **Multiply.** Because $\mathrm{ct}_{s_i}$ encrypts a constant,
   multiplying $p_i(X) \cdot \mathrm{ct}_{s_i}$ is just *scalar-by-vector*:
   coefficient $j$ of the result gets $A_j[0][i] \cdot s_i$. Summing over
   $i \in [n]$ yields coefficient $j = \sum_i A_j[0][i] s_i = (A_j \vec s)_0$
   — the partial decryption of clue $j$, all $D$ clues in one homomorphic
   pass. The final range check $|\cdot| \le r$ is deferred to the recipient
   in the clear, avoiding a deep circuit.

This is why the paper can set $q = t$ (footnote 6 of §5): the RLWEenc
modulus doubles as the BFV plaintext modulus, so the partial decryption
arithmetic lands directly in $\mathbb{Z}_t$ with no rescaling until the
deliberate modulus switch $Q \to Q'$.

## A worked packing, in miniature (four clues)

To see the packing itself in action, continue the toy world of §3.2 with a
toy BFV ring $\mathbb{Z}_{17}[X]/(X^4+1)$: dimension $D = 4$, so **four
"slots" = four clues processed at once**. The recipient is the same one as
before, so the secret is still $s = 1 - X^3$, i.e.
$\vec s = (1, 0, 0, -1)$ (weight 2 at positions 0 and 3).

Four senders produce four clues, each with its own public $a$-part:

$$
\begin{aligned}
a^{(0)} &= 1 + 2X + 5X^3, & a^{(1)} &= 4 - 3X + 2X^2, \\
a^{(2)} &= -2 + X - 4X^3, & a^{(3)} &= 6 - X + 3X^2 + 2X^3 .
\end{aligned}
$$

**Step 1 — unfold to scalars.** The constant coefficient of a negacyclic
product pairs position $i$ with $(0-i) \bmod 4$, negating when it wraps
(§3.2's linear form):

$$
(a s)_0 = a_0 s_0 \;-\; a_3 s_1 \;-\; a_2 s_2 \;-\; a_1 s_3
        \;=\; \sum_{i=0}^{3} A[0][i]\, s_i .
$$

**Step 2 — encrypt the secret as constants.** The detector holds
$\mathrm{ct}_{s_i} = \textsf{BFV.Enc}(\text{constant } s_i)$ for each $i$;
on plaintexts these are simply

$$
\mathrm{ct}_{s_0} \mapsto 1, \quad \mathrm{ct}_{s_1} \mapsto 0, \quad
\mathrm{ct}_{s_2} \mapsto 0, \quad \mathrm{ct}_{s_3} \mapsto -1 .
$$

(In the real construction each is a full BFV ciphertext, carrying BFV noise
that never mixes with the RLWEenc noise until the merging step.)

**Step 3 — pack the four clues into coefficients.** For each secret
position $i$, clue $j$'s scalar $A_j[0][i]$ goes into *coefficient $j$* of
one plaintext polynomial:

$$
\begin{aligned}
p_0(X) &= 1 + 4X - 2X^2 + 6X^3 &&\text{(the } a^{(j)}_0 \text{)}\\
p_1(X) &= -5 + 0X + 4X^2 - 2X^3 &&\text{(the } -a^{(j)}_3 \text{)}\\
p_2(X) &= 0 - 2X + 0X^2 - 3X^3 &&\text{(the } -a^{(j)}_2 \text{)}\\
p_3(X) &= -2 + 3X - X^2 + X^3 &&\text{(the } -a^{(j)}_1 \text{)}
\end{aligned}
$$

For example $p_1$ reads down the $a_3$ column — clue 0 has $a_3 = 5$, so
coefficient 0 gets $-5$; clue 2 has $a_3 = -4$, so coefficient 2 gets $+4$.

**Step 4 — homomorphically multiply and sum.** The detector computes
$\sum_{i} p_i(X) \cdot \mathrm{ct}_{s_i}$. Since each $\mathrm{ct}_{s_i}$ is
a constant, this is four scalar-broadcasts; on plaintexts:

$$
1\cdot p_0 + 0\cdot p_1 + 0\cdot p_2 + (-1)\cdot p_3 = p_0 - p_3
= 3 + X - X^2 + 5X^3 .
$$

**Read off the answers.** Coefficient $j$ of the result is clue $j$'s
partial decryption $(a^{(j)} s)_0$ — all four at once:

| coefficient | value | direct check $(a^{(j)}s)_0$ |
|---|---|---|
| 0 | $3$ | $1 + 2\cdot(-s_3) = 1 + 2 = 3$ $\checkmark$ |
| 1 | $1$ | $4 + (-3) = 1$ $\checkmark$ |
| 2 | $-1$ | $-2 + 1 = -1$ $\checkmark$ |
| 3 | $5$ | $6 + (-1) = 5$ $\checkmark$ |

The recipient then finishes each clue in the clear: for message $j$, take
$b^{(j)}[0]$ minus coefficient $j$, and range-check $|\cdot| \le r$ to
recover the bit (§3.2's decoding).

The cost accounting is the whole point: $n = 4$ plaintext-ciphertext
multiplications handled $D = 4$ clues. In the real scheme, $n = 900$
such multiplications process $D = 2048$ clues in one SIMD pass — and for
$\ell > 1$ the whole procedure simply repeats per decoded bit $j$.

> **Why send all of $a$ when only the constant coefficient of $a\cdot s$
> is kept?** Three reasons. *(1) The whole of $a$ is used.* The kept
> coefficient is a linear form in *every* coefficient of $a$ — in the real
> scheme, $(as)_0 = \sum_{i=0}^{1023} \pm a_i\, s_{(0-i) \bmod n'}$ — so
> the $\ell$ decoded bits are compressions of $a$, not selections from it.
> *(2) $a$ is the armored form of the ephemeral key, and armor doesn't
> compress.* $a = \alpha y + x'$ is the per-message secret $y$ wrapped in
> RLWE; if $y$ leaked, anyone could compute $\mathrm{pk}[1]\,y$, strip $b$
> down to $t + x'' \approx t$, and read the message. So the sender cannot
> ship $y$ (~58 bytes) or a seed for it — only its full-dimensional,
> uniform-looking image, which is exactly what RLWE security requires. The
> 2816 bytes of $a$ are the ciphertext's entropy, not waste. *(3) The
> scheme already trims what can be trimmed.* The ciphertext type is
> $(a, \vec b) \in R_q \times \mathbb{Z}_q^{\ell}$: $b$ is truncated to its
> $\ell$ needed coefficients (~3 bytes for Param1), which costs no
> security, and the detector's reply is compressed to $\ell$ elements of
> $\mathbb{Z}_{Q'}$ per message by the merging step — that is where the
> 20×–1080× communication win of UnifOMR comes from, not from shrinking
> $a$.

## Mechanizing it: `script/packing.sage`

Both worked examples above are mechanized in
[`script/packing.sage`](script/packing.sage), which replays them with assertions against
the doc's numbers and then re-verifies everything on fresh random draws.
This subsection walks through the script, its output, and what each piece
certifies.

### Setup and centered rendering

```python
q = 17
nring = 4
P = PolynomialRing(Zmod(q), "X")
X = P.gen()
R = P.quotient(X^nring + 1, "x")

def centered_coeffs(f):
    cs = f.lift().coefficients(sparse=False)
    cs += [0] * (nring - len(cs))
    return [ZZ(c) - q if ZZ(c) > q // 2 else ZZ(c) for c in cs]
```

Sage stores coefficients as residues $0..16$, but the scheme's numbers are
naturally small signed integers, so every readout passes through
`centered_coeffs`, which lifts each coefficient into $(-q/2, q/2]$ — this
is also exactly the decoder's point of view, since $|d_i| \le r$ is only
meaningful for centered values. One consequence to keep in mind when
comparing outputs with §3.2's by-hand arithmetic: $9 \equiv -8$ and
$10 \equiv -7 \pmod{17}$, so the script prints
$\mathrm{pk}[1] = -8 + 2X - X^2 - 2X^3$ where the hand computation wrote
$9 + 2X - \dots$ — same element, different representative.

### Replaying §3.2: KeyGen / Enc / Dec

The script pins the exact draws of §3.2 and asserts every intermediate:

```python
pk1 = alpha * s + x
assert centered_coeffs(pk1) == [-8, 2, -1, -2]

a = alpha * y + x1
b = pk1 * y + t + x2
assert centered_coeffs(a) == [0, -6, 5, 4]
assert centered_coeffs(b) == [-1, 0, -7, 4]

d = b - a * s
assert centered_coeffs(d) == [5, 1, 1, 0]
```

Run output (verbatim):

```
=== worked KeyGen / Enc / Dec (section 3.2) ===
alpha                : 3 + 5X - 2X^2 + X^3
s                    : 1 - X^3   (weight 2)
pk[1] = alpha*s + x  : -8 + 2X - X^2 - 2X^3
a = alpha*y + x'     : -6X + 5X^2 + 4X^3
b = pk[1]*y + t + x''  : -1 - 7X^2 + 4X^3
d = b - a*s          : 5 + X + X^2
noise = d - t        : -3 + X + X^2   per-coeff [-3, 1, 1, 0]
decode |d_0|=5 vs r=3  : m_0 = 1   (sent 1)
```

Every line lands on the by-hand values of §3.2 (modulo the
representative choice above), including the punchlines: the noise
per-coefficient $(-3, 1, 1, 0)$ — the toy instance of the
$\sqrt{2h+1}\,\sigma$ accumulation — and the decode
$|d_0| = 5 > r = 3 \Rightarrow m_0 = 1$.

### Replaying §3.4: the unfolding, in code

The linear-form unfolding is one function, implementing
$A[k][i] = \pm a_{(k-i) \bmod n'}$ with the sign negative exactly when the
product wraps past $X^{n'} = -1$:

```python
def lin_form_scalars(a_coeffs, k=0):
    scalars = []
    for i in range(nring):
        val = a_coeffs[(k - i) % nring]
        scalars.append(ZZ(val) if i <= k else -ZZ(val))
    return scalars
```

(For the general decoded coefficient $k$: position $i \le k$ pairs with
$a_{k-i}$ unwrapped; $i > k$ pairs with $a_{k-i+n'}$ and carries the
negacyclic minus sign. The toy decodes $k = 0$ only.)

Then packing is a transpose — secret position $i$ indexes the polynomial,
clue $j$ indexes the coefficient:

```python
A = [lin_form_scalars(clues[j]) for j in range(nring)]
p = [R([A[j][i] for j in range(nring)]) for i in range(nring)]
total = sum(p[i] * R(s_vec[i]) for i in range(nring))
```

Run output (verbatim):

```
=== worked four-clue packing (section 3.4) ===
clue 0: a^(0) = 1 + 2X + 5X^3
clue 1: a^(1) = 4 - 3X + 2X^2
clue 2: a^(2) = -2 + X - 4X^3
clue 3: a^(3) = 6 - X + 3X^2 + 2X^3
p_0(X) = 1 + 4X - 2X^2 + 6X^3
p_1(X) = -5 + 4X^2 - 2X^3
p_2(X) = -2X - 3X^3
p_3(X) = -2 + 3X - X^2 + X^3
sum_i p_i * s_i      : 3 + X - X^2 + 5X^3
coefficient 0 = +3   direct (a^(0)*s)_0 = +3   [ok]
coefficient 1 = +1   direct (a^(1)*s)_0 = +1   [ok]
coefficient 2 = -1   direct (a^(2)*s)_0 = -1   [ok]
coefficient 3 = +5   direct (a^(3)*s)_0 = +5   [ok]
```

The four `[ok]` lines are the script cross-checking the SIMD result against
brute force: for each clue $j$ it also computes the honest ring product
$a^{(j)} \cdot s$ in $R$ and compares its constant coefficient with
coefficient $j$ of $\sum_i p_i s_i$. Packing provably loses nothing.

### Randomized checks

Finally the script abandons the fixed draws and loops:

* **200 full Enc/Dec cycles** with fresh $\alpha, s, y$ and Gaussian noise
  (the same `GeneralDiscreteDistribution`-based $\chi_\sigma$ as
  `script/sample_errors.sage`), asserting the decoded bit equals the sent bit;
* **50 packing runs** with random secrets and random clue $a$'s, asserting
  the coefficient-extraction identity holds for every clue.

```
=== randomized checks ===
random trials        : 199/200 decodes correct
random packing trials: 50/50 coefficient-extraction checks passed
```

The packing check is exact algebra — it must always pass, and any failure
would mean a bug in the unfolding. The decode trials are statistical: with
toy noise std $\sqrt{2h+1}\,\sigma = \sqrt5 \cdot 0.6 \approx 1.34$ against
$r = 3$, the tail $\Pr[|\text{noise}| > r] \approx 2\%$ — the toy's own
false-negative budget, i.e. the 1-in-200 miss seen above is the miniature
version of exactly the $\mathrm{erfc}$ analysis of §8.3 (where the real
parameters push the tail to $2.9\times10^{-10}$).

## Communication, after merging

The raw partial decryption would be a BFV ciphertext per clue, but after
the coefficient-packing trick and modulus switching, each message costs
only $\ell$ elements of $\mathbb{Z}_{Q'}$ (a few bytes each) — versus
hundreds of bytes for shipping the whole RLWEenc ciphertext to the
recipient. This factor is the "20×–1080×" headline improvement of UnifOMR.

# Parameter Set 1

From Table 1, Section 7.1:

| symbol | value | meaning |
|---|---|---|
| $n$ | 900 | RLWEenc secret key dimension |
| $n'$ | 1024 | ring dimension ($n$ padded up, footnote 10) |
| $q$ | 4169729 | ciphertext modulus |
| $\sigma$ | 0.6 | RLWEenc error standard deviation |
| $h$ | 80 | secret key Hamming weight |
| $r$ | 48 | error range bound |
| $\ell$ | 1 | bits encrypted per ciphertext |
| $\epsilon_n$ | $2^{-30}$ | false negative rate |
| $\epsilon_p$ | $2^{-15}$ | false positive rate |

(BFV-side parameters $D, Q, h_{\text{BFV}}, Q', \sigma_{\text{BFV}}, r'$ are
discussed in §10; the sampler does not touch them.)

# Step 1 — the cyclotomic (negacyclic) ring

## Why a power-of-two dimension

The secret dimension is $n = 900$, not a power of two. Footnote 10 of the
paper: $n$ is zero-padded up to the smallest power of two $n' \geq n$, here
$1024$, so that the ring is a power-of-two cyclotomic ring.

The reason for insisting on a power of two is **fast multiplication via the
convolution theorem**: multiplying two polynomials is a convolution of their
coefficient vectors, and the DFT/FFT (over $\mathbb{C}$) or the NTT (the
number-theoretic transform, its exact finite-field analogue over
$\mathbb{Z}_q$) turns convolution into pointwise multiplication —

$$
a \cdot b \;=\; \mathrm{NTT}^{-1}\big(\mathrm{NTT}(a) \odot \mathrm{NTT}(b)\big),
$$

reducing ring multiplication from $O(n^2)$ schoolbook work to
$O(n \log n)$. Radix-2 FFT/NTT algorithms need the transform length to be a
power of two (here $2n' = 2048$ for the negacyclic convolution, §5.3) to get
their divide-and-conquer splitting; a dimension like 900 would force either
a generic (slower) transform or zero-padding anyway. Since RLWEenc
encryption, decryption *and* the paper's homomorphic circuit are dominated
by ring multiplications, this is not a minor optimization — it is the
difference between practical and impractical.

```python
n = 900
nring = 1024
m = 2 * nring          # 2048
q = 4169729
```

## The ring in Sage

```python
P = PolynomialRing(Zmod(q), "X")
X = P.gen()
phi = P(cyclotomic_polynomial(m))
R = P.quotient(phi, "x")
```

The 2048-th cyclotomic polynomial satisfies
$\Phi_{2048}(X) = X^{1024} + 1$ (the script asserts this at runtime), so

$$R = \mathbb{Z}_q[X]/(X^{1024}+1).$$

## What "negacyclic" means

The defining relation is $X^{n'} = -1$: powers of $X$ wrap around the top with
a **sign flip**, $X^{n'+k} = -X^k$. Polynomial multiplication therefore becomes
a cyclic convolution in which the wrapped terms are *negated*:

$$
\Big(\sum_i a_i X^i\Big)\Big(\sum_j b_j X^j\Big)
\;\Rightarrow\;
c_k \;=\; \sum_{i+j=k} a_i b_j \;-\!\!\sum_{i+j=k+n'} a_i b_j
\pmod q .
$$

Compare *cyclic* reduction modulo $X^{n'}-1$ (i.e. $X^{n'} = 1$), where the
wrapped terms are added instead. Negacyclic rings are the standard choice in
lattice crypto because $\Phi_{2n'}$ is the sparse polynomial $X^{n'}+1$, and
multiplication mod $X^{n'}+1$ maps to an NTT (negacyclic convolution) with no
twist overhead.

The demo in §8.4 exercises exactly this: it multiplies in $R$ and checks the
wraparound/sign behaviour by recovering the error term exactly.

# Step 2 — the error distribution $\chi_\sigma$

## Definition

**What it is.** A *discrete Gaussian* is the continuous bell curve
restricted to the integers: take the Gaussian density evaluated at each
integer, and renormalize so the values sum to 1. The standard notation in
the lattice-crypto literature is $D_{\mathbb{Z}, \sigma}$:

$$
D_{\mathbb{Z},\sigma}(k) \;=\; \frac{\rho_\sigma(k)}{\rho_\sigma(\mathbb{Z})}
\;=\; \frac{\exp(-k^2/2\sigma^2)}{\displaystyle\sum_{j \in \mathbb{Z}} \exp(-j^2/2\sigma^2)},
\qquad k \in \mathbb{Z},
$$

where the symbols mean:

* $\mathbb{Z} = \{\dots, -2, -1, 0, 1, 2, \dots\}$ — the integers, the
  support of the distribution (each coefficient of the error *is an
  integer*, because ring elements have integer coefficients);
* $k$ — the integer being sampled;
* $j$ — a dummy index running over all of $\mathbb{Z}$ in the normalizing
  sum;
* $\rho_\sigma(x) = e^{-x^2/2\sigma^2}$ — the (unnormalized) Gaussian mass
  function, and $\rho_\sigma(\mathbb{Z}) = \sum_j \rho_\sigma(j)$ its total
  mass over the integers, which replaces the continuous Gaussian's
  $\sqrt{2\pi}\,\sigma$ normalization constant;
* $\sigma$ — the standard deviation: with the $2\sigma^2$ in the exponent,
  the distribution really has spread $\sigma$.

(Some references, e.g. the original RLWE paper of Lyubashevsky–Peikert–Regev,
use the "$\pi$-convention" $D_{\mathbb{Z},s}(k) \propto e^{-\pi k^2/s^2}$
with parameter $s$; the two are the same distribution with
$s = \sigma\sqrt{2\pi}$. The UnifOMR paper's $\sigma = 0.6$ is in the
standard-deviation convention used above.)

**Why such a thing exists.** The error must be *discrete* — ciphertext
coefficients live in $\mathbb{Z}_q$, so an error term has to be an integer
vector to be addable to them. But it must also be *Gaussian-shaped*, for
two reasons:

1. **Security.** The RLWE hardness results (worst-case to average-case
   reductions from approximate lattice problems) are proven for (discrete)
   Gaussian errors. Spherical-ish Gaussian noise is also the "least
   structured" noise, giving the best concrete security per unit of
   variance — attacks rely on errors being small, and for a fixed variance
   the Gaussian is the shape that maximizes entropy.
2. **Correctness calculus.** Sums and products of many small Gaussians
   stay statistically predictable (§8.3): variance adds linearly, which is
   what lets the paper write down a closed-form failure probability.

**What it is used for.** Every "small" random quantity in LWE-style
cryptography: the errors $x, x', x''$ here, RLWE secrets, lattice-trapdoor
samplers, FHE noise. In this script, one draw $k \leftarrow D_{\mathbb{Z},0.6}$
per coefficient builds an error polynomial $x \leftarrow \chi_\sigma$.

An error polynomial $x \leftarrow \chi_\sigma$ has $n' = 1024$ i.i.d.
coefficients, each from $D_{\mathbb{Z}, \sigma}$ with $\sigma = 0.6$.

## Numbers, and the shape

With $\sigma = 0.6$:

$$Z = 1.5065,\qquad
\Pr[0] = \tfrac{1}{Z} = 0.6638,\quad
\Pr[\pm1] = \tfrac{e^{-1/0.72}}{Z} = 0.1655 \text{ each},\quad
\Pr[\pm2] = 0.00257 \text{ each}.$$

So $\chi_\sigma$ is **nearly ternary**: about 99.5% of coefficients land in
$\{0, \pm1\}$. The exact variance is
$\sum_k k^2\Pr[k] = 0.3516$, slightly below $\sigma^2 = 0.36$ — this is why
the script's empirical std reads ~0.593 rather than 0.600 (§8.2).

The distribution, to scale (one `#` $\approx$ probability 0.01):

```
  k        weight w(k)      P(k)                   shape
                            (w(k)/Z)
  -3       3.7e-6           0.0000025
  -2       3.9e-3           0.00257            ##
  -1       0.2494           0.16552            ################
   0       1.0000           0.66382            ##################################################################
   1       0.2494           0.16552            ################
   2       3.9e-3           0.00257            ##
   3       3.7e-6           0.0000025

             normalization  Z = 1.50644
             exact variance 0.35162   (vs sigma^2 = 0.36)
```

## Implementation, worked out

The weights are computed as $w(k) = e^{-k^2/2\sigma^2} = e^{-k^2/0.72}$,
step by step:

| $k$ | exponent $-k^2/0.72$ | $w(k) = e^{(\cdot)}$ | $P(k) = w(k)/Z$ |
|---|---|---|---|
| 0 | $-0$ | $1.0$ | $0.66382$ |
| $\pm1$ | $-1.3889$ | $0.24935$ | $0.16552$ |
| $\pm2$ | $-5.5556$ | $0.003866$ | $0.00257$ |
| $\pm3$ | $-12.5$ | $3.73\times10^{-6}$ | $2.5\times10^{-6}$ |
| $\pm8$ | $-88.9$ | $2.8\times10^{-39}$ | — (cutoff) |

```python
PREC = RealField(200)
support = list(range(-ceil(12 * sigma), ceil(12 * sigma) + 1))
weights = [PREC(exp(-x^2 / (2 * sigma^2))) for x in support]
chi_sigma = GeneralDiscreteDistribution(weights)

def chi_coeff():
    return ZZ(support[chi_sigma.get_random_element()])
```

Remarks:

* The support is cut off at $\pm\lceil 12\sigma\rceil = \pm 8$; the mass
  beyond is $e^{-72} \sim 10^{-31}$, unreachable in a lifetime of sampling.
* `GeneralDiscreteDistribution` takes **unnormalized** weights — it
  normalizes internally by the sum $\approx 1.5065$, so the raw weights
  summing to 1.5 is harmless.
* Sage quirk: its `get_random_element()` returns the **index** into the
  weight list, not the value; `chi_coeff` maps the index back through
  `support`.
* An error polynomial is then

```python
def sample_error():
    return R([chi_coeff() for _ in range(nring)])
```

**A genuine sample.** One run of `sample_error()` drew 1024 coefficients
with counts $\{-2{:}\,2,\; -1{:}\,165,\; 0{:}\,692,\; +1{:}\,162,\; +2{:}\,3\}$
— 332 nonzero coefficients out of 1024, matching the probabilities above.
The resulting ring element is the sparse polynomial

$$
x = -X^2 - X^3 - X^4 + X^5 + X^7 - X^{13} + X^{15} - X^{18} + X^{20}
    - X^{21} + X^{26} - \dots - X^{1023} \;\in\; R,
$$

i.e. about a third of the coefficients are $\pm1$, a handful are $\pm2$,
and the rest vanish — a concrete picture of "near-ternary noise".

# Step 3 — the secret distribution and $\alpha$

The secret is ternary with fixed Hamming weight $h = 80$ over the first
$n = 900$ positions (zero-padded to the ring dimension), and $\alpha$ is
uniform:

```python
def sample_secret():
    coeffs = [0] * nring
    for i in Subsets(range(n), h).random_element():
        coeffs[i] = choice([-1, 1])
    return R(coeffs)

def sample_uniform():
    return R([Zmod(q).random_element() for _ in range(nring)])
```

The sparse ternary secret is an efficiency choice (fast multiplication, small
keys); hardness of sparse-secret RLWE is standard, and the paper certifies the
concrete parameters with the **LWE estimator** at 110-bit security.

## Why $\pm1$ coefficients (ternary), and why weight $h = 80$

**Significance of the $\pm1$ entries.** The secret being ternary — every
nonzero coefficient exactly $\pm1$ — pays off four separate ways:

1. **Multiplication becomes addition.** Computing $a \cdot s$ needs no
   integer multiplications at all: for each nonzero $s_i = \pm 1$, add or
   subtract one rotated copy of $a$ (the negacyclic rotation happens for
   free inside the ring). With $h = 80$ nonzeros out of $1024$ positions,
   that is 80 vector additions instead of a $1024 \times 1024$
   multiply-accumulate — a $\sim\!12\times$ saving on the exact operation
   decryption performs, mirrored in the detector's homomorphic version.
2. **Bounded noise growth.** The accumulated decryption noise is a sum of
   $2h + 1$ terms of the form $x_i \cdot s_j$ (§8.3.1). Because
   $s_j = \pm 1$, each term is exactly as large as $x_i$ itself, giving
   effective std $\sqrt{2h+1}\,\sigma$. Had the secret been, say, uniform in
   $\mathbb{Z}_q$, each term would carry a factor $\sim q/4$, the noise
   would be $\sim\!10^5$ times wider, and no feasible $r$ could rescue
   correctness. Small secret entries are what make the whole
   noise-calculus of §8 work.
3. **Small keys.** The key is fully described by the $h$ positions and their
   signs: $\log_2\!\big(\binom{900}{80} 2^{80}\big) = 465$ bits
   $\approx 58$ bytes, versus 2816 bytes for a dense ring element (§2).
4. **Sufficient entropy all the same.** 465 bits of key-space entropy
   (min-entropy of the distribution) is far beyond the 110-bit security
   target, so guessing is hopeless; the *real* threat for sparse secrets is
   specialized lattice/combinatorial attacks, which is exactly what the
   LWE estimator certifies for this shape.

**Why $h = 80$ specifically.** The weight sits on a two-sided trade-off:

* *Correctness/communication pushes $h$ down.* The noise std grows like
  $\sqrt{2h+1}$, and everything downstream — $r$, then $r' = (Q'/t)r +
  \mathrm{eModSW} + \lceil Q'/t \rceil$, then the false-positive bound
  $\epsilon_p$ — scales with it. Concretely, holding all else fixed:

  | $h$ | noise std $\sqrt{2h+1}\,\sigma$ | min $r$ | resulting $\epsilon_p$ |
  |---|---|---|---|
  | 80 | 7.61 | 47 | $2^{-15.05}$ $\checkmark$ |
  | 400 | 16.98 | 104 | (worse) |
  | 900 (dense) | 25.46 | 156 | $2^{-13.6}$ $\times$ |

  A fully dense ternary secret would blow the $\epsilon_p = 2^{-15}$ budget
  by $\sim\!3\times$, forcing a larger $t = q$ or smaller $Q'$ — i.e. more
  communication, the metric UnifOMR exists to minimize.
* *Security pushes $h$ up.* The sparser the secret, the more leverage the
  known sparse-secret attacks get; $h$ must stay comfortably above the
  range where hybrid guessing or BKW-style combinatorics beat lattice
  attacks. $h = 80$ with $n = 900$ ($\sim\!9\%$ density) keeps the
  estimator's verdict at 110 bits.
* *Matching the BFV side.* The merged-error analysis (§10) needs
  $\mathrm{eModSW} \le h_{\text{BFV}}/2$ to stay comparable to the RLWEenc
  error $\Theta(h)$; the paper sets both weights linear in the security
  parameter — $h = 80$ against $h_{\text{BFV}} = 400$ — so neither noise
  source dominates after merging.

So $h = 80$ is the calibrated sweet spot: enough mass for 110-bit concrete
security, little enough that $r = 48$ fits the false-positive budget, and
formally aligned with the BFV secret's weight.

> **What is an LWE estimator?** The *lwe-estimator* (Albrecht, Player,
> Scott) is the standard open-source tool for *concrete* security of
> (Ring-)LWE instances. You feed it the instance parameters — dimension,
> modulus $q$, error width $\sigma$, secret distribution, number of samples
> available — and it estimates the computational cost of the cheapest known
> attack. The attack families it models are lattice-reduction based: the
> "primal" (uSVP embedding) and "dual" (BDD/short dual vector) attacks,
> which reduce to running BKZ with some block size $\beta$ on a related
> lattice, plus combinatorial/code-based attacks (BKW, decoding, arithmetic
> progressions). Its output is a work factor; "110-bit security" means the
> best attack it finds costs $\approx 2^{110}$ operations. It is an
> estimate *against known attacks* — better algorithms would lower it — but
> it is the community's yardstick, and the number quoted in Table 1.

# Step 4 — verification, block by block

## Group 1: the ring sanity check

```
ring                  : Univariate Quotient Polynomial Ring in x over Ring of integers modulo 4169729 with modulus X^1024 + 1
defining relation     : Phi_2048(X) = True (X^1024 + 1)
```

Confirms $\Phi_{2048}$ collapses to $X^{1024}+1$, i.e. we really are in the
negacyclic ring of Definition 4.6.

## Group 2: statistics of $\chi_\sigma$

2000 error polynomials $\times$ 1024 coefficients = 2,048,000 draws:

```
empirical mean        : -0.00122      (centered, as expected)
empirical std         : 0.59351       (target 0.6)
empirical P(0)        : 0.66351       (theory 0.6638)
empirical P(+1)=P(-1) : 0.16506       (theory 0.1655)
max |coeff|           : 3             (bound r = 48)
```

The std sits at 0.593, not 0.600, purely because of the discretization
discussed in §6.2 (variance $0.3516$ vs $0.36$). The largest of 2M raw
coefficients is 3 — individual error terms are tiny; it is their
*accumulation* through encryption that must respect $r$.

## Group 3: the correctness budget — and what $\mathrm{erf}$ is doing

This is the heart of the script's analysis, so let us derive everything.

### Where the effective noise comes from

Decrypt an honest ciphertext:

$$
d \;=\; b - a\,s \;=\; (\alpha s + x)y + t + x'' - (\alpha y + x')s
    \;=\; t + \underbrace{x\,y + x'' - x'\,s}_{\text{noise}} .
$$

The message $t$ is what decoding thresholds on; everything else is noise.
Now count the Gaussian contributions per coefficient.

**What we are doing here.** Each of the $n'$ coefficients of the noise
polynomial is a *sum of many independent random terms* — one per way the
individual errors can reach that coefficient through the ring
multiplications. We want the distribution of such a sum, because decoding
correctness depends on whether the sum stays inside $[-r, r]$. The tool
that lets us merge the terms into one number is variance additivity for
*independent* random variables:

$$
\operatorname{Var}[U_1 + \dots + U_m] = \operatorname{Var}[U_1] + \dots + \operatorname{Var}[U_m]
\quad \text{(independent terms)}.
$$

**Term by term:**

* $x\,y$: take one coefficient, say the constant coefficient. It is
  $\sum_i x_i\, y_{(0-i) \bmod n'}$ — the products of *all* error
  coefficients with the ternary entries of $y$ that land on this position
  (with negacyclic signs, which don't matter by symmetry, see below). Of
  the $n'$ terms, exactly $h = 80$ have $y_j = \pm1$ (nonzero); the other
  $n' - h$ terms have $y_j = 0$ and vanish. Each surviving term
  $x_i \cdot (\pm 1)$ is distributed exactly like $x_i$, because
  $D_{\mathbb{Z},\sigma}$ is **symmetric**: multiplying a symmetric random
  variable by $-1$ does not change its distribution. So this coefficient is
  a sum of $h$ independent $\sigma$-Gaussian variables.
* $x'\,s$: identically, a sum of $h$ independent $\sigma$-Gaussian terms
  ($s$ is ternary with the same weight $h$).
* $x''$: one fresh $\sigma$-Gaussian term.

(The three error polynomials $x, x', x''$ and the secrets are drawn
independently — that is what licenses adding variances rather than
worrying about correlations.)

**Adding up:**

$$
\operatorname{Var}[\text{noise}]
= h\,\sigma^2 + h\,\sigma^2 + 1\,\sigma^2
= (2h+1)\,\sigma^2
\quad\Rightarrow\quad
\text{effective std} = \sqrt{2h+1}\,\sigma .
$$

With $h = 80$: $\sqrt{161} \cdot 0.6 = 12.6886 \cdot 0.6 = 7.6131$.

**Why treating the sum as Gaussian is valid.** Strictly, a sum of discrete
Gaussians on $\mathbb{Z}$ is not exactly Gaussian. But the central limit
theorem applies to sums of many independent, mean-zero, light-tailed terms
(with the Berry–Esseen theorem quantifying the gap: the error decays like
$1/\sqrt{m}$ relative to the total spread). Here $m = 2h+1 = 161$ terms of
comparable size — comfortably many — so the accumulated noise is Gaussian
to high accuracy, which is what the tail computation below assumes. (It is
also mildly *conservative* for the far tail: the true discrete tail decays
at least as fast as the Gaussian fit.)

### What $\mathrm{erf}$ is, and what it is for

**The problem it solves.** We will need to answer: "a Gaussian random
variable with std $s$ — what is the probability it lands further than $c$
from its mean?" That is an integral of the Gaussian density,
$\int_c^\infty e^{-z^2/2s^2}dz$, and this integral **cannot be expressed in
elementary functions** (this is a theorem, not a failure of effort — the
antiderivative of $e^{-z^2}$ is provably not expressible by finitely many
basic operations). Mathematics' answer is to *name* the integral and tabulate
it: the **error function**

$$
\mathrm{erf}(u) \;=\; \frac{2}{\sqrt\pi}\int_0^u e^{-t^2}\,dt
\qquad\text{(signed area under } e^{-t^2} \text{ from 0 to } u\text{)},
$$

together with its complement
$\mathrm{erfc}(u) = 1 - \mathrm{erf}(u) = \frac{2}{\sqrt\pi}\int_u^\infty e^{-t^2}\,dt$,
the *upper tail* — exactly the quantity a failure probability needs.

**Purpose in this script.** We hold two numbers: the accumulated noise std
$\sqrt{2h+1}\,\sigma = 7.6131$ and the bound $r = 48$. The correctness
budget is the tail probability $\Pr[|\text{noise}| > r]$, and $\mathrm{erf}$
is the machinery that converts "(std, bound)" into that probability. We are
computing how much probability mass of a bell curve sticks out beyond
$\pm 48$ when the bell has width $7.61$ — about 6.3 standard deviations out.

**Derivation of the formula.** For a normal $G \sim \mathcal{N}(0, s^2)$:

$$
\Pr[|G| \le c]
= \frac{1}{\sqrt{2\pi}\,s}\int_{-c}^{c} e^{-z^2/2s^2}\,dz
\overset{t = z/\sqrt2 s}{=\!=\!=\!=}
\frac{2}{\sqrt\pi}\int_0^{c/\sqrt2 s} e^{-t^2}\,dt
= \mathrm{erf}\!\Big(\frac{c}{\sqrt2\, s}\Big),
$$

where the substitution $t = z/\sqrt2 s$ (so $dz = \sqrt2 s\,dt$) absorbs the
normalization constant and matches the integrand to $e^{-t^2}$, $\mathrm{erf}$'s
definition. Hence

$$
\Pr[|G| > c] = \mathrm{erfc}\!\Big(\frac{c}{\sqrt{2}\,s}\Big)
= 1 - \mathrm{erf}\!\Big(\frac{c}{\sqrt{2}\,s}\Big).
$$

**This is why the argument has exactly the form $r / (\sqrt2 \cdot \sqrt{2h+1} \cdot \sigma)$**:
the random variable being bounded is the accumulated decryption noise with
$s = \sqrt{2h+1}\,\sigma$, the threshold is $c = r$, and the $\sqrt2$ is the
substitution constant that converts a Gaussian CDF into $\mathrm{erf}$.

### The condition

Definition 4.6 sets $r$ so that a decryption fails — noise crosses $r$ — with
probability at most $\epsilon_n/\ell$:

$$
\Pr[\,|\text{noise}| > r\,]
= 1 - \mathrm{erf}\!\Big(\frac{r}{\sqrt2\,\sqrt{2h+1}\,\sigma}\Big)
\;\le\; \frac{\epsilon_n}{\ell} = 2^{-30} = 9.31\times10^{-10}.
$$

(The paper prints the condition as
"$\mathrm{erf}(\ldots) \le \epsilon_n/\ell$", which as written is inverted —
$\mathrm{erf}$ is $\approx 1$ here; the intended statement is the erfc/tail
form above, as the numbers confirm.)

Numerically for Parameter Set 1 ($\sqrt2 \cdot 7.6131 = 10.765$):

| $r$ | $\Pr[\,\lvert\text{noise}\rvert > r\,]$ | vs budget $9.31\times10^{-10}$ |
|---|---|---|
| 44 | $7.49\times10^{-9}$ | over |
| 45 | $3.40\times10^{-9}$ | over |
| 46 | $1.52\times10^{-9}$ | over |
| **47** | $6.68\times10^{-10}$ | **ok (minimum)** |
| 48 | $2.88\times10^{-10}$ | ok ($3.2\times$ margin) |

So $r = 47$ already meets the budget; the table's $r = 48$ buys a further
$\sim\!2.3\times$ tightening ($3.2\times$ inside the budget overall).

Script output:

```
effective noise std   : 7.6131 = sqrt(2h+1)*sigma
P(|noise| > r) approx : 2.884e-10   (need <= 9.313e-10)
```

## Group 4: end-to-end KeyGen demo

```python
s = sample_secret()
alpha = sample_uniform()
pk = alpha * s + sample_error()
noise_coeffs = (pk - alpha * s).lift().coefficients(sparse=False)
```

`pk = alpha*s + x` is computed in the genuine quotient ring — negacyclic
wraparound, mod-$q$ reduction and all — and subtracting $\alpha s$ again
recovers $x$ *exactly*, largest coefficient 2. If Sage's ring arithmetic were
off by even one wrapped sign, this check would explode to $\sim q$-sized
values instead of 2.

# The complete script

```python
#!/usr/bin/env sage
# (see script/sample_errors.sage; abridged listing — the file carries the
#  full commented version with last-run values)

from random import choice
from sage.all import (Zmod, PolynomialRing, RealField, RR, ZZ,
                      erf, ceil, sqrt, cyclotomic_polynomial, Subsets)
from sage.probability.probability_distribution import GeneralDiscreteDistribution

n, nring, m, q = 900, 1024, 2048, 4169729
sigma, h, r, ell, eps_n = 0.6, 80, 48, 1, 2^-30

P = PolynomialRing(Zmod(q), "X"); X = P.gen()
phi = P(cyclotomic_polynomial(m)); R = P.quotient(phi, "x")

PREC = RealField(200)
support = list(range(-ceil(12*sigma), ceil(12*sigma)+1))
weights = [PREC(exp(-x^2/(2*sigma^2))) for x in support]
chi_sigma = GeneralDiscreteDistribution(weights)

def chi_coeff():
    return ZZ(support[chi_sigma.get_random_element()])

def sample_error():
    return R([chi_coeff() for _ in range(nring)])

def sample_secret():
    coeffs = [0]*nring
    for i in Subsets(range(n), h).random_element():
        coeffs[i] = choice([-1, 1])
    return R(coeffs)

def centered(coeffs):
    return [min(ZZ(c), q - ZZ(c)) for c in coeffs]

verify()
```

# Parameter derivations (Table 1, Parameter Set 1)

The paper selects Parameter Set 1 to "minimize computation and communication
while guaranteeing BFV correctness and $\epsilon_n = 2^{-30}$", trading
$\epsilon_p$ up to $2^{-15}$, at 110-bit computational security per the
LWE-estimator. Here is what can be *derived* versus what was *chosen by
search*.

## RLWEenc side (what the script uses)

**$n = 900$, ring dimension $1024$.** The secret-key dimension is a search
parameter balancing security against ciphertext/key size; footnote 10 pads it
to $1024$ for the cyclotomic ring.

**$q = 4169729$ is prime with $q - 1 = 2^{13}\cdot 509$.**
Consequences:

* $\log_2 q \approx 22.0$ — a 22-bit modulus keeps RLWEenc ciphertexts small
  (communication is a first-class cost in OMR).
* $2m = 4096 = 2^{12}$ divides $q-1$, so $\mathbb{Z}_q$ contains $2m$-th
  (in particular $4096$-th) roots of unity: multiplication in
  $R_q = \mathbb{Z}_q[X]/(X^{1024}+1)$ admits a plain NTT implementation
  (negacyclic convolution via a $2m$-point transform). This is the standard
  "NTT-friendly prime" shape.

**$h = 80$, $\sigma = 0.6$.** Sparse ternary secret plus near-ternary noise.
Note $\sigma < 1$ deliberately: it keeps raw errors inside $\{0,\pm1\}$-ish,
which keeps the *accumulated* noise (and hence $r$, and hence the false
negative/positive rates) small. Security with these shapes was certified at
110 bits with the LWE estimator — this is a search result, not derivable in
closed form.

**$r = 48$ — derived.** As computed in §8.3: $r = 47$ is the smallest integer
with $\Pr[|\text{noise}| > r] \le 2^{-30}$; the table's $48$ adds margin
($2.88\times10^{-10}$, about $3.2\times$ inside budget).

**$\ell = 1$.** One payload bit per RLWEenc ciphertext; chosen for the
detection use case (the clue is a flag), and it makes the budget
$\epsilon_n/\ell = \epsilon_n$.

## BFV side (context — not exercised by the sampler)

The detection circuit evaluates a BFV homomorphic decryption, whose noise
budget drives the remaining parameters:

**$D = 2048$ — BFV ring dimension.** The BFV plaintext ring is
$\mathbb{Z}_t[X]/(X^{2048}+1)$, and $D$ is its dimension. Three demands fix
it:

1. **SIMD width.** Because $t \equiv 1 \pmod{2D}$, the ring CRT-decomposes
   into $D$ slots, each carrying one $\mathbb{Z}_t$ value. The detector
   processes $D$ clues *in parallel* (one per slot, §3.3), so $D$ directly
   sets the homomorphic throughput: the amortized per-message cost of the
   partial decryption scales like $n \cdot (\text{BFV ops}) / D$.
2. **The $D \geq n$ constraint** (Algorithm 1 of the paper). The partial
   decryption is the sum $\sum_{i \in [n]} p_i(X) \cdot \mathrm{ct}_{s_i}$
   over all $n = 900$ secret coefficients, with each $p_i(X)$ of degree
   $< D$; the construction needs the ring wide enough to hold these $n$
   coefficient-packed scalar plaintexts, i.e. $D \geq n$. With the
   power-of-two padding $n \to 1024$, this means $D \geq 1024$.
3. **Power-of-two shape.** Same NTT story as §5.1: negacyclic
   multiplication in dimension $D$ runs as an $O(D \log D)$ transform, and
   the slot decomposition itself needs $t \equiv 1 \pmod{2D}$ with $D$ a
   power of two.

Netting out: $D \geq 1024$ from (2), power of two from (3), and doubling
again to 2048 for (1) — twice the SIMD lanes for one extra transform level.
The consistency of the whole design also shows up in the modulus: since
$q = t$ (Algorithm 1) and $q - 1 = 2^{13}\cdot 509$ is divisible by
$2D = 2^{12}$, the *same* prime simultaneously satisfies the RLWEenc NTT
condition ($2m \mid q-1$) and the BFV slot condition ($t \equiv 1
\bmod 2D$) — one modulus, two FFT-friendly roles.

**$Q \approx 2^{60}$** — initial BFV ciphertext modulus, sized for one level
of plaintext multiplication plus $\sim\!900$ additions.

**$h_{\text{BFV}} = 400$, $\sigma_{\text{BFV}} = 3.2$** — BFV secret weight
and error std (3.2 is the customary SEAL-style default; rounded/subgaussian
sampling in practice).

**$Q' = 67104769 = 2^{26} - 2^{12} + 1$ (prime).** The modulus after modulus
switching; $Q' - 1 = 2^{12}\cdot 3\cdot 43\cdot 127$, again NTT-friendly.
$Q'/t \approx 16$ with the plaintext modulus $t$ (see below).

**$r' = 61\cdot16 + 16 = 992$ — derived.** The final error bound after
"merging" via modulus switching. Algorithm 1 of the paper sets

$$
r' \;=\; \frac{Q'}{t}\,r + \mathrm{eModSW} + \lceil Q'/t \rceil ,
$$

Here **eModSW is the error introduced by modulus switching**: before the
final decoding step, the BFV ciphertext is rescaled from modulus $Q$ down
to $Q'$; each coefficient is multiplied by $Q'/Q$ and *rounded* to an
integer, and each rounding is off by at most $0.5$. When the switched
ciphertext is (decryption-)multiplied by the BFV secret, these per-coefficient
rounding errors accumulate over the secret's $h_{\text{BFV}}$ nonzero
positions, giving the bound $\mathrm{eModSW} \le h_{\text{BFV}}/2$ — up to
half a unit of error per nonzero secret entry.

Plugging in Parameter Set 1, with $t = q$ (Algorithm 1 sets the RLWEenc
modulus equal to the BFV plaintext modulus, so $t = 4169729$):

$$
r' \;=\; \frac{67104769}{4169729}\cdot 48 + \frac{400}{2} + \lceil 16.10 \rceil
\;=\; 772.5 + 200 + 17 \;\approx\; 990 \;\approx\; 992 ,
$$

matching the table's $61\cdot16 + 16$ up to the paper's rounding
conventions ($Q'/t \approx 16$ is the "16" in the printed form).

**$\epsilon_p \le 2^{-15}$ — reproduces.** The false-positive rate bound of
Theorem 6.1 is
$\epsilon_p \le \big(2\lfloor t\,r'/Q' \rfloor + 1\big)/t\big)^{\ell}$.
With $r' = 992$: $t\,r'/Q' = 61.65 \Rightarrow \lfloor\cdot\rfloor = 61$,
so

$$
\epsilon_p \;\le\; \frac{123}{4169729} \;=\; 2.950\times10^{-5}
\;=\; 2^{-15.05},
$$

matching the table's $2^{-15}$. (The cruder bound quoted in §5 of the
paper, $((2r+1)/t)^\ell = 97/4169729 = 2^{-15.4}$, is consistent.)

**BFV correctness.** The circuit is one plaintext multiplication
(growth factor $\log(t)\cdot D$) plus $n+1$ additions (factor $\sqrt{n+1}$),
so the final error std is

$$
\log(t)\cdot D \cdot \sqrt{n+1}\cdot\sigma_{\text{BFV}}
\;\approx\; 2^{20}\text{--}2^{22},
$$

which must stay within the post-decoding noise budget
$\lfloor\log Q\rfloor - \lceil\log t\rceil - 1 \approx 60 - 22 - 1 = 37$ bits
(the paper prints the budget as $\lfloor\log Q\rfloor/\lceil\log t\rceil - 1$;
read as a literal ratio it gives a nonsensical 1.7-bit budget, so the
intended reading is almost certainly the subtraction — which the numbers
then satisfy with room to spare).

# Pointers into the papers

* UnifOMR paper: Definition 4.6 (RLWEenc), Section 4.2 (RLWE assumption),
  Theorem 6.1 (error analysis), Section 7.1 + Table 1 (parameters),
  footnote 10 ($n \to$ power of two).
* Lyubashevsky–Peikert–Regev, *On Ideal Lattices and Learning with Errors
  over Rings*: the RLWE foundation — canonical embedding, Gaussian measures
  over $K_{\mathbb{R}}$, and why the error distribution is what it is.
