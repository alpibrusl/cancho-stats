# cancho-stats numerics (milestone S0, issue #1)

**Status: sections 1 to 7 are the design and the gates of S0, written and committed before any code of this milestone. Section 8 holds what the scripts of this repository
printed; a number or a verdict is in it only if a script printed it, with machine, load and the minimum of N.** Tags: **[read]** read from a document or a
source file at a named revision, **[measured]** printed by a script here, **[choice]** a decision the maintainer can overturn, **[open]** not decided.
Written 2026-10-15.

The questions S0 answers, in the order of `docs/design.md` section 5 and issue #1:

1. Which stores does `stats` consume, at which revisions, and how is "pinned, never copied" checked (section 2).
2. How accurate are the special functions the p-values rest on, and what do they do in the far tails (sections 3 and 4).
3. How wide must an accumulator be so that `sum(x^2)` and `sum(x*y)` of arbitrary doubles are exact before the one rounding (section 5).
4. How are the bytes of printed results compared across platforms and thread counts (section 6).

No verb is built here. The code of S0 is a library (`lib/special`) and a driver (`tools/special-check`) that exist to be judged by the gates of section 7.

## 1. Conventions

* **Inputs are doubles, exactly.** A test case is a tuple of binary64 values; the oracle evaluates the function at *those* values in `mpmath` with 60 significant
  digits (`mp.dps = 60`), so the reference is not itself rounded by an earlier step.
* **Error measures.** For a computed `y` and the true value `f` at the same inputs: relative error `|y - f| / |f|`; absolute error `|y - f|`. A
  *bound* is `B(f) = B0 + B1 * |ln f|` for the probability functions (sections 3.2 and 3.3, the second term is the stated far-tail behaviour: the argument of the last
  `exp` carries a rounding error that grows with the size of that argument) and `B0 * max(1, |f|)` for `lgamma`.
* **Out of range is a refusal, never a trap.** Every public function answers `(value, code)`; a `code` other than 0 is a refusal with a rule tag `stats.<name>`
  (section 3.4) and the value is NaN. No input, whatever its bits, reaches a trap, a panic or an unbounded loop.
* **What `exp`, `log`, `log1p`, `expm1` and `pow` promise** **[read]**: `docs/float-math.md` section 9 of the compiler at `a4572ea`: fdlibm's algorithms written in cancho,
  not the platform `libm`; measured against glibc: `exp` 1 ulp, `log` 1 ulp, `expm1` 2 ulp, `log1p` 2 ulp, `pow(x, y)` 2 ulp for `y` of 3.7 or 0.3 and up to 4 ulp for
  `pow(1.7, x)` at `x = +-700`; `sqrt` is a builtin instruction. They are used as they are; their error is part of every bound below and is not hidden.
  There is no fused multiply-add in the language (a product is `a * b`, rounded) **[read]**, so every algorithm below is written for plain IEEE double arithmetic.

## 2. The pins and the pin check

**[choice]** One `cancho.toml` with one `[dependencies.<label>]` block per store, one `rev` per repository, `path = ".cancho-vcs/<module>"`, never a copy.

| repository | stores | rev | state |
|---|---|---|---|
| `alpibrusl/cancho` | the compiler, `cancho = "<rev>"` of `[package]` | `a4572ea3b506e4bcec1c4be43deda1043c11ba9f` | **[read]** the rev `cancho-ml` pins at the start of S0 |
| `alpibrusl/cancho-table` | `dec`, `facc`, `flt`, `pow5`, `reader` | `202b214372a03ec6ed8ca841357f3d2ed4de8aa0` | **[read]** merge commit of cancho-table#44, the rev `cancho-ml` pins (its `docs/reader.md` section 1.3) |
| `alpibrusl/cancho-tools` | the contract package: `cli`, `describe`, `fail`, `limit`, `lines`, `out`, `path`, `place`, `text` (`toolbox.<name>`) | `8dab9176e521660333914dc768a80a29330a2a8e` | **[read]** the rev `cancho-ml` pins; the nine it lists |
| `alpibrusl/cancho-ml` | the 23 numeric modules of cancho-ml#46 (`ml.refuse`, `ml.limits`, `ml.rng`, `ml.exact`, `ml.fbits`, `ml.fsq`, `ml.matrix`, `linalg.{order,tensor,dot,axpy,gemv,syrk,gemm,chol,par}`, `metrics.{rules,ratio,classify,sums,regress,rank,splits}`) | `2de565c99eb5e8f8c28ed984fbcbea3136b5b2f0` | **cancho-ml#46 is not merged and the branch is not pushed** (open at the start of S0, no pull request); the rev is the head of the local branch `v1-stores`, pinned by a `file://` URL of the maintainer's clone. Until it is pushed and merged the pin cannot be fetched from another machine: see 2.3 |

The label of a block is the module name with the dot written as a dash (`ml-rng`, `linalg-chol`), as `cancho-ml`'s `docs/stores.md` section 3 says **[read]**; the toolbox
labels are the module's short name, as in `cancho-ml`'s `cancho.toml`.

### 2.0 The pinned set, as the check reads it

cancho rev: `a4572ea3b506e4bcec1c4be43deda1043c11ba9f`
cancho-table rev: `202b214372a03ec6ed8ca841357f3d2ed4de8aa0`
cancho-tools rev: `8dab9176e521660333914dc768a80a29330a2a8e`
cancho-ml rev: `2de565c99eb5e8f8c28ed984fbcbea3136b5b2f0`

<!-- pins:begin -->
| cancho-tools | cli | toolbox.cli |
| cancho-tools | describe | toolbox.describe |
| cancho-tools | fail | toolbox.fail |
| cancho-tools | limit | toolbox.limit |
| cancho-tools | lines | toolbox.lines |
| cancho-tools | out | toolbox.out |
| cancho-tools | path | toolbox.path |
| cancho-tools | place | toolbox.place |
| cancho-tools | text | toolbox.text |
| cancho-table | dec | dec |
| cancho-table | facc | facc |
| cancho-table | flt | flt |
| cancho-table | pow5 | pow5 |
| cancho-table | reader | reader |
| cancho-ml | ml-refuse | ml.refuse |
| cancho-ml | ml-limits | ml.limits |
| cancho-ml | ml-rng | ml.rng |
| cancho-ml | ml-exact | ml.exact |
| cancho-ml | ml-fbits | ml.fbits |
| cancho-ml | ml-fsq | ml.fsq |
| cancho-ml | ml-matrix | ml.matrix |
| cancho-ml | linalg-order | linalg.order |
| cancho-ml | linalg-tensor | linalg.tensor |
| cancho-ml | linalg-dot | linalg.dot |
| cancho-ml | linalg-axpy | linalg.axpy |
| cancho-ml | linalg-gemv | linalg.gemv |
| cancho-ml | linalg-syrk | linalg.syrk |
| cancho-ml | linalg-gemm | linalg.gemm |
| cancho-ml | linalg-chol | linalg.chol |
| cancho-ml | linalg-par | linalg.par |
| cancho-ml | metrics-rules | metrics.rules |
| cancho-ml | metrics-ratio | metrics.ratio |
| cancho-ml | metrics-classify | metrics.classify |
| cancho-ml | metrics-sums | metrics.sums |
| cancho-ml | metrics-regress | metrics.regress |
| cancho-ml | metrics-rank | metrics.rank |
| cancho-ml | metrics-splits | metrics.splits |
<!-- pins:end -->

(repository, label of the block, module). `scripts/check_pins.sh` holds the same list and PN7 fails when the two differ.

### 2.1 What the check refuses (`scripts/check_pins.sh`)

Each is a failure with its own message, exit 1:

| id | fails when |
|---|---|
| PN1 | a store of the table above has no `[dependencies.<label>]` block, or a block names a store that is not in the table (an unlisted dependency) |
| PN2 | a block is listed twice (the same label, or the same `path` of the same `git`) |
| PN3 | a `rev` is not 40 hex digits, or a repository's blocks do not all carry one `rev`, or the `rev` is not the one of the table in this document (the documented line `<repository> rev: <hash>`) |
| PN4 | a `git` is not the repository's GitHub URL; the one exception is `cancho-ml`, whose `file://` URL is accepted **only while** this document says in section 2.3 that the branch is not pushed |
| PN5 | a `path` is not `.cancho-vcs/<module>` of the block's module |
| PN6 | a `.cho` file in this tree (outside `build/`) declares `module <name>;` for one of the pinned modules (a copy), or a `sources` entry reaches outside the tree |
| PN7 | the documents disagree: the table of this section and the modules the script lists |
| PN8 | `build/deps`, if present, lacks an installed module of the list |
| PN9 | `CHECK_REMOTE=1`: the rev does not have the store on GitHub (skipped, not failed, offline or for the `file://` repository) |

`scripts/mutants_pins.sh` shows each id able to fail: one mutant per row (a block removed, an unlisted block added, a block duplicated, a `rev` of 39 digits, a `rev` that
differs in one block, a `git` pointing elsewhere, a `path` of another module, a module declared in `lib/`, `build/deps` missing a module, the document's `rev` line changed),
each run on a copy, each must exit 1 with the message of its id. The unmutated copy is the control and must exit 0.

### 2.2 Which stores the code of S0 imports

`lib/special` imports `std.math` and its own modules only. The driver `special-check` imports `ml.fbits` (the exact bits-to-double decode of its input), `acc-check` imports `facc`, `flt`
and `ml.fsq` (cancho-table's exact sum and the wide accumulator of cancho-ml; the Dekker product it tests is `lib/special/xprod.cho`, this repository's), and the sweeps are generated in Python
(a SplitMix64 stream with no call to a platform `libm`). The other pinned modules are pinned now because the verbs of S1 to S7 import them and the
issue asks for the pins in S0; a `consume` run (section 7, gate G9) proves that the whole set installs and that a program reaching each family builds.

### 2.3 The `file://` pin of the numeric stores **[open until cancho-ml#46 is pushed and merged]**

The `cancho-ml` URL in `cancho.toml` is `file:///Users/alfonso/Workspace/alpibrusl/lexsys/cancho-ml-stores`. It installs on the machine that has that clone (**[measured]**
2026-10-15 on the Mac: `cancho install` fetched 39 files for all the blocks). On any other machine `cancho install` fails for these 23 blocks. Consequences, stated so that nothing
is hidden: gram builds with a `git config url.<clone>.insteadOf <the URL>` pointing at a clone of the same commit (section 6); the CI workflow runs the pin check and every gate that
needs no store and probes the pin (`git fetch` of the rev) before the gates that do, and reports "skipped: the rev of cancho-ml is not fetchable" instead of passing silently. When
cancho-ml#46 is merged the URL becomes `https://github.com/alpibrusl/cancho-ml`, the rev the merge commit, PN4's exception is removed from the script and this subsection is deleted
(in one commit).

## 3. Special functions: design and accuracy target

**This section is written before the first line of `lib/special`.** Algorithms are a **[choice]**; the bounds are the contract the gates enforce.

### 3.1 What is implemented, and from what

| function | public name | method |
|---|---|---|
| `ln Gamma(x)`, `x > 0` | `special.lgamma` | Stirling series for `x >= 10` (coefficients `B2k / (2k (2k-1))`), upward recurrence `lgamma(x) = lgamma(x + n) - ln(x (x+1) ... (x+n-1))` below, reflection never (the domain is positive) |
| `erf`, `erfc` | `special.erf`, `special.erfc` | `erf(x) = sign(x) P(1/2, x^2)`, `erfc(x) = Q(1/2, x^2)` for `x >= 0`, `erfc(-x) = 2 - erfc(x)` |
| `P(a, x)`, `Q(a, x)` (regularised incomplete gamma) | `special.gamma_pq` | power series for `x < a + 1`, Lentz's continued fraction otherwise; the prefactor `x^a e^-x / Gamma(a+1)` by Loader's saddle-point form (`stirlerr`, `bd0`) for `a >= 15`, by `exp(a ln x - x - lgamma(a+1))` below |
| `I_x(a, b)` (regularised incomplete beta), both tails | `special.beta_ix` | Lentz's continued fraction on the side where it converges (`x < (a+1)/(a+b+2)`, else the swapped arguments); input is the pair `(x, 1-x)` so neither is rounded; the prefactor `x^a (1-x)^b / B(a, b)` by Loader's binomial form for `a, b >= 10`, by `exp(a ln x + b ln(1-x) - betaln(a, b))` below, `betaln` by the Stirling-difference form when one argument is large |
| normal quantile | `special.qnorm` | Acklam's rational approximation as a start (**[choice]**, its published relative error is about `1.15e-9`), then three Newton steps on `ln Phi(z) - ln p` (the log domain, so that a subnormal `p` does not lose its bits) |
| normal, Student t, chi-square, F, binomial (cdf and survival) | `special.dist` | the functions above: `Phi(z) = erfc(-z/sqrt2)/2`, `t` and `F` through `beta_ix`, chi-square through `gamma_pq`, binomial `P(X <= k) = I_{1-p}(n-k, k+1)` |

Both tails of every distribution are computed on their own: a survival function is never `1 - cdf`, so it does not lose a tail. `beta_ix` and `gamma_pq` answer
the pair `(lower, upper)`; the one that comes from the continued fraction is computed directly, the other is `1 - it` and is the larger of the two (at least 1/2 up to
the switch point, so it carries no cancellation).

### 3.2 Accuracy target (the bounds the gates enforce)

`B(f) = B0 + B1 |ln f|` relative to `f` unless the metric says absolute. `B1 = 4e-16` everywhere it appears (about two ulps of `|ln f|`: the rounding of the
argument of the last `exp`). Ranges are of the *inputs*, sampled as in 7.1.

| function | input range | metric | `B0` | notes |
|---|---|---|---|---|
| `lgamma(x)` | `x` in `[1e-300, 1e15]` | `abs`: `B0 * max(1, abs f)` | `1e-14` | near `x = 1` and `x = 2`, where `f = 0`, the bound is the absolute `1e-14` |
| `erf(x)` | `abs x` in `[1e-300, 6]` | relative | `2e-13` | for `abs x` below `1e-8`, `erf(x)` is `2x/sqrt(pi)` to the last bit of the product |
| `erfc(x)` | `x` in `[-6, 26.5]` | relative | `2e-13` | `B1` term applies; from `x = 26.5` the result is below `1e-306` and the bound is waived, see 3.3 |
| `P(a, x)`, `Q(a, x)` | `a` in `[1e-2, 1e4]`, `x` in `[1e-300, 1e5]` | relative | `2e-13` | the smaller tail included (any `P` or `Q` above `1e-300`) |
| same | `a` in `(1e4, 1e7]`, `x` within `a +- 12 sqrt(a)` | relative | `5e-12` | beyond 12 standard deviations the relative sensitivity to the rounding of `n * p`-type products (about `a eps` times the deviation) is of the size of the bound; the sweep stops there and the document says so |
| same | `a` in `[1e-300, 1e-2)` | relative | `2e-12` | |
| `I_x(a, b)` both tails | `a`, `b` in `[1e-2, 1e3]` | relative | `5e-13` | |
| same | `a` or `b` in `(1e3, 1e6]`, the other in `[1e-2, 1e6]` | relative | `1e-11` | `x` within `a/(a+b) +- 12 sd` of the mean |
| normal `cdf`, `sf` | `z` in `[-38, 38]` | relative | `2e-13` | |
| Student `t` cdf, sf | `nu` in `[0.1, 1e6]` and `nu = 1e8`, `t` such that the tail is at least `1e-300` | relative | `1e-12` | |
| chi-square cdf, sf | `nu` in `[1e-2, 1e6]` | relative | `2e-13` for `nu <= 1e4`, `5e-12` above (`x` within `nu +- 12 sqrt(2 nu)`) | |
| `F` cdf, sf | `d1`, `d2` in `[0.1, 1e5]` | relative | `1e-12` | |
| binomial cdf, sf | `n` in `[1, 1e6]`, `p` in `[1e-6, 1 - 1e-6]`, integer `k` | relative | `5e-12` | |
| `qnorm(p)` | `p` in `[1e-300, 1)`, and subnormal `p` down to `5e-324` | absolute plus relative | `1e-15 + 1e-13 * abs x` | `p > 1/2` uses `-qnorm(1-p)`, exact subtraction |

### 3.3 Far tails and the edges, stated before measuring

* **The error grows like `ln f`, not like `f`.** The `B1` term is the rule: at `f = 1e-300` (`|ln f| = 690`) the bound is `B0 + 2.8e-13`. This is the honest cost of
  one `exp` of a rounded argument; it is **[open]** whether a split `exp(-x^2)` (a high part with an exact square) would remove it for `erfc`; not done in S0 (a choice: the bound above is
  met without it or the bound is corrected in place).
* **Underflow.** Every bound carries an absolute allowance of `2e-323` (four units of the smallest subnormal): `|y - f| <= B(f) |f| + 2e-323`. It matters only below the smallest normal double
  (`2.2250738585072014e-308`), where a result has fewer bits; above it the allowance is below the bound by a factor of at least `1e10`. A true value that is zero answers zero exactly, and a true
  value below `2.5e-324` answers 0 (not NaN, not a refusal: the correct answer rounded to the representable). *Corrected in place:* the first version of this bullet waived the relative bound below the
  smallest normal and allowed only four units; the first sweep that reached `1e-309` showed a relative error of `1.8e-13` there (36 units), inside `B(f)`, so the rule is the relative bound plus the allowance.
* **Edges.** `P(a, 0) = 0`, `Q(a, 0) = 1`; `x = +inf` gives `P = 1`, `Q = 0` (a finite `a`); `I_0 = 0`, `I_1 = 1`; `qnorm(0.5) = 0` exactly and `qnorm(p) = -qnorm(1-p)` bit for bit for `p` in `[1/2, 1)`;
  `erf(0) = 0`, `erf(+-inf) = +-1`; NaN in any argument is a refusal `stats.domain-nan`; a parameter that must be positive and is not is `stats.domain-param`; a probability outside `(0, 1)` is `stats.domain-prob`
  (`qnorm(0)` and `qnorm(1)` are refusals, not infinities).
* **Non-convergence.** The continued fractions and the series have an iteration cap derived from the arguments (`2000 + 100 * sqrt(a)` terms for the incomplete gamma, `2000 + 100 * sqrt(max(a, b))` for the incomplete beta; each is a named function `cap`); reaching it is the
  refusal `stats.noconverge` with value NaN, never a loop that does not end. The edge gate reached it once (`gamma_p` with a tiny `a` and `x` near the largest double, where a subnormal first
  term kept the fraction from settling); the prefactor underflow test removes that input, and no input of the gates reaches the refusal now (section 8.3). The mutant that sets the cap to 1 shows
  that the path is wired.
* **Ties between branches.** The series/continued-fraction switch at `x = a + 1`, the swap of `beta_ix` at `x = (a+1)/(a+b+2)`, Loader against direct prefactor at `a = 15` or `a, b = 10`, the Stirling switch at `x = 10`,
  Acklam's tail switch at `p = 0.02425`: the sweeps contain the boundary value, its predecessor and its successor (`nextafter` by integer arithmetic on the bits) and the bound must hold on all three; this is where a branch
  that is wrong by a constant shows.

### 3.4 Rule tags **[choice]**

`stats.domain-nan`, `stats.domain-param`, `stats.domain-prob`, `stats.domain-range` (an argument beyond the range of 3.2 for which the functions do not promise the bound: a refusal, e.g. `a > 1e7`),
`stats.noconverge`. A refusal never carries input data in its tag. Every refusal the driver prints is `refused <tag>`; the gate of 7.2 checks that each tag in `special.rules` other than `stats.noconverge` is reached by at least one case and
that no other tag is printed; `noconverge` is the one tag no input is known to reach (see 3.3). The refusal limits (`domain-range`): `a` above `1e8` for the incomplete gamma, `a` or `b` above `1e8` for the
incomplete beta, `nu` or a degree of freedom above `2e8` for t, chi-square and F, `n` above `1e8` for the binomial, `x` above `1e300` for `lgamma`. Between the limits of 3.2 (where the bounds are measured)
and these, the functions answer and the gates do not bound the answer.

## 4. How the special functions are judged

The reference is `scripts/special_ref.py` and the gate `scripts/gate_special.py` (a venv inside the clone: `mpmath`, `numpy`, `scipy`; the venv is not committed, `ci/requirements.txt` lists the versions). `mpmath`'s own `gammainc` and `betainc` fail to converge or take minutes for parameters of 1e4 and above **[measured, 2026-10-15, mpmath 1.4.1]**, so for those two the reference evaluates the defining series and continued fraction in `mpmath` arithmetic at 100 digits (more when the exponent is large); that reference is itself gated by `special_ref.py --selfcheck`: it must agree with `mpmath`'s native functions to 1e-30 relative wherever the native value is at least 1e-20 (below that `mpmath` itself returns noise, even a negative number, **[measured]**), and with `scipy` in G3. `numpy` and `scipy` appear in
**one** place: the *independence check* of 7.1 item 5, a second reference to catch an error shared by `mpmath`'s own evaluation at a few thousand points; **no sentence of this repository says a result "matches" scipy beyond the bound
fixed here and measured.**

The sweeps (7.1) are generated by the oracle from a seed (SplitMix64 with the derivation of `docs/rng.md` of cancho-ml, rewritten in Python, **[read]**), written as the bit patterns of the inputs, evaluated by the driver
`tools/special-check`, and the driver's output (bit patterns of the results and the codes) is compared by the oracle. The comparison is by `mpmath` at 60 digits, so the measured maximum is a true error and not a difference of two libraries.

## 5. Exact products for variance and covariance

**The question (open in `docs/design.md` section 5):** `facc` (cancho-table, **[read]**: 72 limbs of 32 bits, an integer multiple of `2^-1074`, 73 integers per accumulator, 2^27 additions between carries, nothing rounds
until `finalize_key`) holds any sum of finite doubles exactly. The variance needs `sum(x^2)` and the covariance `sum(x*y)`; a product of two doubles is an integer of up to 106 bits times `2^(ex + ey)`, and
`facc` has no room for `2^-2148`. What must be added, and when is it needed?

**The candidates:**

* **A. Two-product into `facc`.** `x * y = p + e` exactly (Dekker, `std.math`'s `two_product`, no fma), then `facc.add_bits(p)` and `facc.add_bits(e)`. Exact **iff** the pair is exact: no overflow in the split
  (`abs x < 2^996`), `p` finite, and `e` not below the unit of `facc` (`e` is a multiple of `2^(ex + ey - 104)`, which must be at least `2^-1074`). This needs a **range rule**, and a refusal for what is outside it.
* **B. A wide limb accumulator.** The exact product laid into limbs of an integer accumulator of `2 * 2098 + 106` bits in units of `2^-2148` (the width `ml.fsq` already uses for `(y - p)^2`: 528 integers, **[read]**) with a
  cross-product variant. Exact for every pair of finite doubles whose product is finite, with no range rule, and costs more per element.
* **C. A running mean (Welford) or two passes** is excluded by `docs/design.md` section 3 (the last bits would depend on the row order); recorded here only so the measurement shows how far it is from exact.

**The measurement (written first):** `scripts/measure_acc.py` generates, from a seed, vectors of doubles in nine families (uniform in `[0,1)`, wide exponent (`2^e`, `e` uniform in `[-1074, 1023]`), large (`2^e`, `e` in `[500, 1000]`),
small (`e` in `[-1074, -480]`), mixed signs with cancellation (`x`, `-x`, `1e-30`), a constant column, subnormals only, near-overflow products, and the cases of `ml.fsq`'s own tests), computes
`sum(x^2)` and `sum(x*y)` **exactly** with Python's `fractions.Fraction` and rounds once, and compares against the cancho results of candidate A, candidate B (through `ml.fsq` for squares) and
candidate C. The tables report, per family and candidate, the number of vectors exact to the bit, the number refused, and the largest relative error of the rest. **The rule that follows is written in section 8.2
after the run**, as one of: "A with the range rule R" (R stated as a closed interval of exponents, with the test that it is exactly the exactness boundary: a vector one step outside must be shown inexact), "B", or
"both: A inside R, B outside". The one-rounding statement is fixed now: **the result is `round_to_nearest_even(exact sum)`, one rounding, at the end, and no intermediate rounding of any term or partial sum.**
The gate G4 (7.1) fails if any vector inside the rule's range is not the exactly rounded value.

## 6. Determinism

**What is compared:** the *bytes printed* by the driver: for each function and each sweep, one line per case of the form `<case id> <hex bits of the result> <code>`, hashed. `scripts/xplat_fp.py` prints the
**fingerprint**: one tab-separated line per (function, sweep) `name  count  sha256-of-the-lines`, plus the accumulator experiments' results as exact hex bits, in the style of cancho-ml's `scripts/xplat_check.py` **[read]**.
The **golden** is `docs/xplat_fp.txt`, written on the maintainer's Mac (arm64) and read back on every other platform; the check is *equal to the golden*, so equal to each other.

**Platforms:** macOS arm64 (this Mac); gram, Linux x86-64 (`ssh gram`, a clean clone under `~/Workspace/xplat-ml`, **`taskset -c 0-6` only**: one core there is busy with someone else's work); the
Raspberry Pi (`ssh pi`) **only if the maintainer says so**. The CI workflow (section 7, G7) runs the same check on `ubuntu-24.04` and `macos-15`-class named images (the exact names are in the workflow, `-latest` never).
What the golden shows is stated with its conditions in 8.4: *equal bytes on these platforms for these sweeps*, nothing about others.

**Threads:** the functions are pure and read no shared state, so the thread count cannot change a bit by construction; this is **checked, not assumed**: the driver's `--threads T` mode cuts the case list into `T`
contiguous blocks, evaluates each block on its own thread (the `spawn`/`join` pattern of `linalg.par` **[read]**), and hashes the results in case order. `T` in `1..16`: the fingerprint must equal the `T = 1` line for
every `T`. A mutant (a result slot written in completion order) shows the check able to fail. On gram the threads run on cores 0 to 6 (7 cores, so `T = 8..16` oversubscribe: the claim is about bytes, not speed).
The accumulators of section 5 are sequential in S0; their row-order invariance is tested by shuffling the rows (exact sums: every permutation gives the same bytes).

## 7. The gates (written before the code; each can fail)

### 7.1 Accuracy gates

| id | gate | fails when |
|---|---|---|
| G1 | **Sweeps against `mpmath`**: per function and range of 3.2, at least 20,000 cases (2,000 for the functions whose `mpmath` evaluation is slow: `betainc` with a parameter above `1e3`, the large-`a` gamma), half log-uniform in the parameters and half concentrated on the transition region (`x` near `a`, near the mean of the beta, near the Stirling and Acklam switches, and the `+-1 ulp` neighbours of every branch tie of 3.3) | any case exceeds `B(f)`; the report prints the maximum error per range and the case that attains it |
| G2 | **Controls**: the sweep is run against a driver built with a deliberately wrong algorithm: a Stirling coefficient `1/12` changed to `1/12.001`, the series cut at 12 terms, the Lentz stopping tolerance loosened from `2.5e-16` to `2.5e-6`, `Q = 1 - P` in the survival function, the beta swap threshold moved by 5%, Acklam's constant changed in the 4th digit **and the Newton steps removed**. For each control the gate G1 must **fail**, in the range the control is for | a control stays inside the bound: the bound is too loose or the sweep does not reach the code (the report says which) |
| G3 | **Independence**: 4,000 points per function against `scipy.special` (`gammaln`, `erf`, `erfc`, `gammainc`, `gammaincc`, `betainc`, `ndtr`, `ndtri`, `stdtr`, `chdtr`, `fdtr`, `bdtr`) and `mpmath`, reported separately; the sweep reports the disagreement of the two *references*, so an `mpmath` evaluation error cannot be hidden by the driver agreeing with it | the references disagree beyond `1e-10` relative on a point (the point is reported and removed from G1 only with a written reason) |
| G4 | **Exact accumulation**: section 5, per family, per rule | an exactly rounded value is not reproduced bit for bit inside the rule's range |
| G5 | **Mutants** for every function (7.3) | a survivor without an explanation |

### 7.2 Domain, edges and fuzz

| id | gate | fails when |
|---|---|---|
| D1 | the edges of 3.3 for every function: `0`, `1`, `2^-1074`, `2^-1022`, the largest double, `+-inf`, a NaN bit pattern (every NaN payload: the driver decodes `bits` itself), `-0.0`, the neighbours of each branch tie | the process traps, hangs (a 10 s timeout per case batch, the driver is killed) or answers a value where a refusal is specified (or the reverse) |
| D2 | fuzz: 200,000 random 64-bit patterns per argument position, each function, uniformly over the *bit patterns* (so exponents are uniform: huge, tiny, subnormal, inf, NaN all appear) | any trap or timeout; any value outside the contract of the function (cdf outside `[0, 1]`, `cdf + sf` not within `1e-12` of 1 where neither is below `1e-10`, the quantile not monotone on a random sorted set) |
| D3 | every tag of `special.rules` is reached; no other tag is printed; no refusal prints input data | a tag unreached or unknown |

### 7.3 Mutants

`scripts/mutants_special.py`, in the style of `mutlib.py` of cancho-ml: a mutant is `(name, file, old, new)`, `old` occurs exactly once, applied to a copy, the driver rebuilt, the **fast** sweeps (G1 at 600 cases per range,
deterministic) rerun; killed = the gate fails, or the driver traps. At least five mutants per public function (an operator or a constant, a branch condition off by one on a tie, an iteration cap, a dropped term, a sign) and
the survivors **each explained in the file** (an equivalent mutant: e.g. a change in a cap that no sweep reaches is not equivalent, it is a hole: the sweep is extended, not the survivor explained).

### 7.4 Repository gates

| id | gate | fails when |
|---|---|---|
| G6 | `scripts/check_pins.sh` and its mutants (section 2) | see 2.1 |
| G7 | `.github/workflows/ci.yml`, modeled on `alpibrusl/cancho-ml`'s: named runner images (never `-latest`), actions pinned by full commit SHA, `permissions: contents: read`; `scripts/check_workflow.py` and its mutants | any `uses:` without a 40-hex SHA, any `-latest`, any permission other than read |
| G8 | the fingerprint of section 6: equal to the golden on this Mac and on gram; `T = 1..16` equal | a differing line |
| G9 | `scripts/consume.py`: a scratch package (outside the repository) with every block of `cancho.toml`, `cancho install`, and a program that imports and calls one function of each family (the numeric stores, `facc`, `flt`, `reader`, the contract package) and prints known values | an install, build or value failure; mutants: a block removed, a store byte changed |
| G10 | no file over 2,000 lines; docs have no changelog narration (`scripts/check_docs.sh`: the sentences of the document are in the present tense, no "was added", no dates in prose other than the status line) | a file at 2,001 lines |

## 8. Results

*(filled in only from what the scripts printed; empty until they run)*
