# `stats describe`: design (milestone S1)

**Status: design, written before any code of S1. Nothing here is measured for `describe`; what S0 measured is in `docs/numerics.md` and is cited by section.** Tags: **[read]** read from a document or source,
**[measured]** printed by a script (S0's), **[choice]** a decision the maintainer can overturn, **[open]** not decided. Written 2026-10-15.

## 1. What it does

`stats describe FILE [--by COL] [--cols A,B] [--quantile-type 1..9] [--ddof 0|1] [--bins N]` reads a CSV with cancho-table's reader (the pinned `reader`, `flt`, `facc` stores, **[read]** `docs/numerics.md` 2) and prints one JSON
document (`schema stats.v1`, `docs/design.md` section 6): per numeric column `count`, `missing`, `distinct`, `min`, `max`, `mean`, `variance`, `sd`, `skewness`, `kurtosis`, the quantiles (default 0, 0.25, 0.5, 0.75, 1) and an
exact histogram; per non-numeric column `count`, `missing`, `distinct` and the most frequent levels. With `--by`, the same per group. Every number is the decimal text of a correctly rounded double (the printer of cancho-table, **[read]**).

## 2. Definitions **[choice]**

* **mean** = `S1 / n`, one rounding. **variance** = `(n S2 - S1^2) / (n (n - ddof))`, exact numerator, one rounding; default `ddof = 1` (the sample variance, as R and pandas; numpy's default is 0), refused with
  `stats.sample-too-small` when `n - ddof < 1`. **sd** = the correctly rounded square root of the *rounded* variance is a second rounding; instead sd is `sqrt` of the exact variance rounded once to 106 bits
  (double-double) and then rounded: **[open]** whether the single-rounding claim covers sd (the square root of a rational is not a rational; the answer is correctly rounded only if the 106-bit value is not within 2^-106 of a tie, which is checked and, if it happens, refused as `stats.sd-tie`).
* **skewness** = `m3 / m2^(3/2)` and **excess kurtosis** = `m4 / m2^2 - 3`, with `m_k = (1/n) sum (x - mean)^k` (the biased moments; the adjusted G1 and G2 are flags `--skew-kind`): the central moments are exact rationals
  from the power sums `S1..S4` (`n^k m_k` is an integer combination of `S_j S1^(k-j)`), so the rounding happens in the root and the quotient only, and the answer says which (`rounding: "exact-then-2"`).
* **quantiles**: the nine Hyndman-Fan types. Types 1 to 3 select an order statistic exactly; types 4 to 9 interpolate `x[j] + g (x[j+1] - x[j])` with `j` and the rational `g` from integers (`n`, `p` given as a decimal
  rational): evaluated exactly, rounded once. Default type 7 (R, numpy) **[choice]**. Order statistics by an exact selection (no float comparison issue: a total order on the doubles, `-0` before `0` is not distinguished).
* **histogram**: `--bins N` equal-width bins over `[min, max]` with edges that are the rounded rationals `min + i (max - min) / N` and counts by exact comparison against the *rational* edges (not the rounded ones), the last bin closed.

## 3. Exactness plan

* **Sums.** `S1 = sum x` by `facc` (2304 bits, **[read]**). `n` by a count, missing cells excluded and counted.
* **Power sums.** `S2 = sum x^2`: `ml.fsq` via `add_square(acc, bits, 0)`, exact for every finite double (**[measured]** `docs/numerics.md` 8.2). `S3` and `S4` are not provided by a pinned store: a module `acc.psum` in this
  repository lays the integer mantissa `m` (53 bits) and its power `m^k` (at most 212 bits for k = 4) into limbs at the bit offset `k e` (the same scheme as `ml.fsq`'s, which is the reference for its limb size and carry rule),
  width `k * 2098 + 53 k + 64` bits per accumulator, no carry until the final pass. **[open]** 26-bit limbs as `fsq` or 32-bit as `facc`; decided by the measurement of section 5, gate S1-A.
* **Final arithmetic.** A small unsigned big-integer module `acc.bigint` (multiply of two limb numbers, subtraction, comparison, shift, division by a number of up to 128 bits, rounding to 53 bits with a sticky bit):
  `std.bignum` of the pinned compiler has none of these (**[read]**, `docs/numerics.md` 8.2). Everything here is integers and shifts: no float operation, so the bytes are the same on every platform and thread count.
* **Order.** All accumulators are sums of integers: row order and thread partition do not change a bit; per-thread accumulators are merged limb by limb (exact, `facc`'s own argument, **[read]**).
* **Ranges.** Exact for every finite double; a sum beyond the largest double is `stats.sum-overflow`. No rule R is needed for `describe` (that is for cross products, S2).

## 4. Refusals (every one has a `stats.` tag, an exit code, a repair where safe)

`stats.sample-too-small` (n - ddof < 1), `stats.variance-zero` (skewness, kurtosis of a constant column: `m2 = 0`), `stats.sum-overflow`, `stats.non-numeric` (a column asked for as numeric has a cell that is not a number: the row and
column are named, no cell value is printed), `stats.quantile-type`, `stats.bins-limit`, `stats.rows-limit`, `stats.cols-limit`, plus the reader's, the limits' and the arguments' tags of the contract package. Bounded resources: rows, columns,
distinct levels and bins each have a flag, a default and a ceiling, checked before anything sized by the input.

## 5. Gates (written before the code, each able to fail)

| id | gate | fails when |
|---|---|---|
| S1-A | the accumulators: `acc.psum` and `acc.bigint` against `fractions.Fraction` on vectors of nine families (the families of `scripts/measure_acc.py`, plus constant, subnormal, huge, alternating-sign and 1e6 rows), the width decided by this measurement | one bit differs after the single rounding |
| S1-B | oracle bytes: count, min, max, mean, variance, skewness numerator, the quantiles of types 1 to 9, the histogram counts, distinct, missing: the printed decimal text equals `Fraction`-computed text for 20,000 generated columns | one byte differs |
| S1-C | quantile types against an independent Python implementation of the Hyndman-Fan definitions (written from the paper's formulas, not numpy's code) and, reported separately, against `numpy.quantile(method=...)`; the bound is exact equality of the rounded rational | a type differs from the definitions |
| S1-D | row-order and thread invariance: 50 random permutations and 1 to 16 threads give the same bytes | one byte differs |
| S1-E | determinism: a fingerprint of 40 columns extends `docs/xplat_fp.txt`; equal on the platforms of the workflow | a line differs |
| S1-F | fuzz: hostile CSVs (empty, one row, all missing, constant, huge, subnormal, NaN-like and inf-like cells, ragged rows, 1e6 rows): no trap, no hang, every refusal tagged | a trap or an untagged refusal |
| S1-G | mutants for every statistic and every quantile type, killed or each survivor explained | an unexplained survivor |

## 6. Open points for the maintainer

1. `ddof` default 1 (sample variance) and biased moments for skewness and kurtosis: confirm or change (section 2).
2. Whether `sd` may be a second rounding (rounding the exact variance first) instead of the tie check of section 2: simpler, but then "one rounding" is not true for `sd`, and the answer would say `rounding: "exact-then-2"`.
3. Whether `acc.psum` and `acc.bigint` should be published as stores for S2 to pin, or stay in this repository until a second consumer exists (the rule of `cancho-ml/docs/stores.md`: a module stays out until a consumer needs it).
4. The histogram's edge convention (closed last bin, exact comparison against rational edges).
