# `stats describe`: design (milestone S1)

**Status: the design was written before any code of S1. Section 7 says what the first slice of the verb does and where it differs from the design; sections 1 to 6 are the design as written, and a paragraph that the code has overtaken says so. What S0 measured is in `docs/numerics.md` and is cited by section.** Tags: **[read]** read from a document or source,
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

## 7. The first slice (built) and what is not yet

`stats describe FILE` exists (`tools/stats/main.cho`, built by `cancho build --bin stats`). It prints, per column: `kind` (numeric, text, empty), `count`, `missing`, and for a numeric column `min`, `max`, `mean`,
`variance`, `sd`, `skewness`, `kurtosis` (excess) and `quantiles`, one JSON document of the contract of cancho-tools' `toolbox` (`stats.v1`, `schemas/stats.v1.json`), with `introspect` and `skill` from the same flag table and the authority report the compiler derives
(`manifests/stats.authority.json`, bounded; `args`, `dir_read`, `err_write`, `file_read`, `fs_read`, `heap`, `io_write`; no `conc`, `clock`, `net` or `ffi`).

**Where the code differs from sections 2 and 3.**
* **`S1` and `S2` both come from `acc.psum`** (k = 1 and 2, own accumulators, gated in `scripts/gate_psum.py`), not from `facc` and `ml.fsq`: the numerator `n S2 - S1^2` needs both sums as limb numbers in
  compatible units (2^-2148), and `ml.fsq` and `facc` publish only a rounded result. `facc` and `ml.fsq` stay pinned for the reader and the measurement of `docs/numerics.md` 8.2.
* **Skewness and excess kurtosis are exact rationals rounded once** (`lib/acc/central.cho`): with `L2 = n C2`, `K3 = n^2 C3`, `K4 = n^3 C4` computed as exact integer combinations of the power sums (acc.psum k = 1 to 4),
  kurtosis is `(K4 - 3 L2^2) / L2^2` (one rounding) and skewness is `sign(K3) sqrt(K3^2 / L2^3)` through an exact integer square root with a sticky bit (one rounding); the moments are the biased ones of section 2. Both are `null` when the variance is zero or `n` is 1.
* **Quantiles are exact** (`lib/acc/quantile.cho`, `--quantiles`, `--quantile-type`): the nine Hyndman-Fan definitions on the sorted numeric cells (an LSD radix sort of order-preserving integer keys, the pinned `metrics.rank`), with the
  probability an exact decimal rational of at most 9 decimals; types 1 to 3 pick an order statistic, types 4 to 9 are `(1 - h) x[j] + h x[j + 1]` evaluated as one exact rational over `24 10^k` and rounded once. The default is type 7
  and the probabilities 0, 0.25, 0.5, 0.75, 1. The cells of the columns are held in memory (`--max-cells`, default 20,000,000, ceiling 100,000,000, one more is `limit.too-many-cells`); `p` is echoed as the text given.
* **`sd` is the square root of the rounded variance**: two roundings, said in the answer (`rounding`), the simpler of the two options of section 2 and the one of open point 2; the tie check of section 2 is not built.
* **`ddof` is 0 or 1**; the variance is `null` when `n - ddof < 1` or when the exact variance is beyond the largest double.
* **A cell that is not a number makes its column `text`**: a quoted cell, an empty-or-not cell the reader of cancho-table (`flt.parse_float`) refuses (NaN, infinity, a bare sign, a space). No cell value is printed.
* **A blank line is skipped, except in a file of one column, where it is one missing cell.** A record that continues on the next line (a quoted newline) is refused (`stats.csv-multiline`), not read.
* Rules of this tool's own: `stats.csv-ragged`, `stats.csv-multiline`, `stats.csv-bad-quote`, `stats.no-header`, `stats.bad-verb` (exit 2), `limit.too-many-rows`, `limit.too-many-columns`, each with exit 8 except the verb.

**Gates of this slice** (all in CI and `scripts/check.sh`): `scripts/gate_describe.py` (generated files of 14 column shapes, every number (including skewness and kurtosis) against `Fraction`, bit for bit, with `--ddof 1` and `0` and a shuffled-row rerun that must print
the same bytes), `scripts/gate_describe_refusals.py` (25 table cases with the exact tag and exit code, and fuzz of hostile files and arguments), `scripts/mutants_stats.py` (26 mutants), `scripts/manifest.py --check`, and the
fingerprint lines `describe` of `docs/xplat_fp.txt`. Measured results are in section 8.

**Not built yet:** the histogram, `distinct` and the most frequent levels of a text
column, `--by`, `--cols`, `--bins`, `--format text`, the tie check of `sd`. The four open points of section 6 keep their proposed defaults and are the maintainer's to change.

## 8. Results of the first slice

*(printed by the scripts named; the machine is the Mac of docs/numerics.md section 8, load as stated there)*

* **Against exact fractions** (`gate_describe.py --files 1000`, Mac): 1,000 files, 54,244 values (count, missing, kind, min, max, mean, variance, sd, skewness, excess kurtosis and the quantiles, default and random type and probabilities, of columns of the shapes uniform, integers, an offset of 1e8 with noise
  below 1, wide exponents (2^-300 to 2^300), tiny (about 2^-1000), huge (2^600 to 2^1000), constant, gappy, empty, text, mixed, quoted, exponent notation), compared bit for bit with `Fraction` rounded once: **0 differ**;
  every file also rerun with its rows shuffled and with `--ddof 0`: the shuffled bytes are equal in every file. The variance of a column of huge values is `null` exactly when the exact variance is beyond the largest double.
  On gram (Linux x86-64, `taskset -c 0-6`): 400 files, 21,540 values, 0 differ. The skewness and kurtosis columns were added after the first measurement of this slice (which had 40,078 values); the numbers above are for the code with them.
* **Refusals and fuzz** (`gate_describe_refusals.py`): 36 table cases, each with the exact rule tag and exit code, and 5,000 fuzz files (Mac; 2,000 on gram): **0 problems**: no trap, no hang, one line of valid JSON, every `error.rule` one
  of the tags of `introspect`, every exit code that of its rule.
* **Mutants** (`mutants_stats.py`): **26 of 26 killed** (the sign and the unit of the mean, the unit of the variance, the divisor, the numerator's subtraction, the limb shift, the comparisons of min and max, the first value, the counts,
  the quoted-cell rule, the square root, the guards, the power of the second accumulator, the blank-line rule, the ragged check, the limit edges, the empty kind, the verb, the delimiter, the `ddof` ceiling, the key of a negative number,
  the multi-line record, the overflow code). `mutants_central.py` (skewness, kurtosis and the exact square root; 29 mutants: the factors and signs of the three integer combinations, the subtraction of `3 L2^2`, the zero-variance guard, the sign and the zero case of the skewness, the number of bits of the square root and of the quotient, the equality case of the root, the sticky bit, the exponent, the shift, the accumulators of the third and fourth powers): **29 of 29 killed**. `mutants_quantile.py` (the nine definitions' constants and branches, the clamps, the exact blend, the sentinel and gather of the cells, the sort range, the probability grammar, the cells limit and the default type; 28 mutants): **28 of 28 killed**. `mutants_psum.py` and `mutants_bigint.py` (S1's accumulators): 19 of 20 and 20 of 21 killed, the survivors explained in the files.
* **The same bytes** (`xplat_fp.py`, 12 `describe` files added to the golden, now 60 lines): equal to the golden, written on the Mac, on gram, and the thread check of the other 48 lines holds on both; the fingerprint does
  not exercise threads for `describe` because the program has no `conc`: its authority has none, and `manifests/stats.authority.json` derived on gram is the one committed from the Mac (`manifest.py --check`).
* **On GitHub** (run 38811179034, `ubuntu-24.04`, `ubuntu-24.04-arm`, `macos-26`): the `describe` gates (200 files against Fractions, 1,000 fuzz files, `manifest.py --check`) and the fingerprint of 60 lines with its thread check
  all succeeded, none skipped.
* **Quantiles against an independent reference** (`gate_quantile.py`, the nine definitions written from the paper's formulas in `Fraction`, the cancho side in the `(a, b)` form so that the two share no algebra): 6,000 cases (Mac and gram) over all nine types, data of the families of `measure_acc.py` with ties and duplicates, 1 to 400 values, `p` from 0 to 1 including the neighbours of the jumps: **0 differ**, bit for bit. Reported only, not gated: 1,303 of 1,500 answers equal `numpy.quantile(method=...)` bit for bit (numpy rounds in floats).
* **Not measured:** speed (no time is asserted); files above 1,000 rows were not generated beyond 400 rows per column in these gates (`--max-rows` defaults to 1,000,000: the memory of the accumulators is 508 integers per column and
  does not grow with rows, but that bound is argued here, not measured at a million rows).
