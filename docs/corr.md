# `stats corr`: design (milestone S2)

**Status: design only, written 2026-10-22 before any code of S2. Nothing in this document has been measured for `corr`; what it reads from the code of S0 and S1 is marked [read], and every bound below is a [choice] fixed before the first run and corrected in place if a run shows it false.** Tags: **[read]** read from a document or source, **[measured]** printed by a script (none yet), **[choice]** a decision the maintainer can overturn, **[open]** not decided. Issue #3.

## 1. What it does

`stats corr FILE [--cols A,B,...] [--method pearson|spearman|kendall] [--missing pairwise|listwise]` reads a CSV with the reader `describe` uses and prints, for every pair of the chosen numeric columns, the coefficient, the number of rows used, the test statistic, the degrees of freedom where there are any, the two-sided p-value, the method of the p-value, and the ties, as one JSON document of the contract (`stats.v1`, `docs/design.md` section 6). The matrix of a table is the list of its pairs `(x, y)` with `x` before `y` in the order of the file: the diagonal is 1 and the lower triangle is the upper one, so neither is printed.

It does not adjust p-values for the number of pairs (`adjustment: "none"` is in the answer; the corrections are `stats test --correct`, S3), and it does not print an interval for the coefficient (Fisher's `z` interval is `stats ci correlation`, S5, where the level and the bootstrap live).

## 2. Definitions **[choice]**

With `n` the rows used for a pair, `Sx = sum x`, `Sy = sum y`, `Sxx = sum x^2`, `Syy = sum y^2`, `Sxy = sum x y`, all exact:

* **Pearson.** `r = N / sqrt(Lx Ly)` with the integers `N = n Sxy - Sx Sy`, `Lx = n Sxx - Sx^2`, `Ly = n Syy - Sy^2` (the `L2` of `describe`, **[read]** `lib/acc/central.cho`). `r` is `sign(N) sqrt(N^2 / (Lx Ly))`, one rounding of an exact rational through the exact integer square root that `describe` already uses for the skewness (`central.sqrt_ratio_bits`, **[read]**). `r` is `null` when `Lx = 0` or `Ly = 0` (a constant column), when `n < 2`; its absolute value is exactly 1 when `N^2 = Lx Ly`, which is decided by the integers, never by the rounded `r`.
* **Spearman.** Pearson's `r` of the ranks, ties given the average rank. A rank doubled is an integer below `2 n`, so the doubled ranks are exact doubles and go through the same exact accumulators (no new arithmetic). The ranks are those of the rows used for the pair.
* **Kendall tau-b.** `tau = S / sqrt((n0 - n1)(n0 - n2))`, `S = C - D` (concordant minus discordant pairs), `n0 = n(n-1)/2`, `n1 = sum t(t-1)/2` over the groups of tied `x`, `n2` likewise for `y`. `S` is an exact integer found in `O(n log n)` by counting inversions with a merge sort after sorting by `x` then `y` (Knight), every count an integer; `tau` is one rounding of `S / sqrt(...)` by the same integer root. A pair tied in both `x` and `y` counts in neither `C` nor `D`.
* **p-values, a [choice] of method, named in the answer (`p_method`).** Pearson and Spearman: `t = r sqrt((n-2)/(1-r^2))` on `n - 2` degrees of freedom, two-sided, `t` found from the exact rational `t^2 = N^2 (n-2) / (Lx Ly - N^2)` rounded once, `p = 2 t_sf(|t|, n-2)` (**[read]** `special.t_sf`). For Pearson this is exact under bivariate normality; for Spearman it is the usual approximation and is named `"t-approximation"`. Kendall: `z = S / sqrt(V)` with Kendall's variance corrected for ties, `V = (v0 - vt - vu)/18 + v1/(2n(n-1)) + v2/(9n(n-1)(n-2))`, where `v0 = n(n-1)(2n+5)`, `vt = sum t(t-1)(2t+5)`, `vu` likewise for `y`, `v1 = (sum t(t-1))(sum u(u-1))`, `v2 = (sum t(t-1)(t-2))(sum u(u-1)(u-2))`; `V` is an exact integer ratio, `z` one rounding, `p = 2 norm_sf(|z|)`, named `"normal-approximation"`. **No exact permutation p-value in S2**: it needs the exact-test machinery of S4 (`stats.exact-too-large`, enumeration limits), and the answer says which approximation it is rather than leaving it to be assumed. **[open]** whether S4 then adds `--exact` here.
* **`null` and why.** An undefined statistic is `null` (as `describe` does) with `reason` in the pair: `"constant-column"`, `"too-few-rows"` (`n < 2` for `r`, `n < 3` for a p-value, `n < 3` for Kendall's `V`), `"all-tied"`.
* **`|r| = 1`.** `t` is infinite; the p-value is `0` and `p_exact_zero: true`, not a number from the tail function. Reached only from the integers.

## 3. Exactness plan

* **One new accumulator, `acc.xsum`**: `sum x y` for two finite doubles. `x y = +-(mx my) 2^(fx + fy - 2148)` with 53-bit mantissas, so the product is a 106-bit integer laid at the bit offset `fx + fy` in units of `2^-2148`: **the same unit as `acc.psum` k = 2**, so `Sxy`, `Sxx` and `Syy` combine in `N` and `Lx Ly` without conversion (**[read]** `lib/acc/psum.cho`: `k = 2` is `m^2` at offset `2 f`, units `2^(-1074 k)`). Positive and negative accumulators, 26-bit limbs, 64 bits of headroom, carry at once: the order of the rows cannot change a bit (the argument of `psum`). This is candidate B of `docs/numerics.md` section 5 (the "cross-product variant"), which that section left for S2.
* **`Sx`, `Sy`** are `psum` k = 1 accumulators of the rows used for the pair; `Sxx`, `Syy` its k = 2. For `n` columns the single sums are per column when listwise and per pair when pairwise (the rows used differ), so pairwise costs `O(n_cols^2 n)` accumulator work; the limit is section 5.
* **Final arithmetic** by `acc.bigint` (**[read]**: multiply, subtract, compare, shift, round with a sticky bit). `Lx Ly` is a product of two numbers of about 4,300 bits.
* **Ranks.** The order is that of `metrics.rank` (the pinned store `describe`'s quantiles use, **[read]**): an order-preserving integer key per double (`-0` equals `0`), an LSD radix sort. The average rank of a group of `t` ties starting at position `p` (1-based) is `p + (t-1)/2`; twice that is an integer.
* **Sizes.** Rows are at most `--max-rows` (ceiling 10,000,000); `n^3` terms of Kendall's `V` reach 1e21 and are carried in `acc.bigint`, not in 64-bit integers.

## 4. Flags, output, refusals **[choice]**

* `--cols A,B,...` (as `describe`): the columns, in the order of the file; default every numeric column. `--method` default `pearson`. `--missing` default `pairwise` (a row counts for a pair if both cells are numbers; `n` is printed per pair); `listwise` drops a row from every pair if any chosen column is missing there. A cell that is not a number makes its column text, as in `describe`.
* Pair: `{"x","y","n","r"` (the key is `rho` for Spearman and `tau` for Kendall`)`,"statistic","df","p_value","p_method","ties":{"x":groups of size above 1,"y":...,"rule":"average"|"tau-b"},"reason"}`. Numbers are decimal text of correctly rounded doubles.
* `--max-pairs` (default 1,000, ceiling 32,640 = 256 columns): more is `limit.too-many-pairs`, checked before any pair is computed. The cell limit of `describe` (`--max-cells`) applies, since ranks and the sort need the cells in memory.
* Refusals with rule tags: `stats.corr-columns` (fewer than two numeric columns chosen; the repair names `--cols`), `stats.non-numeric` (a column named by `--cols` is text or empty; the column is named, no value is printed), `stats.no-such-column` and `stats.ambiguous-column` (as built), the reader's tags, `limit.too-many-pairs`, `limit.too-many-cells`.
* The authority of `stats` does not change: no `conc`, so thread invariance holds by construction and the gate says so (section 5, C-F).

## 5. Gates (written before the code, each able to fail)

| id | gate | fails when |
|---|---|---|
| C-A | **Exact accumulation**: `acc.xsum` against `Fraction` on vectors of the nine families of `scripts/measure_acc.py` (plus constant, subnormal, huge, alternating-sign), pairs of different families, 0 differ after the single rounding and after shuffling | one bit differs |
| C-B | **Oracle bytes**: `r` (Pearson), `rho` of the average ranks, `tau-b` and `S`, the tie counts, `t`/`z` as rounded, each from `Fraction` and the brute-force `O(n^2)` count for Kendall (an independent method, not Knight's), on generated tables of the families of C-A with ties and missing cells, both `--missing` modes | one byte differs; `|r| = 1` decided by the rounded value |
| C-C | **p-values against `mpmath`** at 50 digits for the t and normal tails, from the exact rational `t^2` and the exact `V`. **The bound, fixed now:** relative `1e-12 + (n - 1) 2^-53` for the p-value of Pearson and Spearman (the rounding of `t` is amplified by at most `df + 1` through the tail), `1e-12` plus the rounding of `z` times `z^2` for Kendall, for p-values down to `1e-300` | outside the bound; **a control** (`n` instead of `n - 2` degrees of freedom, a one-sided tail, the untied variance of Kendall on tied data) must fall outside it |
| C-D | **scipy**, reported with the bound fixed here: `pearsonr`, `spearmanr`, `kendalltau(variant="b")` agree on `r` to `1e-12` absolute and on the p-value to `1e-9` relative on well-conditioned families (not the offset-`1e8` family, where scipy's own float sums lose digits: that family is printed as a difference, not gated) | outside the bound on a well-conditioned family; the controls (Welford, a pairwise-deleted rank taken before deletion, average rank replaced by ordinal rank) must fall outside |
| C-E | **Calibration**: under independent normal columns, seeded, 20,000 replicates at `n` = 10, 30, 100, 1,000, the rejection rate at 0.05 inside the binomial bounds and the p-values uniform (KS test) for Pearson at every `n`, for Spearman and Kendall at `n` >= 30; for `n` = 10 the Spearman and Kendall deviations are **printed as measured and documented, not gated** (the approximation is known to be rough there) | a wrong tail or degrees of freedom: the C-C controls must also fail this gate |
| C-F | **Power**: bivariate normal with `rho` = 0.5, `n` = 50, the rate at which Pearson rejects equals the value from the exact density of `r` (Fisher, integrated with `mpmath`) within binomial bounds | outside the bounds |
| C-G | **Invariance**: 50 row permutations give the same bytes; reversing `--cols` gives the mirrored pairs with equal values (`r(x, y) = r(y, x)` byte for byte); no `conc` in the authority (`manifest.py --check`) | one byte differs |
| C-H | **Fuzz**: constant columns, all ties, one row, two rows, all missing, `n = 3`, huge and subnormal values, hostile files, hostile arguments: no trap, no hang, every refusal tagged | a trap or an untagged refusal |
| C-I | **Mutants** for the accumulator, the three coefficients, the tie rules, the p-value branches and every limit: killed or each survivor explained; the real count is printed (`mutants_corr.py` refuses to count kills if the unmutated build fails the gates, as the others do) | an unexplained survivor |
| C-J | **Determinism**: `corr` lines in `docs/xplat_fp.txt`, equal on the runners of the workflow | a line differs |

## 6. Slices (one pull request each)

1. **S2.0, no behaviour change.** `tools/stats/main.cho` is 1,534 lines of the 2,000 allowed (`CLAUDE.md`): the reader loop, `discover` and `read_file` (the header scanner is written twice, which the mutants of 2026-10-22 showed) move to a module shared by `describe` and `corr`. Gate: the existing gates and the fingerprint unchanged.
2. **S2.1** `acc.xsum` and C-A. 3. **S2.2** Pearson with its p-value (C-B, C-C, C-D, C-E, C-F for Pearson). 4. **S2.3** ranks and Spearman. 5. **S2.4** Kendall. 6. **S2.5** mutants, fuzz, fingerprint, the page and the skill text.

## 7. Open points for the maintainer

1. `--missing pairwise` as the default (the usual choice of `pandas.DataFrame.corr` and `cor(use="pairwise")`), against `listwise` (one `n` for the table, ranks computed once).
2. Spearman by the t-approximation and Kendall by the normal approximation in S2, the exact permutation p-values waiting for S4's machinery: confirm, or move S4's enumeration earlier.
3. `--max-pairs` default 1,000 and the pair list instead of a square matrix: a square `matrix` of `r` could be added for `--format text` readers.
4. Whether `acc.xsum`, `acc.psum` and `acc.bigint` become published stores now that a second consumer exists (point 3 of `docs/describe.md` section 6, which had waited for one): proposed yes, after S2.1, so `cancho-ml` can pin them.
5. Whether a Fisher-`z` interval is also printed by `corr`, or only by `stats ci correlation` (proposed: only by `ci`).
