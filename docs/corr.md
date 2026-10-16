# `stats corr`: design (milestone S2)

**Status: written before any code of S2. Section 7 says which slice is built and where it differs.** Tags: **[read]** read from a document or source, **[measured]** printed by a script, **[choice]** a decision the
maintainer can overturn, **[open]** not decided. Written 2026-10-16. What S0 and S1 measured is in `docs/numerics.md` and `docs/describe.md`.

## 1. What it does

`stats corr FILE [--cols A,B,...] [--method pearson,spearman,kendall] [--missing pairwise|listwise] [--max-vars N] [--format json|text]` reads a CSV and prints, for each method asked, the matrix of correlations of the
selected numeric columns, the number of rows each pair used, and the two-sided p-value of the null of no association, in one JSON document (`stats.v1`, the contract of `describe`).

## 2. Definitions **[choice]**

* **Pairs.** `--missing pairwise` (default): a row counts for the pair `(a, b)` if both cells are numbers; `n_ab` is the count and is printed. `listwise`: only rows where every selected column is a number. A cell that is not a number in a
  SELECTED column makes the column non-numeric: `stats.non-numeric` (the column is named; no cell value is printed). `--cols` omitted: every column that is numeric in the file (at least 2, else `stats.corr-too-few-columns`).
* **Pearson** `r = (n Sxy - Sx Sy) / sqrt((n Sxx - Sx^2)(n Syy - Sy^2))` from the EXACT sums of the doubles (`acc.psum`, with a cross-product accumulator added in S2a): the numerator and both factors of the denominator are exact
  integers, so `r` is `sign sqrt(num^2 / (Dx Dy))`, ONE rounding (the exact integer square root with a sticky bit of S1). No rule R: exact for every finite pair of doubles whose sums are finite. `r` is `null` when `Dx` or `Dy` is 0
  (a constant column: the correlation is undefined, and the notes say which pair), or when `n_ab < 2`.
* **p-value, Pearson and Spearman** (the t approximation, said in the answer as `p_method: "t"`): `t^2 = r^2 (n - 2) / (1 - r^2)`, so with `nu = n - 2` the tail is the beta tail `I_x(nu/2, 1/2)` with `x = nu / (nu + t^2) = 1 - r^2` and
  `y = r^2`. Both `x` and `y` are EXACT rationals here (`(Dx Dy - num^2) / (Dx Dy)` and `num^2 / (Dx Dy)`), each rounded once, and enter `special.beta_xy` as the pair it was built for, so the p-value carries the forward bound of S0
  (5e-13 relative) and nothing from `1 - r^2` cancelling when `|r|` is near 1. `p` is `null` when `n_ab < 3`; it is 0 when `|r| = 1` and `n_ab > 2`.
* **Spearman** is Pearson of the average ranks within the rows of the pair. The rank of a value is the mean of the ranks of its ties; `2 rank` is an integer below `2 n`, a double, so the same exact machinery applies to the ranks.
  The tie rule is `average` and is printed. Ranks are computed per pair over its own rows (for `pairwise`), so a pair's coefficient is a function of that pair's rows only.
* **Kendall tau-b** `tau_b = (C - D) / sqrt((n0 - n1)(n0 - n2))` with `C`, `D` the concordant and discordant pairs, `n0 = n(n-1)/2`, `n1` and `n2` the numbers of tied pairs in each column, by Knight's `O(n log n)` count over
  integers (exact; `n <= 1e7` keeps every count below `2^63`); `tau_b` is `sign sqrt((C - D)^2 / ((n0 - n1)(n0 - n2)))`, one rounding, by the same exact root. p-value by the normal approximation with the standard tie-corrected variance,
  `p_method: "normal"`; **[open]** the exact permutation distribution for small `n` without ties (not in the first slice, said in the answer).
* **Accuracy claims [choice].** `r`, `rho`, `tau_b`: correctly rounded (one rounding of an exact quantity). The p-values: the S0 bound of the function they call (`t` tail: 1e-12 relative), measured again here against `mpmath` on the EXACT rational inputs.

## 3. Exactness plan

* **Cross-product accumulator.** `psum.add_cross` lays `mx my 2^(fx + fy)` (the integer product of the two 53-bit mantissas, 106 bits) at the bit offset `fx + fy` in the same limb layout as the squares of `psum` (units of `2^-2148`), the sign
  by the signs of the two factors, so `Sxy` has exactly the unit of `n Sxy` terms of `Sx Sy` and the subtraction needs no rescaling. The squares `Sxx`, `Syy` are the `k = 2` accumulators of S1.
* **Per pair memory.** Pairwise needs `Sx`, `Sy`, `Sxx`, `Syy`, `Sxy` over the pair's rows: about 1,350 integers per pair; `--max-vars` (default 32, ceiling 128) bounds the pairs, and the memory is `vars (vars - 1) / 2` times that.
* **Order.** All sums are of integers: row order and thread partition change no bit; the gates shuffle the rows.
* **Time [not measured yet].** `O(pairs x rows)` accumulations for Pearson, `O(pairs x rows log rows)` for Spearman and Kendall; measured before any claim.

## 4. Refusals

`stats.non-numeric`, `stats.corr-too-few-columns`, `stats.no-such-column`, `stats.ambiguous-column` (as `--cols` of `describe`), `limit.too-many-columns` (more selected variables than `--max-vars`), `stats.bad-method` (a method not in
the list), the reader's and the limits' tags of `describe`.

## 5. Gates (written before the code, each able to fail)

| id | gate | fails when |
|---|---|---|
| C-A | the cross-product accumulator against `fractions.Fraction`: vectors of the families of `measure_acc.py` including wide exponents, subnormals and opposite signs, the accumulator's sum and `sum x y` rounded once | one bit differs |
| C-B | Pearson oracle: generated files (columns of the shapes of `gate_describe.py`, pairs of every combination, missing cells in each pattern, constants, perfect correlation, anticorrelation, huge and tiny scales) against `Fraction`: `r` bit for bit, `n`, `p` against `mpmath` at the exact rational `(x, y)` with the bound of S0 | one byte differs |
| C-C | Spearman against a Python implementation of average ranks and the same `Fraction` Pearson | one bit differs |
| C-D | Kendall against a brute-force `O(n^2)` count in Python and the `Fraction` root; p-value against the closed form | one bit differs |
| C-E | row order and column order: shuffled rows give the same bytes; the matrix is symmetric; permuting `--cols` permutes the answer | a byte differs |
| C-F | refusals and fuzz (the table and the hostile files of `describe`, plus `corr` arguments) | a trap, a hang, an untagged refusal |
| C-G | mutants of the accumulator, the exact formulas, the p-value pair, the missing-value rule, the ranks, the Knight count | an unexplained survivor, with the control of the unmutated build first |
| C-H | fingerprint lines for `corr`, equal on the platforms of the workflow, thread check | a line differs |

## 6. Open points for the maintainer

1. `--missing pairwise` as the default (pandas' and R's `pairwise.complete.obs`), against `listwise` (one set of rows, a positive semidefinite matrix).
2. The p-value of Pearson and Spearman by the t approximation (what `scipy.stats.pearsonr` and `cor.test` print) and of Kendall by the normal approximation, with the exact small-sample distributions left for a later slice.
3. `--max-vars` default 32.

## 7. Slices

* **S2a: Pearson** (this slice): the cross-product accumulator, the exact formulas, `corr --method pearson`, the p-value pair, `--missing`, `--cols`, `--format text`; gates C-A, C-B, C-E, C-F, C-G, C-H for Pearson.
* **S2b: Spearman**, **S2c: Kendall tau-b**: each its own slice with its gates (C-C, C-D).
