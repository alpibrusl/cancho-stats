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
`variance`, `sd`, `skewness`, `kurtosis` (excess), `quantiles`, `distinct` and a `histogram`, one JSON document of the contract of cancho-tools' `toolbox` (`stats.v1`, `schemas/stats.v1.json`), with `introspect` and `skill` from the same flag table and the authority report the compiler derives
(`manifests/stats.authority.json`, bounded; `args`, `dir_read`, `err_write`, `file_read`, `fs_read`, `heap`, `io_write`; no `conc`, `clock`, `net` or `ffi`).

**Where the code differs from sections 2 and 3.**
* **`S1` and `S2` both come from `acc.psum`** (k = 1 and 2, own accumulators, gated in `scripts/gate_psum.py`), not from `facc` and `ml.fsq`: the numerator `n S2 - S1^2` needs both sums as limb numbers in
  compatible units (2^-2148), and `ml.fsq` and `facc` publish only a rounded result. `facc` and `ml.fsq` stay pinned for the reader and the measurement of `docs/numerics.md` 8.2.
* **Skewness and excess kurtosis are exact rationals rounded once** (`lib/acc/central.cho`): with `L2 = n C2`, `K3 = n^2 C3`, `K4 = n^3 C4` computed as exact integer combinations of the power sums (acc.psum k = 1 to 4),
  kurtosis is `(K4 - 3 L2^2) / L2^2` (one rounding) and skewness is `sign(K3) sqrt(K3^2 / L2^3)` through an exact integer square root with a sticky bit (one rounding); the moments are the biased ones of section 2. Both are `null` when the variance is zero or `n` is 1.
* **Quantiles are exact** (`lib/acc/quantile.cho`, `--quantiles`, `--quantile-type`): the nine Hyndman-Fan definitions on the sorted numeric cells (an LSD radix sort of order-preserving integer keys, the pinned `metrics.rank`), with the
  probability an exact decimal rational of at most 9 decimals; types 1 to 3 pick an order statistic, types 4 to 9 are `(1 - h) x[j] + h x[j + 1]` evaluated as one exact rational over `24 10^k` and rounded once. The default is type 7
  and the probabilities 0, 0.25, 0.5, 0.75, 1. The cells of the columns are held in memory (`--max-cells`, default 20,000,000, ceiling 100,000,000, one more is `limit.too-many-cells`); `p` is echoed as the text given.
* **Distinct counts and the histogram of a numeric column are exact** (`lib/acc/histogram.cho`, `--bins`, default 10, 0 for none, ceiling 1000): `distinct` is the number of different values (`-0` and `0` are one);
  the histogram has `bins` equal-width bins over `[min, max]`, `edges` that are `((N - i) min + i max) / N` as exact rationals, each printed as the nearest double (one rounding), and `counts` found by an exact comparison of every value
  with the RATIONAL edge (a signed combination of three doubles with small integer coefficients, formed exactly), the last bin closed. A column whose values are all equal has one bin, `edges` [min, max] and the count `n`.
  The levels of a TEXT column are in section 9.
* **`sd` is the square root of the rounded variance**: two roundings, said in the answer (`rounding`), the simpler of the two options of section 2 and the one of open point 2; the tie check of section 2 is not built.
* **`ddof` is 0 or 1**; the variance is `null` when `n - ddof < 1` or when the exact variance is beyond the largest double.
* **A cell that is not a number makes its column `text`**: a quoted cell, an empty-or-not cell the reader of cancho-table (`flt.parse_float`) refuses (NaN, infinity, a bare sign, a space). No cell value is printed.
* **A blank line is skipped, except in a file of one column, where it is one missing cell.** A record that continues on the next line (a quoted newline) is refused (`stats.csv-multiline`), not read.
* Rules of this tool's own: `stats.csv-ragged`, `stats.csv-multiline`, `stats.csv-bad-quote`, `stats.no-header`, `stats.bad-verb` (exit 2), `limit.too-many-rows`, `limit.too-many-columns`, each with exit 8 except the verb.

**Gates of this slice** (all in CI and `scripts/check.sh`): `scripts/gate_describe.py` (generated files of 14 column shapes, every number (including skewness and kurtosis) against `Fraction`, bit for bit, with `--ddof 1` and `0` and a shuffled-row rerun that must print
the same bytes), `scripts/gate_describe_refusals.py` (25 table cases with the exact tag and exit code, and fuzz of hostile files and arguments), `scripts/mutants_stats.py` (26 mutants), `scripts/manifest.py --check`, and the
fingerprint lines `describe` of `docs/xplat_fp.txt`. Measured results are in section 8.

**Not built yet:** `--cols`, `--format text`, the tie check of `sd`, records that continue on the next line. The four open points of section 6 keep their proposed defaults and are the maintainer's to change.

## 8. Results of the first slice

*(printed by the scripts named; the machine is the Mac of docs/numerics.md section 8, load as stated there)*

* **Against exact fractions** (`gate_describe.py --files 1000`, Mac): 1,000 files, 64,394 values (count, missing, kind, min, max, mean, variance, sd, skewness, excess kurtosis, the quantiles, the distinct count and the histogram (edges and counts), default and random type, probabilities and bins, of columns of the shapes uniform, integers, an offset of 1e8 with noise
  below 1, wide exponents (2^-300 to 2^300), tiny (about 2^-1000), huge (2^600 to 2^1000), constant, gappy, empty, text, mixed, quoted, exponent notation), compared bit for bit with `Fraction` rounded once: **0 differ**;
  every file also rerun with its rows shuffled and with `--ddof 0`: the shuffled bytes are equal in every file. The variance of a column of huge values is `null` exactly when the exact variance is beyond the largest double.
  On gram (Linux x86-64, `taskset -c 0-6`): 400 files, 26,173 values, 0 differ. The skewness and kurtosis columns were added after the first measurement of this slice (which had 40,078 values); the numbers above are for the code with them.
* **Refusals and fuzz** (`gate_describe_refusals.py`): 39 table cases, each with the exact rule tag and exit code, and 5,000 fuzz files (Mac; 2,000 on gram): **0 problems**: no trap, no hang, one line of valid JSON, every `error.rule` one
  of the tags of `introspect`, every exit code that of its rule.
* **Mutants** (`mutants_stats.py`): **26 of 26 killed** (the sign and the unit of the mean, the unit of the variance, the divisor, the numerator's subtraction, the limb shift, the comparisons of min and max, the first value, the counts,
  the quoted-cell rule, the square root, the guards, the power of the second accumulator, the blank-line rule, the ragged check, the limit edges, the empty kind, the verb, the delimiter, the `ddof` ceiling, the key of a negative number,
  the multi-line record, the overflow code). `mutants_central.py` (skewness, kurtosis and the exact square root; 29 mutants: the factors and signs of the three integer combinations, the subtraction of `3 L2^2`, the zero-variance guard, the sign and the zero case of the skewness, the number of bits of the square root and of the quotient, the equality case of the root, the sticky bit, the exponent, the shift, the accumulators of the third and fourth powers): **29 of 29 killed**. `mutants_histogram.py` (the signs and coefficients of the exact comparison and of the edges, the equality case, the unit, the search, the alignment of exponents, the constant column, the distinct counter, the differences of the counts, the last edge, the default and the ceiling of `--bins`; 14 mutants): **14 of 14 killed**. `mutants_quantile.py` (the nine definitions' constants and branches, the clamps, the exact blend, the sentinel and gather of the cells, the sort range, the probability grammar, the cells limit and the default type; 28 mutants): **28 of 28 killed**. `mutants_psum.py` and `mutants_bigint.py` (S1's accumulators): 19 of 20 and 20 of 21 killed, the survivors explained in the files.

  **Correction (2026-10-16).** The four counts above for `mutants_stats.py`, `mutants_central.py`, `mutants_histogram.py` and `mutants_quantile.py` (26, 29, 14 and 28 "killed") were not measurements: those scripts passed `--driver` to `gate_describe.py`, which only knows `--bin`, so the gate exited non-zero on its arguments for every mutant and every one counted as killed. Found while adding the levels mutants; fixed, and each script now first runs the same gates on the UNMUTATED build and refuses to count kills if they do not pass. Two patterns had also gone stale after the refactors (`count-limb-shift`, `second-power-accumulator`) and were moved or rewritten. The real counts, run on the Mac on 2026-10-16: `mutants_stats.py` **25 of 25**; `mutants_histogram.py` **14 of 14**; `mutants_quantile.py` **26 of 28**, the two survivors explained in the file (`blend-equal-keys`, `type3-negative-guard`: equivalent as far as the gates can show); `mutants_central.py` **27 of 30**, the three survivors explained in the file (`skewness-zero-numerator` equivalent; `sqrt-equality` and `sqrt-sticky` matter only at a rounding tie of the square root that no vector reaches: argued, not tested). The first real runs had 8 survivors; they were closed by gates that were missing (columns with the exact skewness -1.5 and with a zero third moment, columns exactly at `--max-cols`, a one-column file with blank lines, and `gate_quantile.py` added to the quantile mutants, which therefore also rebuild `quantile-check`).
* **The same bytes** (`xplat_fp.py`, 12 `describe` files added to the golden, now 60 lines): equal to the golden, written on the Mac, on gram, and the thread check of the other 48 lines holds on both; the fingerprint does
  not exercise threads for `describe` because the program has no `conc`: its authority has none, and `manifests/stats.authority.json` derived on gram is the one committed from the Mac (`manifest.py --check`).
* **On GitHub** (run 38811179034, `ubuntu-24.04`, `ubuntu-24.04-arm`, `macos-26`): the `describe` gates (200 files against Fractions, 1,000 fuzz files, `manifest.py --check`) and the fingerprint of 60 lines with its thread check
  all succeeded, none skipped.
* **Quantiles against an independent reference** (`gate_quantile.py`, the nine definitions written from the paper's formulas in `Fraction`, the cancho side in the `(a, b)` form so that the two share no algebra): 6,000 cases (Mac and gram) over all nine types, data of the families of `measure_acc.py` with ties and duplicates, 1 to 400 values, `p` from 0 to 1 including the neighbours of the jumps: **0 differ**, bit for bit. Reported only, not gated: 1,303 of 1,500 answers equal `numpy.quantile(method=...)` bit for bit (numpy rounds in floats).
* **Not measured:** speed (no time is asserted); files above 1,000 rows were not generated beyond 400 rows per column in these gates (`--max-rows` defaults to 1,000,000: the memory of the accumulators is 508 integers per column and
  does not grow with rows, but that bound is argued here, not measured at a million rows).

## 9. Levels of a text column (design, written before the code of this slice)

A column that has a cell that is not a number is `text` (section 7). This slice gives it `distinct` and its most frequent `levels`. Tags as above.

* **Two passes [choice].** Whether a column is text is known only at the end of the file (one bad cell makes it text), so the levels are counted in a second pass over the same file, reopened, only when some column is text and `--levels` is not 0. A single pass would count the levels of every numeric column too and spend the bound below on numbers.
* **What a level is [choice].** The bytes of the cell. A quoted cell is its bytes between the quotes (`"Paris"` and `Paris` are one level); a doubled quote inside stays doubled, it is not unescaped (the reader of cancho-table keeps it so). An empty unquoted cell is missing, not a level; a quoted empty cell `""` is the level of zero bytes. Bytes are compared as bytes: no case folding, no trimming, no Unicode normalisation.
* **Order [choice].** The `levels` are the first `--levels N` by count descending, then bytes ascending: a total order, so a tie is broken the same way on every machine. A value that is not UTF-8 is printed as `{"b64": ...}` (the `text.put` of the pinned `toolbox.text`).
* **Bounds.** `--levels N` (default 10, 0 for none, ceiling 100) and `--max-levels N` (default 100,000, ceiling 1,000,000): the most distinct (column, level) pairs held in all text columns together. A pair that would be the one more is not added and its column is `levels_capped: true` with `distinct: null` and `levels: null`: an answer that would be wrong (a capped count is not the count) is not printed, and the cap is data, not a refusal, because the numeric part of the answer is still right. The pairs are added in file order, so which column is capped does not depend on anything but the file and the flags.
* **Exact.** The counts are integers; `distinct` is the number of different byte strings. No rounding. The same bytes at any platform: no float is involved.
* **The file changed between the passes.** If the second pass sees a different number of rows, a record the first pass accepted that it refuses, or a line the first pass read that is now too long, the answer is `stats.file-changed` (exit 8), not a wrong count.
* **Memory [read, to be measured].** The keys of the held pairs, at most `--max-levels` of them, each at most `--max-line-bytes` long, plus 4 integers per pair; the selection of the top `N` is `N` passes over the pairs, at most 100,000,000 steps at both ceilings.
* **Gates.** `gate_describe.py` compares `distinct` and `levels` of generated text columns (many and few levels, quoted cells, UTF-8, ties, missing cells, a column of all-unique levels) with a Python count; the cap is checked on a file with one more pair than the bound; `gate_describe_refusals.py` gets the flag edges and the `stats.file-changed` case; mutants of the key, the order, the cap, the missing rule and the sort.

**What this slice built (2026-10-16).** `--levels N` and `--max-levels N` as designed; a text column prints `distinct`, `levels` (`value`, `count`) and `levels_capped`; `docs/example/readings.csv` has a `site` column to show it; `scripts/update_replay.py` rewrites the replay of the landing page from real runs and `--check` (in `check.sh` and CI) fails if the page differs from the binary. The pass over the file is the same reader with the same blank-line rule; the key of a pair is two bytes of column number then the cell, so up to 65,536 columns are told apart and the gate has a file of 300 text columns for it.

**Gates of this slice.** `gate_describe.py`: text shapes (`levels`, `qtext`, `utf8`, `uniq`) compared with a Python count for the default, for `--levels` 0, 1, 2, 3, 10 and 100, and for `--max-levels` 1 to 1000 against a simulation of the order in which pairs are added; fixed columns; a 300-column file; a one-column file with blank lines. `gate_describe_refusals.py`: the flag edges. `mutants_levels.py`: **22 of 24 killed**; the two survivors (`changed-ignored-rows`, `changed-ignored-ragged`) are in `stats.file-changed`, which is reached only if the file changes between the two passes: no gate can do that deterministically, so that rule is argued here and never exercised. The fingerprint golden `docs/xplat_fp.txt` was rewritten on purpose: the `describe` lines changed (new column shapes, and text columns now print levels); the 48 lines of the other sweeps are unchanged.

## 10. `--by COL` (design, written before the code of this slice)

`stats describe FILE --by COL` prints the description of section 1 once per value of the column `COL`. Tags as above.

* **What a group is [choice].** The rows that have the same bytes in column `COL` (the bytes of the cell, as a level in section 9: a quoted cell is its bytes between the quotes; no folding or trimming). A row whose `COL` cell is empty (missing) forms the group `by: null`, printed last; the other groups are ordered by the bytes of their value, ascending (a total order, the same on every machine). `COL` is a header name and must match exactly one column: no match is `stats.no-such-column`, two or more is `stats.ambiguous-column`, both exit 2. The `by` column is described like the others, inside every group (one value in each group, so its `distinct` is 1).
* **Exact by construction [choice].** One group at a time: the description of a group is the description of the file reduced to the rows of that group, computed by the code of section 7 on exactly those rows, so every number is the number of describing that sub-file by itself, bit for bit. No per-group accumulator is kept alongside the others; memory is the memory of one description, whatever the number of groups.
* **Cost [read].** A pass 0 finds the groups (the distinct values of `COL`, bounded), then one pass of the file per group, each reading every line and keeping the rows of its group: `G + 1` reads of the file for `G` groups. Time is `O(G * file)`, so `--max-groups` (default 50, ceiling 1000) bounds it; one more group is `limit.too-many-groups` (exit 8), never a silent drop. The text-levels second pass of section 9 runs inside each group when a column of that group is text.
* **A column's kind is per group [choice].** A column that is numeric in one group and has a text cell in another is `numeric` in the first and `text` in the second: the kind of a column is a property of the rows described. The answer says so (`kind` is inside each group's column).
* **Output.** `{"by": "COL", "groups": [{"value": ..., "rows": n, "columns": [...]}, ...]}` replaces the top-level `rows` and `columns` (`ddof`, `quantile_type` and `rounding` stay at the top). `value` is printed as the levels are (`text.put`). Without `--by` the document is unchanged.
* **Refusals.** `stats.no-such-column`, `stats.ambiguous-column`, `limit.too-many-groups`, plus the existing ones; `stats.file-changed` if a pass sees a different file than pass 0 (argued, not tested, as in section 9).
* **Gates.** The oracle `gate_describe.py` splits each generated file by a generated group column (few and many groups, missing group cells, quoted values, text values with bytes that sort differently as bytes and as Python strings) and compares every group with `expected()` of the rows of that group; groups reordered and rows shuffled give the same bytes; refusals for the three rules and the edges of `--max-groups`; mutants of the filter, the order, the null group, the key and the limit.

**What this slice built (`--by`).** `--by COL` and `--max-groups N` as designed in this section. Built as designed: pass 0 (`discover`) finds the groups, then one filtered pass per group runs the code of section 7 and section 9 on the rows of that group (`take_row` and `level_pass` skip a record whose `by` cell is not the group's; a record that is malformed is refused in every pass). Two differences from the text above: (1) the row and cell limits (`--max-rows`, `--max-cells`) apply to the rows of the group being described, not to the file, so a file may have more rows than `--max-rows` if no group does, while pass 0 still refuses a file with more than `--max-rows` records in all; (2) a record that is malformed AND beyond a limit now reports the malformation first (the limit is checked after the record is parsed, because whether a record counts depends on its group). The introspect list of rules and limits lacked `stats.file-changed` (section 9) and the flags `--levels`, `--max-levels`: added here, together with the three rules of `--by`.

**Gates and results (`--by`).** `gate_describe.py`: for each file an extra column `grp` at a random position (2, 3, 5 or 10 values among: ASCII, a quoted value, a quoted empty value, a value with a space, a non-ASCII value, a numeric-looking one, missing), the order of the groups, `rows` and every key of every column of every group against `expected()` of the rows of that group by themselves, and a shuffled-rows rerun that must give the same bytes. `gate_describe_refusals.py`: 16 table cases (no such column, case sensitive, two columns of one name, malformed file, empty file, the row, group and cell limits and their edges, the default of 50 groups and 51). `mutants_by.py` (the group test, the null group, the filter in both passes, the limits and their edges, the order, the column match; its first run left one survivor, the default of `--max-groups`, closed by the 50 and 51 group cases): **21 of 21 killed** (Mac, 2026-10-16). The fingerprint golden has 4 new `describe-by` lines (60 to 64); the other 60 lines are unchanged.
