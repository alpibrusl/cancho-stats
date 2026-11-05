# `stats did` and `stats cycle`: a regression-adjusted difference, and a cycle found with no period given (milestone S6)

**Status: designed before any code of S6. `stats did` (S6.1) is built, section 12 says where the build differs from this design and what it measured; `stats cycle fit` (S6.2) and `stats cycle detect` (S6.3) are built, section 13.** Tags: **[measured]** printed by a script on the machine named in section 9; **[read]** read from a file of this repository or of the user (the file is named); **[choice]** a decision the maintainer can overturn; **[open]** not decided and not measured.

## 0. Where this comes from, and why it is here

`cancho-consulting/metasearch-bidding` is being moved from Python to cancho (`docs/pure-cancho.md` there). Two of its Python modules cannot move because `stats` lacks the estimation they use **[read]** `service/bidding/analytics.py` (`adjusted_did`, `season_adjusted`) and `service/bidding/seasonality.py` (`detect`, `detect_multi`, `segment_coefs`; 277 lines):

* a **regression-adjusted comparison** of treated and holdout units before and after a switch-on, with a bootstrap interval;
* a **search for one periodic component** in a daily series (a sine and a cosine on a grid of periods, over a polynomial level, a permutation test of the best one), shared by many series, and the fit of that cycle at a known period with a shrunk amplitude per key.

Neither is specific to hotels. The first is the estimator every before-and-after study with a holdout asks for; the second is what anyone asks of a daily series with a calendar in it. Both are the kind of answer `stats` gives: an exact sum where a sum decides, a seed where a draw decides, and a refusal with a rule tag where the data cannot support a number.

This document adds **three verbs**: `stats did`, `stats cycle detect`, `stats cycle fit`. It adds no dependency: least squares uses the exact accumulators already pinned (`facc`, `acc.psum`, `acc.central` **[read]** `docs/numerics.md`) and the Cholesky routine of the `linalg` store of cancho-ml **[read]** `lib/linalg/chol.cho` in cancho-ml; the draws use `ml.rng` as `resample` does **[read]** `docs/resample.md` section 2.

## 1. What the Python does, stated so that it can be judged against

### 1.1 `adjusted_did`

For units with at least `min_pre_days` (14) of listed pre-period days, each unit has a pre level `x = pre_sum / pre_days` and a post level `y = post_sum / post_days` (a unit's days are the days of the window on which it exists). With `treated` a 0/1 indicator, `y = a + b x + tau treated` is fitted by ordinary least squares. The effect is `tau / (mean(y over treated) - tau)` when that denominator is positive and not a number otherwise; the absolute effect is `tau`. A naive figure `sum(y treated) / sum(x treated) - 1` is reported beside it. The interval is the 2.5% and 97.5% percentiles of 1000 (`BOOT`) bootstrap replicates, each drawing the treated units with replacement and the holdout units with replacement (a **stratified** bootstrap over units), refitting, and recomputing both figures; a replicate with no treated or no holdout unit is not a number and is left out of the percentile of the relative effect.

### 1.2 `detect` and `detect_multi`

Series `(t, y, w)`: day numbers (integers, not contiguous), values, weights (inverse variance of a day). Days with a non-finite value or a non-positive weight are dropped; a series with fewer than `3 (order + 3)` days left is dropped. With `order = 2` the level under the cycle is the polynomial `1, t', t'^2` with `t' = (t - lo) / max(t - lo)` (this is `LEVEL_ORDER`). Over the series that remain:

* `span = (last day) - (first day) + 1`; if `span / 2 < 21` there is no answer (`MIN_PERIOD = 21`, `MIN_CYCLES = 2`).
* A grid of frequencies from `1 / (span / 2)` to `1 / 21`, in steps of `1 / (4 span)` (`OVERSAMPLE = 4`), i.e. periods from `span / 2` down to 21.
* For each period, each series' **gain** is the weighted sum of squares explained by adding `sin(2 pi t / p)` and `cos(2 pi t / p)` to its own level, minus what the level alone explains; the statistic of a period is the sum of the gains over series divided by the sum of the series' residual sums of squares after the level alone.
* The same statistic is computed on `499` (`PERMUTATIONS`) shuffles of the values **within every series** (weights and level stay in place); the statistic of a shuffle is its maximum over the grid, which is what makes the test one for "the best period", not for one chosen in advance. `p = (1 + number of shuffles whose maximum is at least the observed best) / 500`.
* The answer is **no cycle** when `p > 0.004` (`ALPHA`: at most one shuffle in 250), or when the best period exceeds `0.8` (`EDGE`) of the longest period tried (a cycle longer than the data resolve is fitted by a shorter, wrong one). Otherwise the best period is refined by a ternary search between its two neighbours on the grid (30 steps) and the answer is `period`, `amplitude` (the root mean square over series of the series' own `hypot(sin coefficient, cos coefficient)`), `p`, and `gain`.

### 1.3 `segment_coefs`

One series at a **known period**: fit level (a polynomial of order 2 in the series' own scaled time) and `a sin + b cos` of the absolute day by weighted least squares; the two cycle coefficients are shrunk towards zero by `max(0, 1 - 2 / F)` where `F` is the Wald chi-square of the two coefficients under the fitted residual variance (the positive-part James-Stein factor). A series with too few days returns `(0, 0)`.

### 1.4 What the three cost in numpy **[measured]**

On an Apple M4 Max, Python 3.14, numpy 2.4.6, one run each (not a benchmark): `detect` on a 270-day series with a cycle of 40, 55 or 90 days found periods 40.13, 54.55 and 91.81 with `p = 0.002`, in 0.01 s each; on 270 days of noise it returned no cycle in under 0.01 s; the grid had 44 periods; `detect_multi` over 20 series of 270 days found 90.11 days in 0.13 s. **Speed is therefore not what this design buys or has to defend.** What it buys is that the module leaves Python; what it has to defend is the number.

## 2. `stats did`

```
stats did --arm A --pre P --post Q [--pre-days D0 --post-days D1] [--min-pre-days 14] [--treated-value 1]
          [--replicates R --seed S] [--level 0.95] FILE
```

*(Corrected by the build, section 12: no `--unit`; `--replicates` and `--seed` replace `--boot`; `--treated-value` is 0 or 1.)* One row per unit. `A` is a column of 0 and 1 (1 = treated; **[choice]** the Python's `arm == 0` means treated and `arm == 1` the holdout; the verb takes `--treated-value V`, 0 or 1, so a client's own coding is not renamed). `P` and `Q` are the unit's pre and post sums; `D0` and `D1` the days they cover. When the day columns are absent the levels are `P` and `Q` as they are. The levels are `P / D0` and `Q / D1` **correctly rounded once** from decimal text, which is what `table --agg sum:` and a division would give **[choice]**.

The answer (`schema stats.v1`, `kind: "did"`): `n`, `n_treated`, `n_holdout`, `eligible` (units that pass `--min-pre-days`), the three coefficients `a`, `b`, `tau` with their standard errors, `effect` (`tau` over the treated units' counterfactual mean, absent with a reason when that mean is not positive), `naive`, and with `--boot` the percentile interval of both `effect` and `tau`, `replicates`, `replicates_defined`, and `mc_error` as `resample` reports it.

### 2.1 Arithmetic

`X'X` and `X'y` for the three-column design `[1, x, treated]` are sums of products of doubles: **[choice]** they are accumulated exactly (`acc.psum` for sums and squares, `acc.xsum` for the cross sum) and rounded once, so they do not depend on the order of the rows; the 3x3 system is then solved **by elimination with partial pivoting in doubles** *(corrected by the build: the design said the Cholesky routine of `linalg`; a 3x3 system does not need a new pinned store, and question 1 of section 11 is answered: no new store for `did`)*; a pivot below 1e-12 of the largest entry is `stats.did.singular`, section 2.3). The coefficients are therefore a function of the multiset of rows, as the means of `describe` are. This is stronger than numpy's `lstsq` (SVD), whose result depends on its own order of operations; it is **not claimed to equal it**. The bound of the oracle gate (G-D2) is fixed in section 7 before the code is written.

### 2.2 The bootstrap

The Python draws the treated units with replacement and the holdout units with replacement, separately. **[choice]** The verb keeps that, with the stream rule of `resample` **[read]** `docs/resample.md` section 2: replicate `b` draws from `derive(derive(seed, "did"), b)`, the treated first and the holdout after, so a replicate is a function of the data, the seed and `b` only, at any thread count. **The draws are not numpy's**, so a replicate is not the Python's replicate and the interval is not the Python's interval, to the last digit: the migration of `analytics` must accept this (section 8). What is kept is the statistic and the scheme.

### 2.3 Refusals **[choice]** (each with a `stats.` tag, an exit code and a repair)

`stats.did.no-treated`, `stats.did.no-holdout` (nothing to compare), `stats.did.too-few-units` (fewer than 4 eligible units, or fewer than 2 in an arm: the standard errors do not exist), `stats.did.singular` (all levels `x` equal, or `treated` collinear with the constant), `stats.did.arm-not-binary` (a value that is neither of the two), `stats.did.days-not-positive`, `stats.seed-required` (with `--boot`; the existing tag), and the table-wide ones the other verbs share. A unit with a missing `P`, `Q` or `A` is a `missing` count in the answer, as in `describe`, and is left out; **[open]** whether it should refuse instead.

## 3. `stats cycle detect`

```
stats cycle detect --day T --y Y [--weight W] --series K [--order 2] [--min-period 21] [--min-cycles 2]
                   [--oversample 4] [--shuffles 499] [--alpha 0.004] [--edge 0.8] [--grid yes] --seed S FILE
```

*(Corrected by the build, section 13: `--shuffles`, not `--permutations`, which `resample perm` already has with another default; `--series` is required in v1; `--grid yes` prints the gain of every period.)*

Long format: one row per (series, day). `--series` absent means one series. The defaults are the Python's constants (section 1.2) and every one is a flag, because they are the Python's **[choice]**s, not the statistics' **[read]**s; the answer reports the values used.

The answer: `found` (bool), `period`, `amplitude`, `p`, `gain`, and, **always**, the evidence a refusal needs: `grid` (number of periods, shortest, longest), `best_on_grid` (the period and gain the search found before the tests), `edge_refused` (bool), `series_used` and `series_dropped` with the count of days dropped for a non-finite value or a non-positive weight. A "no cycle" is `found: false` with `reason` (`too-short`, `p-above-alpha`, `edge`), **not** a refusal: for the consumer of this verb "nothing is there" is an ordinary answer and a nonzero exit would be the wrong interface. (The `stats test` verbs, which refuse on invalid input and answer on valid input, are the model.) Exit 0 for a found cycle and for none; refusals are for requests that cannot be answered (`stats.cycle.no-day-column`, `stats.cycle.day-not-integer`, `stats.cycle.duplicate-day` within a series, `stats.cycle.weight-negative`, `stats.cycle.permutations-range`, `stats.seed-required`).

### 3.1 Arithmetic and cost

For every series the projection of the weighted values on the level's column space is computed once (a QR of the `n x 3` level design in the Python; here the Cholesky of its normal equations in double-double, **[choice]** because `n x 3` with a scaled time has a condition number of the order of 10^2 **[open]** to be measured, and the squaring of it is harmless at that size), and for every period the projection on the level plus the sine and cosine (`n x 5`). The permutation statistic reuses the projections, as the Python does: `499 x 44 x 20` small projections for the 20-series case of section 1.4. **[open]** the cost in cancho; the budget is the numpy 0.13 s times a factor the gate sets (G-C5), because a cost that is not measured is not a claim.

### 3.2 Determinism

The shuffles are `derive(derive(seed, "cycle"), series)` streams (Fisher-Yates with the bounded draw of `ml.rng` **[read]**), one per series per permutation, so the p-value is a function of the data and the seed, at any thread count. The summed statistic is a sum of doubles in a fixed order; **[choice]** the order is series in order of first appearance, then period, and the answer's `gain` and `period` are bit-stable across the three platforms of the CI workflow (a line `cycle-detect` in `docs/xplat_fp.txt`), which is the same promise as every other verb.

**What is not promised:** that the p-value equals the Python's. The Python's shuffles are numpy's PCG64; these are SplitMix64 streams. Two seeds of either implementation give two p-values (the resolution is 0.002); a decision whose p lies within one resolution step of `alpha` can flip with the seed. G-C3 measures how often that happens on the cases of `tests/season_test.py` so that the consumer is told the number rather than a promise.

## 4. `stats cycle fit`

```
stats cycle fit --day T --y Y [--weight W] --series K --period P [--order 2] [--shrink james-stein|none] FILE
```

*(Corrected by the build, section 13: `--predict-days` is not built.)*

For each series (here the **key**: a market, a brand, a destination) the fit of section 1.3 at the given period, one row per key: `a`, `b` (shrunk), `a_raw`, `b_raw`, `f` (the Wald chi-square), `shrink`, `days`, `status` (`fitted` or `too-short`, whose coefficients are `0`). With `--predict-days` the same file's keys are evaluated on the days of the second file (level excluded by default: it is the *cycle* that the consumer adds to a log-odds, `cycle_at(day) = a sin(2 pi day / P) + b cos(2 pi day / P)`; **[choice]** `--with-level` adds the level for the one use that wants it, the season-adjusted before and after of `analytics.season_adjusted`).

The shrinkage is a **[choice]** of the Python (positive-part James-Stein on the two coefficients); the verb states it in its answer so that a reader can see it was applied, and `--shrink none` returns the raw fit. The degrees of freedom of the residual variance are `max(n - 5, 1)`, as the Python has it.

## 5. What the three verbs do not do **[choice]**

No time-series model (no ARIMA, no autocorrelation of the residuals: a daily series with serial correlation has a permutation test that is not valid, **and the detect verb does not check for it**; it says so in `limits`), no harmonics beyond the first, no more than one cycle (the Python looks for one), no missing-day interpolation, no calendar (a week is the consumer's `dow`; periods under `--min-period` are not looked for), no weights that depend on the fit. No automatic choice of `--order`.

**The permutation test assumes exchangeable days under the null, which a trend or a step is not, and the Python knows it:** the level stays in place and only the values move, which is a test on the *residual* shape only if the level is flat. The Python accepts the approximation (the polynomial level absorbs growth and a step at the switch-on, and the measured false-positive rate is in `tests/season_test.py`); the verb keeps the Python's test **[choice]** so that the migration changes the language and not the statistic, and G-C4 measures the false-positive rate of the verb on series with growth and on pure noise rather than assuming it.

## 6. Where it sits

* `lib/stest/did.cho`, `lib/stest/cycle.cho`, `lib/stest/cyclefit.cho` and their answer files, in the pattern of `ttest.cho` / `answer.cho` **[read]** `lib/stest`. No file over 2,000 lines; the answer builders are split as `answer2.cho` ... `answer10.cho` already are.
* One pinned addition: the `linalg` Cholesky store of cancho-ml (`docs/numerics.md` section 2 lists the stores; this is one more block in `cancho.toml` and one more line of `scripts/check_pins.sh`) **[read]**. If the maintainer prefers no new store, the 3x3 and 5x5 systems can be solved inside `lib/special` in double-double **[open]**.
* `stats introspect` gains the three verbs; the tables of rule tags, exit codes, repairs and the OpenAPI-like schema follow the others.
* Authority: unchanged (`Fs` read only, no network, no clock).

## 7. Gates, written before the code, each able to fail

Every gate has a mutant list; a surviving mutant is explained or the gate is rewritten **[read]** CLAUDE.md.

| id | judges | how | fails on |
|---|---|---|---|
| G-D1 | `did` by hand | a fixture of 8 units whose sums, `X'X`, `tau`, `effect` and `naive` are worked on paper in `docs/estimation.md` section 10 and compared byte for byte | any change in the definition |
| G-D2 | `did` against numpy | 500 seeded random tables (sizes 8 to 400, ties in `x`, a heavy unit, negative levels); `tau`, `a`, `b`, `effect` equal numpy's `lstsq` to a **relative bound of 1e-9** **[choice]**, fixed here, before the run; the maximum measured gap is then written in this section | a gap above the bound |
| G-D3 | `did` is order-free | the same table in 20 shuffled row orders: identical bytes | any dependence on order |
| G-D4 | the interval has coverage | simulated studies with a known `tau` (800 per setting, three settings of `n`, one with a heavy unit): the share of intervals containing `tau`, with its binomial standard error, written in the document; **no claim of exact 95%**, a claim of what was measured **[open]** the acceptance band | a coverage outside the band fixed in advance |
| G-D5 | the stratified draw | a mutant that draws from the pooled units must be killed by a table where the arms differ in size | the mutant surviving |
| G-D6 | refusals | each tag of section 2.3 reached by a table, with the exit code | a tag that cannot be reached or a panic (a fuzz of malformed tables) |
| G-C1 | `detect` against the Python | the Python `detect` is re-implemented in `scripts/oracle_cycle.py` from section 1.2 of this document, **not copied from the client's repository**; on 300 seeded series the grid, the statistic of each period and the best-on-grid period equal to a **relative 1e-9** **[choice]** | a difference |
| G-C2 | power | cycles of 40, 55, 90 and 120 days at three amplitudes, 270 days, 10 seeds each: the share found and the period error; the numbers go into the document | below the floors fixed in advance: 10 of 10 at the amplitude of `season_test` **[read]** `tests/season_test.py` in the client's repository |
| G-C3 | seed sensitivity | the decision (found / not found) of the same series under 20 seeds: the share of series whose decision is not the same under all of them | reported, no floor: this is what the consumer is told |
| G-C4 | false positives | 2000 series of noise, 2000 with growth, 2000 with a step; the rate of `found` against `alpha = 0.004`; the number goes into the document, and if it exceeds `3 alpha` the verb does not ship with these defaults | a rate above `3 alpha` |
| G-C5 | cost | the 20-series case of section 1.4, minimum of 5 runs, on a machine named in the document | a cost above the factor fixed in advance **[open]** |
| G-C6 | determinism | thread counts 1 and 8 and the three platforms: identical bytes | any difference |
| G-F1 | `fit` by hand | one series of 12 days with a worked solution | any change |
| G-F2 | `fit` against numpy | 500 series; `a_raw`, `b_raw`, `f` to the bound of G-D2 | a gap above the bound |
| G-F3 | the shrinkage | `F <= 2` gives `(0, 0)`; `F` large gives the raw fit; a mutant with the sign of `2 / F` changed is killed | the mutant surviving |
| G-F4 | too-short series | a series of fewer days than `order + 3`, and of exactly that many | a wrong status |

## 8. What the consumer must accept

A. The **interval of `did`** and the **p-value of `detect`** are not the Python's to the last digit (sections 2.2 and 3.2). The goldens of `metasearch-bidding` that carry them are re-recorded once, after G-D4 and G-C3 have given the numbers, and the commit says so.

B. The **point figures** (`tau`, `effect`, `naive`, `period`, `gain`, the coefficients) are expected to equal the Python's to the bound of G-D2 / G-C1, so they are compared with a tolerance and not byte for byte; the migration's own parity gate states the tolerance.

C. A series with serial correlation is outside what the test covers (section 5).

## 9. Milestones

| | what | done when |
|---|---|---|
| S6.0 | this document, reviewed | the maintainer's answers to section 11 are in it |
| S6.1 | `did`, G-D1 to G-D6, mutants | the gate passes, the mutants are killed |
| S6.2 | `cycle fit`, G-F1 to G-F4 | same |
| S6.3 | `cycle detect`, G-C1 to G-C6 | same |
| S6.4 | the fingerprint lines and the pages | `docs/xplat_fp.txt` equal on the three platforms |

Machine for every **[measured]** here: Apple M4 Max, macOS, Python 3.14, numpy 2.4.6; no cancho figure has been measured.

## 10. Hand-worked fixture for G-D1

Eight units, no day columns (so `x = pre`, `y = post`):

| unit | arm | pre | post |
|---|---|---|---|
| a | 1 | 10 | 16 |
| b | 1 | 12 | 19 |
| c | 1 | 8 | 13 |
| d | 1 | 15 | 22 |
| e | 0 | 11 | 13 |
| f | 0 | 9 | 12 |
| g | 0 | 14 | 15 |
| h | 0 | 10 | 12 |

By hand: `n = 8`, `n1 = 4`; `Sx = 89`, `Sy = 122`, `Sxx = 1031`, `Sxy = 1403`; over the treated units `Sxt = 45`, `Syt = 70`. The normal equations are

```
 8 a +   89 b + 4 tau = 122
89 a + 1031 b + 45 tau = 1403
 4 a +   45 b + 4 tau = 70
```

whose exact solution (fractions, eliminated on paper and checked by `scripts/gate_did.py`) is `a = 205/163`, `b = 174/163`, `tau = 690/163`. The treated units' mean `y` is `70/4 = 35/2`, the counterfactual `35/2 - 690/163 = 4325/326`, so `effect = tau / counterfactual = 276/865` and `naive = Syt / Sxt - 1 = 70/45 - 1 = 5/9`. The program prints `tau = 4.233128834355829`, `effect = 0.3190751445086705` (both the correctly rounded value of the exact fraction), `naive = 0.5555555555555556`, `a = 1.2576687116564766` (the exact value rounded is `1.2576687116564418`: a relative gap of 2.8e-14, inside the bound of G-D2) and `b = 1.067484662576684`.

## 11. Open questions for the maintainer

1. **A new pinned store (`linalg`)** or the small solves inside `lib/special`? (Section 6.)
2. **`--treated-value`**: is renaming the arm in the request the right interface, or should the verb require 0/1?
3. **A missing cell**: leave the unit out and count it (as `describe`), or refuse?
4. **`cycle detect` answering "none found" with exit 0** (section 3): agreed? A nonzero exit would be easier for a shell but is a false report of failure.
5. **The Python's test is kept** although the days of a series are not exchangeable under a trend (section 5). Keep, or add a documented alternative (shuffling blocks of days) behind a flag in S6.3?
6. **The name**: `cycle` (a cycle in a daily series) or `periodogram` (what an agent would search for)? `stats periodogram` is the name `docs/pure-cancho.md` of the client's repository proposed; the verb it describes is a *tested* single-period search, which a periodogram is not.

## 12. What the build of S6.1 (`stats did`) measured, and where it differs from sections 1 to 11

Machine for every figure here: Apple M4 Max, macOS, Python 3.14, numpy 2.4.6, the compiler pinned in `cancho.toml` (`de4e52fd`), one run each.

**Differences from the design** (each corrected in place above): no `--unit`; `--replicates` and `--seed` replace `--boot` (the bootstrap runs when either is given, and `--seed` is then required, `stats.seed-required`); `--treated-value` is 0 or 1; the system is solved by elimination with partial pivoting, not by the Cholesky of `linalg`, so no new pinned store (question 1 of section 11: answered, none for `did`); no standard errors and no `mc_error` in v1, which the client does not use; the stream is `derive(derive(seed, 3), b)` (kinds 1 and 2 are `ci boot` and `resample perm`); a unit with a missing cell is left out and counted in `missing` (question 3: counted, as `describe`); `--treated-value` and the refusal tags are those of `scripts/gate_did.py`. The refusal text of an unknown or non-numeric column still says `--cols`, because it is the shared `rule.refuse_pick`: **[open]** to make it name the flag.

**G-D1** the fixture of section 10 is compared in `scripts/gate_did.py` against exact fractions. **G-D2** 268 random tables (6 to 200 units, day columns or not, negative levels, a heavy unit, ties, a missing cell in a third of them): worst relative gap to the exact solution `a` 3.8e-12, `b` 3.1e-14, `tau` 1.0e-11, `effect` 1.0e-11, `naive` 5.1e-13, against the bound 1e-9 fixed in section 7 before the run; numpy's `lstsq` is within 1.2e-11 of the program on `tau`. **G-D3** the rows in the reverse order give the same bytes in every case; the arm coded the other way round with `--treated-value 0` gives the same bytes. **G-D5** 29 bootstrap tables: every replicate recomputed in exact fractions from the draw rule of section 2.2, the interval endpoints equal within the bound and the percentile positions and the number of defined replicates exactly equal; the three controls (a pooled draw, the stream kind 2, wrong percentile positions) each gave a different interval from the program's on most cases (28, 29 and 29 of 29), so the gate can tell them apart. **G-D6** every refusal of section 2.3 is reached with its tag. `scripts/mutants_did.py`: 41 mutants, 40 killed, 1 survived and explained (`units-total`, equivalent). The first run had 3 survivors, two of them gaps in the gate (a table with one treated unit and four or more units in all; the ceiling of 100000 replicates), which the gate now covers.

**Not measured, so not claimed:** **G-D4** the coverage of the interval (the share of simulated studies whose interval holds the true `tau`) has no script yet; nothing in this document says the interval has 95% coverage. The six `did` lines of `docs/xplat_fp.txt` were produced on this machine only; the Linux runners of the CI matrix have not seen them. The cost of a run was not timed.

## 13. What the build of S6.2 (`cycle fit`) and S6.3 (`cycle detect`) measured, and where it differs

Machine: Apple M4 Max, macOS, Python 3.14, numpy 2.4.6, compiler `de4e52fd`, one run each unless a count says otherwise.

**Differences from the design.** `cycle fit`: no `--predict-days` yet; the least squares is a Householder QR of the weighted design in doubles (not exact sums: the weights and the sine and cosine are not exact quantities), the rows are sorted by day first so the bytes do not depend on the order of the rows; a design whose QR has a diagonal below 1e-6 of the largest is `status: "degenerate"` with zeros, **a difference from the Python, which prints whatever `lstsq` returns**: the first version of the gate found a case (a cycle coefficient of -4.6e6 on data of order 1) where the QR could not honour the bound of 1e-9, and a verb that prints nine digits it cannot stand behind is worse than one that says so. The same holds of `lstsq` in numpy; the Python simply does not say. 4 of 316 series of the gate were reported degenerate, and the gate fails if more than a fifth are. `cycle detect`: `--shuffles` replaces `--permutations`; `--series` is required (a single series needs a constant column); the shuffle of series `s` (in order of the series keys ascending by bytes, among the series long enough) replicate `m` is Fisher-Yates from the end with the stream `derive(derive(derive(seed, 4), s), m)`, draw `j` taking `r = (z >> 1) mod (i + 1)` for `i = n - 1 - j`; the grid, the statistic, the decision and the refinement are section 1.2 exactly; a grid of more than 20,000,000 (periods x shuffles) is `limit.too-much-work`. "Nothing found" is exit 0 with `found: false` and a `reason` (`too-short`, `p-above-alpha`, `edge`), as section 3 proposed (question 4).

**G-F2** (`scripts/gate_cyclefit.py`, 120 files, 316 series, orders 0 to 3, periods 7.3 to 365, weights, gaps, a zero weight): every number against the exact solution of the same weighted design in fractions; worst relative gap `a_raw` 1.9e-10, `b_raw` 2.4e-10, `f` 6.4e-11, `shrink` 2.7e-11, against the bound 1e-9 of section 7 (fixed before the run; the rank threshold above was changed after the gate showed a case that violated it, and that is the whole of the change). The rows reversed give the same bytes; a series of fewer than `3 (order + 3)` days is `too-short` and exactly that many is fitted, at three orders; three controls (the shrinkage sign, the level one degree off, the divisor `n` for `n - p`) differ on 191, 95 and 107 of the series. Mutants (`scripts/mutants_cyclefit.py`): **24 of 24 killed**.

**G-C1/G-C2/G-C6** (`scripts/gate_cycledetect.py`, against `scripts/oracle_cycle.py`, numpy with the shuffles of this document): 30 seeded multi-series files (1 to 8 series, 12 of the 30 with a cycle found), the p-value **exactly** equal (it is a count), the gain of every period of the grid to a relative 1.2e-14, the refined period to 3.6e-16 and the amplitude to 4.0e-15; the rows reversed give the same bytes; three controls (another stream kind, one shuffle fewer, a level of another degree) differ on 7, 14 and 14 of 14 cases (the first run of the gate had controls that did not, because almost every case had a strong cycle and a p-value at its floor; the amplitudes were lowered so that p varies).

**G-C2 power** (`scripts/measure_cycle.py`, 270 days, noise uniform in +-0.6, one series, 10 seeds each): a cycle of 40 days is found 0, 10 and 10 of 10 at amplitudes 0.1, 0.2 and 0.4; of 55 days 2, 10 and 10; of 90 days 5, 9 and 10; **a cycle of 120 days is never found (0 of 10 at every amplitude)**: it is above `0.8 x` the longest period tried (135 days for a span of 270), which is the Python's edge rule (a cycle longer than the data resolve is fitted by a shorter, wrong one), not a lack of power. The largest error of a found period was 11.6% (90 days, amplitude 0.1) and 2.4% or less at amplitude 0.4. **G-C3 seed sensitivity**: 6 of 30 weak-cycle series (amplitude 0.1, 90 days) change decision with the seed among 20 seeds. **G-C5 cost**: 20 series of 270 days, **0.47 s** (the numpy original: 0.13 s on this machine, section 1.4).

**G-C4 false positives, against the rule of section 7 ("if it exceeds `3 alpha` the verb does not ship with these defaults").** 2000 series of 270 days each, one series, `alpha = 0.004`: pure noise 3 found (0.15%, standard error 0.09 points); noise with growth 4 (0.20%, 0.10); **noise with a step at day 135, 38 found (1.90%, standard error 0.31 points), 4.7 times `alpha` and above `3 alpha` = 1.2% by more than two standard errors.** The polynomial level of order 2 does not absorb a step, the days are not exchangeable around it, and the permutation test of section 1.2 does not know. This is a property of the Python's test, which this verb keeps for parity (section 5, question 5). **By the rule of section 7 the verb does not satisfy G-C4 with these defaults**; it is built and documented, and whether it ships with the Python's test, with a block shuffle behind a flag, or with `--order 3`, is the maintainer's decision (question 5). The client's own runs add a step at the switch-on of the service, which is exactly this case.

**Not measured:** the coverage of the interval of `did` (G-D4); the cross-platform bytes of `cycle` on the Linux runners (the cycle lines of `docs/xplat_fp.txt` are from this Mac; the sine and cosine are those of `std.math`, fdlibm in cancho, not the platform's, so no difference is expected, which is a reason and not a measurement); `cycle detect` with more than 20 series or series of a year or more; the cross-platform bytes are the only claim of section 7 that this repository cannot yet check alone.

**Mutants of `cycle detect`** (`scripts/mutants_cycledetect.py`, 34 mutants of the shuffles, the grid, the statistic, the p-value, the decision, the refinement and the flags): **33 killed, 1 survived and explained** (`shift-days`: equivalent, the answer does not depend on where the days start). The first run killed 29 of 33; the four survivors were two gaps in the gate (a series of exactly `3 (order + 3)` days was never in a case, so `>=` against `>` for keeping a series went unseen; one case in a hundred changes the number of periods on the grid, which 10 cases did not reach, now 40 and the grid start in `detect_core` is a mutant of its own), one equivalent mutant, and the weakness of the 10-case gate itself.
