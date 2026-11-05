# `stats did` and `stats cycle`: a regression-adjusted difference, and a cycle found with no period given (milestone S6)

**Status: designed before any code of S6. Nothing here is built.** Tags: **[measured]** printed by a script on the machine named in section 9; **[read]** read from a file of this repository or of the user (the file is named); **[choice]** a decision the maintainer can overturn; **[open]** not decided and not measured.

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
stats did --unit U --arm A --pre P --post Q [--pre-days D0 --post-days D1] [--min-pre-days 14]
          [--boot R --seed S] [--level 0.95] FILE
```

One row per unit. `A` is a column of 0 and 1 (1 = treated; **[choice]** the Python's `arm == 0` means treated and `arm == 1` the holdout; the verb takes `--treated-value V` so a client's own coding is not renamed). `P` and `Q` are the unit's pre and post sums; `D0` and `D1` the days they cover. When the day columns are absent the levels are `P` and `Q` as they are. The levels are `P / D0` and `Q / D1` **correctly rounded once** from decimal text, which is what `table --agg sum:` and a division would give **[choice]**.

The answer (`schema stats.v1`, `kind: "did"`): `n`, `n_treated`, `n_holdout`, `eligible` (units that pass `--min-pre-days`), the three coefficients `a`, `b`, `tau` with their standard errors, `effect` (`tau` over the treated units' counterfactual mean, absent with a reason when that mean is not positive), `naive`, and with `--boot` the percentile interval of both `effect` and `tau`, `replicates`, `replicates_defined`, and `mc_error` as `resample` reports it.

### 2.1 Arithmetic

`X'X` and `X'y` for the three-column design `[1, x, treated]` are sums of products of doubles: **[choice]** they are accumulated exactly (`acc.psum`) and rounded once, so they do not depend on the order of the rows; the 3x3 system is then solved with the Cholesky routine of `linalg` **[read]** (positive definite whenever the design has full rank; a singular design is refused, section 2.3). The coefficients are therefore a function of the multiset of rows, as the means of `describe` are. This is stronger than numpy's `lstsq` (SVD), whose result depends on its own order of operations; it is **not claimed to equal it**. The bound of the oracle gate (G-D2) is fixed in section 7 before the code is written.

### 2.2 The bootstrap

The Python draws the treated units with replacement and the holdout units with replacement, separately. **[choice]** The verb keeps that, with the stream rule of `resample` **[read]** `docs/resample.md` section 2: replicate `b` draws from `derive(derive(seed, "did"), b)`, the treated first and the holdout after, so a replicate is a function of the data, the seed and `b` only, at any thread count. **The draws are not numpy's**, so a replicate is not the Python's replicate and the interval is not the Python's interval, to the last digit: the migration of `analytics` must accept this (section 8). What is kept is the statistic and the scheme.

### 2.3 Refusals **[choice]** (each with a `stats.` tag, an exit code and a repair)

`stats.did.no-treated`, `stats.did.no-holdout` (nothing to compare), `stats.did.too-few-units` (fewer than 4 eligible units, or fewer than 2 in an arm: the standard errors do not exist), `stats.did.singular` (all levels `x` equal, or `treated` collinear with the constant), `stats.did.arm-not-binary` (a value that is neither of the two), `stats.did.days-not-positive`, `stats.seed-required` (with `--boot`; the existing tag), and the table-wide ones the other verbs share. A unit with a missing `P`, `Q` or `A` is a `missing` count in the answer, as in `describe`, and is left out; **[open]** whether it should refuse instead.

## 3. `stats cycle detect`

```
stats cycle detect --day T --y Y [--weight W] [--series K] [--order 2] [--min-period 21] [--min-cycles 2]
                   [--oversample 4] [--permutations 499] [--alpha 0.004] [--edge 0.8] --seed S FILE
```

Long format: one row per (series, day). `--series` absent means one series. The defaults are the Python's constants (section 1.2) and every one is a flag, because they are the Python's **[choice]**s, not the statistics' **[read]**s; the answer reports the values used.

The answer: `found` (bool), `period`, `amplitude`, `p`, `gain`, and, **always**, the evidence a refusal needs: `grid` (number of periods, shortest, longest), `best_on_grid` (the period and gain the search found before the tests), `edge_refused` (bool), `series_used` and `series_dropped` with the count of days dropped for a non-finite value or a non-positive weight. A "no cycle" is `found: false` with `reason` (`too-short`, `p-above-alpha`, `edge`), **not** a refusal: for the consumer of this verb "nothing is there" is an ordinary answer and a nonzero exit would be the wrong interface. (The `stats test` verbs, which refuse on invalid input and answer on valid input, are the model.) Exit 0 for a found cycle and for none; refusals are for requests that cannot be answered (`stats.cycle.no-day-column`, `stats.cycle.day-not-integer`, `stats.cycle.duplicate-day` within a series, `stats.cycle.weight-negative`, `stats.cycle.permutations-range`, `stats.seed-required`).

### 3.1 Arithmetic and cost

For every series the projection of the weighted values on the level's column space is computed once (a QR of the `n x 3` level design in the Python; here the Cholesky of its normal equations in double-double, **[choice]** because `n x 3` with a scaled time has a condition number of the order of 10^2 **[open]** to be measured, and the squaring of it is harmless at that size), and for every period the projection on the level plus the sine and cosine (`n x 5`). The permutation statistic reuses the projections, as the Python does: `499 x 44 x 20` small projections for the 20-series case of section 1.4. **[open]** the cost in cancho; the budget is the numpy 0.13 s times a factor the gate sets (G-C5), because a cost that is not measured is not a claim.

### 3.2 Determinism

The shuffles are `derive(derive(seed, "cycle"), series)` streams (Fisher-Yates with the bounded draw of `ml.rng` **[read]**), one per series per permutation, so the p-value is a function of the data and the seed, at any thread count. The summed statistic is a sum of doubles in a fixed order; **[choice]** the order is series in order of first appearance, then period, and the answer's `gain` and `period` are bit-stable across the three platforms of the CI workflow (a line `cycle-detect` in `docs/xplat_fp.txt`), which is the same promise as every other verb.

**What is not promised:** that the p-value equals the Python's. The Python's shuffles are numpy's PCG64; these are SplitMix64 streams. Two seeds of either implementation give two p-values (the resolution is 0.002); a decision whose p lies within one resolution step of `alpha` can flip with the seed. G-C3 measures how often that happens on the cases of `tests/season_test.py` so that the consumer is told the number rather than a promise.

## 4. `stats cycle fit`

```
stats cycle fit --day T --y Y [--weight W] --series K --period P [--order 2] [--shrink james-stein|none]
                [--predict-days FILE2 --predict-day T2] FILE
```

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

*To be written with S6.1 and before its code: eight units, their `x`, `y`, `treated`, the normal equations, `tau`, `effect`, `naive`.* **[open]**

## 11. Open questions for the maintainer

1. **A new pinned store (`linalg`)** or the small solves inside `lib/special`? (Section 6.)
2. **`--treated-value`**: is renaming the arm in the request the right interface, or should the verb require 0/1?
3. **A missing cell**: leave the unit out and count it (as `describe`), or refuse?
4. **`cycle detect` answering "none found" with exit 0** (section 3): agreed? A nonzero exit would be easier for a shell but is a false report of failure.
5. **The Python's test is kept** although the days of a series are not exchangeable under a trend (section 5). Keep, or add a documented alternative (shuffling blocks of days) behind a flag in S6.3?
6. **The name**: `cycle` (a cycle in a daily series) or `periodogram` (what an agent would search for)? `stats periodogram` is the name `docs/pure-cancho.md` of the client's repository proposed; the verb it describes is a *tested* single-period search, which a periodogram is not.
