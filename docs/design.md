# cancho-stats: design draft

**Status: accepted by the maintainer on 2026-10-15 (new repository, formula grammar in v1, the whole test table in v1, the numeric stores of `cancho-ml` pinned, `exp` and `ln` from cancho). This is the design as accepted; what has been built since (S0, S0b and most of S1) is in `docs/numerics.md`, `docs/describe.md` and the issues, and a paragraph here that the code has overtaken says so.** Tags: **[read]** read from a
repository or document, **[measured]** printed by a script, **[choice]** a decision the maintainer can overturn, **[open]** not decided. A repository is created only
when the maintainer says so. Written 2026-10-15.

## 1. What it is, and what it is not

`stats` is a command line, in the shape of `table` and `ml`, that **describes data and draws inferences from it**: summaries, correlations, hypothesis tests,
confidence intervals, resampling, distribution fits, and regression *with* standard errors and p-values. Its answers are data: every number is printed as the
decimal text of a correctly rounded double, every refusal has a rule tag and a repair, and the same input and seed give the same bytes at any thread count.

It is not a model-fitting tool (that is `ml`), not a query tool (that is `table`), not a plotting tool, and not a server. It stores nothing: no model file, no state.

## 2. Where it sits

| need | owner | evidence |
|---|---|---|
| read a CSV, exact types, group-by, exact `count`, `sum`, `min`, `max`, `mean`, `distinct` | `cancho-table` | **[read]** its README: sums are exact at any width and rounded once |
| fit, store, apply and serve a predictive model; cross-validation; model metrics | `cancho-ml` | **[read]** six models, `ml-predict`, `ml-serve`, issues #42 (splits and CV), #43 (tracking) |
| test, estimate, interpret: p-values, intervals, effect sizes, regression inference, resampling | **`cancho-stats`** | this document |

**[choice]** A rule to keep the line sharp. Anything that produces or scores a *stored model* lives in `ml` (confusion matrices, calibration, a bootstrap interval on a
model's metric, cross-validation). Anything that answers "is this difference or this relationship real, and how large" lives in `stats`. Linear regression *fitting* stays in
`ml`; the *inference* on a regression (standard errors, t statistics, p-values, confidence intervals, leverage, Cook's distance) is `stats lm`, which uses `ml`'s linear
algebra rather than duplicating it.

**[choice] A new repository, not a fourth binary of `cancho-ml`.** The three binaries of `cancho-ml` differ by authority row around one job (a model file); `stats` has a
different job, a different read-only authority, its own verbs and rule tags, and its own issue tracker and release rhythm. It pins `cancho-table` and `cancho-ml` as
dependencies the way `cancho-ml` already pins `cancho-table`.

## 3. Principles it inherits

* **Exactness where a sum decides the result.** A mean, a variance and a covariance are sums; they are accumulated exactly (`facc`, **[read]** in `cancho-table`) and rounded
  once. The usual two-pass or Welford update is *not* used where it would make the last bits depend on the row order.
* **Determinism.** Same input, seed and arguments give the same bytes at any thread count. Resampling draws from one stream per replicate, derived from the seed and the
  replicate's index (the derivation `lib/core/rng.cho` already has, **[read]**), so adding threads can never change a replicate. Cross-machine: `ml`'s fingerprint
  gate shows macOS arm64 and Linux x86-64 agree on 88 lines for 22 fits **[measured, 2026-10-15, cb08d5b]**; `stats` repeats the method with its own fingerprint.
* **Refusals are data.** A rule tag `stats.<name>`, an exit code, a repair where one is safe, the numbers that decided it. No input reaches a panic.
* **Bounded resources.** Rows, columns, replicates, categories and table cells each have a flag, a default and a ceiling, checked before anything sized by the input.
* **Agent-first.** `introspect`, `skill`, one flag table that drives the parser, the schemas and any MCP surface; the contract package of `cancho-tools`, pinned.
* **Gates are written before the code and must be able to fail.** Mutants all killed or each survivor explained.
* **No claim about another library that was not measured here.** "Matches scipy to 1e-12 relative on these 40,000 cases" is a claim; "as accurate as scipy" is not.

## 4. Verbs (the proposed v1) **[choice]**

One binary, `stats`, with these verbs. Each is a milestone with its own gates (section 8).

| verb | does | notes |
|---|---|---|
| `stats describe FILE [--by COL]` | per column: count, missing, distinct, min, max, mean, variance, standard deviation, skewness, kurtosis, quantiles, an exact histogram | quantile definition is a flag (`--quantile-type`, the nine Hyndman-Fan types, default stated) |
| `stats corr FILE --cols A,B,...` | Pearson, Spearman (average ranks), Kendall tau-b; the matrix and the p-values | the tie rule is data in the answer |
| `stats test KIND ...` | `t1`, `t2` (Welch default), `paired`, `anova`, `kruskal`, `mannwhitney`, `wilcoxon`, `ks`, `chisq`, `fisher`, `prop` | every answer: statistic, degrees of freedom, p-value, effect size, an interval, n per group; multiple-comparison correction (`bonferroni`, `holm`, `bh`) over a declared family |
| `stats ci KIND ...` | intervals: mean (t), proportion (Wilson, Clopper-Pearson), difference, correlation (Fisher z), and a seeded percentile or BCa bootstrap of any `describe` statistic | `--replicates R --seed S`, bounded |
| `stats resample KIND ...` | permutation tests and bootstraps for a declared statistic | the statistic is chosen from a closed list, not an expression |
| `stats lm FILE --formula ...` | OLS with standard errors, t, p, intervals, R squared, adjusted R squared, F, residual summary, leverage, Cook's distance, VIF; HC3 robust errors on request | a binary-outcome sibling `stats glm` (logistic, Wald intervals) later |
| `stats fit-dist FILE --col C --family F` | maximum-likelihood fits (normal, lognormal, exponential, gamma, Poisson, ...), with a goodness-of-fit test | families are a closed list |
| `stats power` | sample size and power for the tests above | closed form where it exists, else a seeded simulation |
| `stats compare` | tests on the output of `ml` cross-validation (paired t on folds with the corrected resampled variance, McNemar) | reads `ml`'s JSON, never a model file |

Time series (autocorrelation, decomposition, stationarity tests) and multivariate methods (PCA lives in `ml` as a pipeline step, #41) are deferred.

**[choice, decided 2026-10-15] A closed formula grammar in v1**, with its own parser, fuzz and rule tags:

```
formula  = response "~" terms
response = column | "log(" column ")"
terms    = term { ("+" | "-") term }          -- "- 1" or "+ 0" removes the intercept; the intercept is the default
term     = factor { ":" factor }              -- ":" is an interaction
         | factor "*" factor                  -- a*b means a + b + a:b
factor   = column | "C(" column ")"           -- C() marks a categorical column (treatment coding, the reference level is the
                                              -- bytewise-smallest level unless "C(col, ref=LEVEL)")
         | "log(" column ")" | "sqrt(" column ")" | column "^" 2
column   = name of a header column, or a quoted name "..."
```

Nothing else: no arbitrary expressions, no function calls beyond `log`, `sqrt` and the square, each with a domain refusal (`stats.formula-domain`: a value that makes the transform undefined
is refused with the row, never turned into a missing cell). The expanded design matrix has a bounded width (a limit refused before allocation), the column names are printed
(`C(g)[b]`, `x:z`), and a rank-deficient matrix is `stats.collinear` naming the columns (the repair is the dropped term). Categorical columns are the string table the reader already builds.

## 5. Numerics

* **Exact sums, one rounding.** Count, sum, mean, variance and covariance come from exact accumulation. The variance is `(sum x^2 - (sum x)^2 / n) / (n - ddof)` evaluated exactly in
  an integer or rational form and rounded once; the formula and the `ddof` are stated per verb. **[open]** the width needed for `x^2` of large doubles: `facc` handles sums, and a
  product term needs a wider accumulator or an exact two-product; the design question is settled in S0 before `describe`.
* **Special functions.** The p-values need the normal, Student t, chi-square, F, beta and gamma distributions, which come from the regularised incomplete gamma and beta
  functions, the error function and `lgamma`. **[choice]** implemented here in cancho with a stated target (relative error under 1e-12 in the central range, a documented
  loss in the far tails), gated against `mpmath` at high precision (a script, not a dependency), with the tolerance written before the code. **[read, resolved 2026-10-15]** `exp`, `log`, `log1p`, `expm1`, `pow` and the trigonometric and
  hyperbolic functions are in `std/math.cho` of the compiler, written in cancho (fdlibm's algorithms and constants, bit for bit), not calls into a platform `libm`; their accuracy is
  measured against the C library in the compiler's `docs/float-math.md` section 9 (`log`: within 1 ulp). So they are the same code on every platform, and `ml`'s fingerprint agrees on three
  platforms **[measured]**. What `std.math` lacks and S0 must build: `lgamma`, `erf`/`erfc`, the regularised incomplete gamma and beta functions, the normal quantile; if they prove
  generally useful they can move to `std` later.
* **Ranks and ties.** One tie rule per method, stated in the answer (`average` for Spearman and Mann-Whitney, with the standard tie correction of the variance).
* **Quantiles, order statistics.** Exact selection, no floating-point interpolation until the last step.
* **No silent approximations.** An exact test (Fisher, small-sample Mann-Whitney) refuses above a stated size and names the asymptotic alternative as the repair.

## 6. Output contract

Every verb prints one JSON document: `{"ok":true,"command":"stats","schema":"stats.v1","data":{...}}` or the refusal, as `table` and `ml` do. A test answer carries
`{"test","statistic","df","p_value","effect_size":{"name","value"},"interval":{"level","low","high"},"n":{...},"method","assumptions":{...}}`. `assumptions` holds
*facts* (sample sizes, variance ratio, number of ties, expected counts below 5), never a verdict. Numbers are decimal text; `--format text` prints a table.

Rule tags **[choice]** (illustrative): `stats.sample-too-small`, `stats.variance-zero`, `stats.table-sparse`, `stats.exact-too-large`, `stats.family-too-large`,
`stats.formula`, `stats.collinear`, `stats.replicates-limit`; `describe` as built uses none of the first two (an undefined statistic is `null`, see `docs/describe.md` section 4), plus the reader's, the limits' and the arguments' tags shared with `table`.

## 7. Authority and surface

`stats`: args, heap, file read, standard output and error; threads (`conc`) only if resampling is threaded; no file write, no network, no listener, no exec, no clock. **[choice]** No
model file, no server. An MCP surface is derived from the flag table later, as `ml`'s is planned (#45). Output to a file is the shell's job.

## 8. Gates (written before the code, each able to fail)

1. **Independent oracles in Python**, using `fractions` for exact sums and `mpmath` for special functions, comparing the *bytes of the printed decimal* where the arithmetic is exact (counts, sums,
   means, medians, ranks, contingency tables) and a tolerance fixed in advance where it is not (p-values, intervals).
2. **scipy/statsmodels comparisons** on generated data, with the bound fixed before the first run, controls (a deliberately wrong formula must fall outside), and the statement of what
   was compared. No claim of equality beyond the bound.
3. **Calibration gate.** Under a true null with seeded simulation, the rejection rate at level 0.05 of each test lies inside the binomial bounds for N replicates, and the p-values
   are uniform (a KS test on them), for sample sizes from small to large. A wrong tail or a wrong degrees-of-freedom formula fails this even when it passes the point comparisons.
4. **Power/known-effect gate.** A planted effect of known size is detected at the rate the theory gives, within bounds.
5. **Thread and row-order invariance**; **cross-machine fingerprint** on the runners of the CI workflow (the `ml` method, **[measured]** for `ml`).
6. **Fuzz** of arguments and of hostile inputs (huge, tiny, constant, NaN-like cells, all-missing, one row): no trap.
7. **Mutants** for every verb: killed, or each survivor explained.

## 9. Repository, dependencies, licence

* **[choice]** `alpibrusl/cancho-stats`, public, EUPL-1.2, copyright Alpibru S.L. (the maintainer's standing rule for new repositories; `cancho-ml` carries the same text **[read]**).
* **Pins, never copies:** `cancho-table`'s `facc`, `flt`, `reader` stores (as `cancho-ml` pins them, **[read]**); `cancho-tools`' contract package; and from `cancho-ml` the
  modules `lib/core` (RNG, matrix, exact helpers), `lib/linalg` (Cholesky, solves, lane-ordered sums) and `lib/metrics` (ranks, sums). **The last needed the prerequisite issue
  (cancho-ml#46), which is done:** `tools.toml` pins its 23 numeric modules as stores (`docs/numerics.md` section 2, checked by `scripts/check_pins.sh`).
* **Compiler:** the revision the siblings pin, moved only deliberately.
* **CI:** the workflow `cancho-ml` is getting (a branch in progress), with the same fingerprint method.

## 10. Milestones

| | milestone | gate in one line |
|---|---|---|
| S0 | accuracy harness: exact accumulation of products, special functions, the determinism of `exp`/`ln`, the stores consumed | mpmath agreement bounds written first; the platform-math question answered |
| S1 | `describe` | oracle bytes for every exact statistic; quantile types vs a Python implementation |
| S2 | `corr` (design: `docs/corr.md`; done 2026-10-23) | oracle bytes (Pearson from exact sums, ranks); scipy bounds |
| S3 | `test` t family, `chisq`, `prop`, corrections (design and results: `docs/test.md`; done 2026-10-24) | calibration + oracle |
| S4 | `test` rank tests, `anova`, `kruskal`, `ks`, `fisher` | calibration + exact-test oracles |
| S5 | `ci`, `resample` | thread-invariant seeded streams; coverage gate (an interval covers the truth at its nominal rate) |
| S6 | `lm`, then `glm` | oracle for the normal equations in exact arithmetic; statsmodels bounds; leverage identities |
| S7 | `fit-dist`, `power`, `compare` | likelihood oracles; known-power gate |

## 11. GPU **[choice]**

**Not in v1, and the seam is kept honest rather than promised.** `cancho-gpu` **[read]** is a compiler from tile kernels to CUDA C and Metal with a CPU reference interpreter; its GPU-driving layer
is "built, checked against a mock driver; not yet run on a GPU". The reasons this tool does not start there:

* The pitch is exact sums and the same bytes everywhere. A GPU does not do exact accumulation, and Metal has no double precision, so a GPU path is a different numerical tier, not an
  acceleration of the same one. `cancho-ml#47` is the spike that measures whether such a tier is worth building for `ml`; it needs the maintainer's go for every GPU run.
* The data a statistics tool reads is mostly small (thousands to millions of rows); the heavy cases are resampling with many replicates, each cheap on the CPU. The CPU path comes first and is measured.
* If `ml#47` finds a tier worth having, `stats` would adopt it only for the pieces that are *not* exact (a bootstrap of a statistic, with its own tier tag in the answer), never for `describe`,
  `corr` or `test`, whose answers are exact or correctly rounded.

## 12. Decisions of 2026-10-15 and what is still open

Decided by the maintainer: (1) a new repository, `alpibrusl/cancho-stats` (created private, the visibility is the maintainer's to change); (2) the closed formula grammar of section 4 is in v1;
(3) every test in the table of section 4 is in v1 (milestones S3 and S4); (4) the numeric modules of `cancho-ml` are consumed as pinned stores (cancho-ml#46 first), nothing copied;
(5) `exp` and `ln` come from cancho: they are already in `std.math`.

Still open: the order inside S3 and S4; whether the special functions move to `std.math` later; the width needed for exact `x^2` terms (S0); the list of the `stats compare` inputs from `ml cv` (after cancho-ml#42).
