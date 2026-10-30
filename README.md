# cancho-stats

Statistics on tabular data as a command line an agent can read: sums that are exact, answers that are the same bytes everywhere, and a refusal that names its rule instead of a made-up number.

> **Kind:** cli · **Tier:** 1 · **Stage:** built

## One complete flow

`stats describe` on the eight readings of `docs/example/readings.csv`, whose `drift` column is a large offset plus a small change (`100000000.1` to `100000000.5`), and then a test that is refused because a Monte Carlo draw needs a seed. Both outputs below are printed by `build/stats` when this README is generated (the second is an excerpt of the one JSON line it prints); nothing is typed by hand.

<!--run:flow-->
```console
$ stats describe --cols drift --format text docs/example/readings.csv
rows = 8
ddof = 1
quantile_type = 7
rounding = "mean, variance, skewness, excess kurtosis, quantiles: exact, rounded once; sd: the square root of the rounded variance, two roundings"
columns[0].name = "drift"
columns[0].kind = "numeric"
columns[0].count = 8
columns[0].missing = 0
columns[0].min = 100000000.1
columns[0].max = 100000000.5
columns[0].mean = 100000000.275
columns[0].variance = 0.022142858015639467
...

$ stats resample perm --cols bake_min --group oven docs/stories/bakery.csv
{
  "code": "INVALID_ARGS",
  "rule": "stats.seed-required",
  "hint": "stats resample perm --cols A --group G --seed 1 FILE, or --permutations all"
}
$ echo $?
2
```
<!--/run-->

The mean `100000000.275` is the exact mean of the eight numbers, rounded once. A refusal has a stable `rule`, a `hint` that says how to fix the request, and an exit code (here 2, a malformed invocation).

## What is special

| | what | evidence | what it does not cover |
|---|---|---|---|
| **Exact** | sums, squares and products are exact and rounded once; the order of the rows cannot change a bit | gates against `Fraction` and 60-digit references ([numerics](docs/numerics.md)) | the model, the assumptions, the sample, the question: see the [stories](docs/stories.html) |
| **Same bytes** | the same answer on Linux x86-64, Linux aarch64 and macOS arm64, at 1 to 16 threads | a committed fingerprint, [`docs/xplat_fp.txt`](docs/xplat_fp.txt), checked on every push | only what is built is fingerprinted |
| **Refusals are data** | a rule tag, an exit code and a repair, never a trap or an invented number | `stats.v1`, one JSON line; a test per rule | |
| **Bounds first** | every function has an error bound committed before its code | the bound, a reference that tries to break it, and the correction when one broke | the functions of S0; the rest have their own sections |
| **Authority is a fact** | no network, no clock, no write by path | [`manifests/stats.authority.json`](manifests/stats.authority.json), checked in CI | |

![Box plots of the baking time of two ovens: B is faster, with a long tail of slow batches.](docs/figures/ovens-boxes.svg)

## For humans, for agents

| for humans | for agents |
|---|---|
| the [project page](https://alpibrusl.github.io/cancho-stats/), the [stories](docs/stories.html) (by topic), [docs/design.md](docs/design.md) | `stats skill` (a guide, from the binary), `stats introspect` (flags, limits, rules and authority as one JSON document), [AGENTS.md](AGENTS.md), [llms.txt](llms.txt) |

## Quick start

The compiler is pinned in `cancho.toml` (by commit); `cancho install` fetches the pinned dependencies, and nothing is copied into this repository.

```console
$ cancho install && cancho build --bin stats
$ build/stats describe docs/example/readings.csv
```

## Examples

Welch's t test of the baking time of two ovens (`docs/stories/bakery.csv`, invented, 80 batches), the interval for the difference, and the exact-or-Monte-Carlo permutation test that needs no assumption about the shape of the data. Each output is excerpted from a real run.

<!--run:t2-->
```console
$ stats test t2 --cols bake_min --group oven --format text docs/stories/bakery.csv
test = "t2"
group = "oven"
levels.a = "A"
levels.b = "B"
rows = 80
level = 0.95
variance = "welch"
adjustment = "none"
tests[0].column = "bake_min"
tests[0].statistic = 1.8222469060619957
tests[0].df = 73.1757032479013
tests[0].p_value = 0.0725033707662692
...
```
<!--/run-->

<!--run:ci-->
```console
$ stats ci diff --cols bake_min --group oven --format text docs/stories/bakery.csv
ci = "diff"
group = "oven"
levels.a = "A"
levels.b = "B"
rows = 80
level = 0.95
variance = "welch"
intervals[0].column = "bake_min"
intervals[0].estimate = 1.9775000000000005
intervals[0].interval.of = "difference of means"
...
```
<!--/run-->

<!--run:perm-->
```console
$ stats resample perm --cols bake_min --group oven --seed 1 --permutations 9999 --format text docs/stories/bakery.csv
resample = "perm"
statistic = "mean-diff"
group = "oven"
levels.a = "A"
levels.b = "B"
rows = 80
p_method = "monte-carlo"
seed = 1
assumptions.shuffle = "Fisher-Yates of the pooled values (group a then group b) from the end: i = m - 1 - d, z = mix(s_b + (d + 1) gold), j = (z >> 1) mod (i + 1), swap x[i] and x[j]; group a* is x[0 .. n_a)"
assumptions.stream = "s_b = derive(derive(seed, 2), b), b from 0"
tests[0].column = "bake_min"
tests[0].statistic = 1.9775000000000005
...
```
<!--/run-->

Longer cases, with data, in the [stories](docs/stories.html), by topic. Two of them, as a picture:

![Conversion lift of an ad: the dashboard's 5.3 points against the randomized holdout's 1.4, with quiet and busy users in between; only the holdout excludes the selection.](docs/figures/dashboard-vs-holdout.svg)

![Three endpoints of a drug trial, raw and Holm-adjusted p-values: the diastolic result is below 0.05 alone and above it after the correction.](docs/figures/holm-pvalues.svg)

A case can span packages: the marketing case uses `cancho-table` for the ledger and `cancho-ml` for the uplift.

## Use it with an AI agent

`stats` is a command line, so an agent uses it through any tool that runs a command. It describes itself (`stats skill`, `stats introspect`), answers and refuses in one JSON schema, and a refusal says what to change. A real session, scripted:

<!--run:agent-->
```console
# 1. the agent asks the tool what it is, not this file
$ stats skill
---
name: stats
description: Describe, correlate and test the columns of a CSV file: statistics from exact sums rounded once (describe), Pearson, Spearman and Kendall coefficients (corr), the t family, one-way analysis of variance, Kruskal-Wallis, Mann-Whitney, Wilcoxon, Kolmogorov-Smirnov, chi-square and proportion tests with intervals, effect sizes and corrections for several tests (test); the same bytes on every platform and at any thread count.
---
...

# 2. it runs the test it wants, without a seed: a refusal, as data
$ stats resample perm --cols bake_min --group oven docs/stories/bakery.csv
{
  "code": "INVALID_ARGS",
  "rule": "stats.seed-required",
  "hint": "stats resample perm --cols A --group G --seed 1 FILE, or --permutations all"
}
$ echo $?
2

# 3. it applies the hint (a seed) and runs again
$ stats resample perm --cols bake_min --group oven --seed 1 --permutations 9999 --format text docs/stories/bakery.csv
resample = "perm"
statistic = "mean-diff"
group = "oven"
levels.a = "A"
levels.b = "B"
rows = 80
p_method = "monte-carlo"
seed = 1
assumptions.shuffle = "Fisher-Yates of the pooled values (group a then group b) from the end: i = m - 1 - d, z = mix(s_b + (d + 1) gold), j = (z >> 1) mod (i + 1), swap x[i] and x[j]; group a* is x[0 .. n_a)"
...
```
<!--/run-->

What the agent cannot do is part of the design: the authority row has no network, no clock and no write by path, so a call cannot send the data anywhere or leave a file behind. There is no MCP server for `stats` yet (`cancho-table` and `cancho-ml` have one; it is a gap, not a decision).

## Compared with

| | cancho-stats | the usual alternatives (SciPy, statsmodels, R) | where cancho-stats is worse |
|---|---|---|---|
| sums | exact, rounded once; the order of the rows cannot change a bit | floating-point sums; whether a routine is order-independent was not measured here | the exactness has a cost in time that has not been measured |
| output | one JSON schema, named refusals, repairs | library calls and printed text | no plotting, no model objects to keep between calls |
| reproducibility | the same bytes on three CI platforms, at 1 to 16 threads | not measured here | only what is built is covered by the fingerprint |
| accuracy of the special functions | each has a stated bound and a 60-digit reference; `scipy.special` served as a second reference and, on the binomial, is the outlier by up to 3.4e-9 ([docs/numerics.md](docs/numerics.md)) | as measured there | that comparison covers the functions of S0 only |
| breadth | the verbs listed under *What exists* | far more methods (not counted here) | **much narrower; no regression yet (planned)** |

Nothing here is a claim about speed or about agreement with another library beyond the bound that was fixed in advance and measured in `docs/numerics.md`.

## What exists, what cannot yet

Built and measured: `describe`, `corr`, `test` (`t2`, `paired`, `mannwhitney`, `wilcoxon`, `anova`, `kruskal`, `ks`, `fisher`, `chisq`, `prop`) with Holm and other corrections, `ci` (`prop`, `mean`, `diff`, `corr`, `boot` with percentile and BCa), and `resample perm` for the difference of means, Monte Carlo or exact. Not built: other statistics for `resample perm`, paired sign flips and the correlation permutation (S5.4b, S5.4c), and regression (`lm`, `glm`). The design of each verb, with its status and measured results, is in [docs/](docs/); [docs/design.md](docs/design.md) has the milestones.

## How fast

**No speed is claimed.** The command line has not been timed (`docs/describe.md` lists speed as not measured). The measured claims are about accuracy and reproducibility, above.

## Learn more

[docs/design.md](docs/design.md) (the whole tool), [docs/numerics.md](docs/numerics.md) (the bounds and the gates), [docs/resample.md](docs/resample.md), [docs/stories.html](docs/stories.html), the [project page](https://alpibrusl.github.io/cancho-stats/). The README, the project page and the stories follow the package kit of [cancho-kit](https://github.com/alpibrusl/cancho-kit).

## Contributing

Read [CLAUDE.md](CLAUDE.md): design before code (each verb gets a section of `docs/design.md` before its implementation), gates written before the code they judge and able to fail, no file over 2,000 lines, every refusal has a rule tag.

## Licence

EUPL-1.2, copyright Alpibru S.L. See [LICENSE](LICENSE) and [NOTICE](NOTICE).
