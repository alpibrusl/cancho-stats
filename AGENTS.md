# cancho-stats, for an agent

`stats` describes, correlates and tests the columns of a CSV file and draws intervals from them. Ask the tool, not this file:

```console
$ stats skill        # a guide, in Markdown, from the binary
$ stats introspect   # one JSON document: the flags, the limits, the refusal rules and the authority
```

* **Answers and refusals are one JSON line**, schema `stats.v1`: `{ok, command, schema, data, error, errors, meta}`. Match on `error.rule` (`stats.<name>`), never on `message`. The `hint` says how to fix the request; `repair`, when present, is an argv to run.
* **Exit codes:** 0 done; 2 a malformed invocation or an unsupported combination; 8 the data refused (for example a short record, `stats.csv-ragged`). `stats introspect` lists them all.
* **A Monte Carlo draw needs `--seed`** (`stats.seed-required`): an unseeded draw is not reproducible, and the tool has no clock. The same data, seed and arguments give the same bytes at any thread count.
* **It reads one file under `--root` and prints to standard output.** It has no network, no clock and cannot write a file; if a task needs one of those, it is not a task for this tool.
* **Say what the number does not show.** An exact p-value is not evidence that the model fits, the sample was random, or the question was the only one asked; correct many tests (`--correct holm`) and say how many were run.
* The stories in `docs/stories.html` are worked cases with data; the design of each verb is in `docs/`.
