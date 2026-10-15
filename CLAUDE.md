# CLAUDE.md -- cancho-stats

Design before code: each verb gets a section in `docs/design.md` before its implementation, and every claim in it is either measured (with machine, threads, load and
minimum of N) or marked as a choice, a read, or open.

- **Agent-first, like cancho-table and cancho-ml.** The CLI, errors and logs follow cancho `docs/agent-toolbox.md` through the `contract/` package of `cancho-tools`,
  pinned by commit and never copied. Errors are data with a stable `stats.<name>` rule tag, an exit code and, where safe, a repair. No payload is ever logged.
- **Exactness where a sum decides the result; determinism.** The same input, seed and arguments give the same bytes at any thread count, and on every platform the CI
  workflow covers. A feature that cannot keep that is refused or documented as the exception, with a test.
- **Dependencies are pinned, never copied:** cancho-table's `facc`, `flt`, `reader` stores, and the numeric stores published by cancho-ml (cancho-ml#46).
- **No file over 2,000 lines. Every refusal has a rule tag. No input reaches a panic.**
- **Gates are written before the code they judge and must be able to fail.** Mutants: all killed, or each survivor explained.
- **Measured claims only.** A claim found false is corrected in place, in the document that made it. Never state that an answer "matches" another library beyond the
  bound that was fixed in advance and measured here.
- **No GPU and no model runs without asking.**
- Licence EUPL-1.2, copyright Alpibru S.L. Commit trailer `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>`; PRs stay draft until the maintainer says merge; no
  force-push to main.
