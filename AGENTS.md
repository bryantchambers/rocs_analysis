# Agent Vocabulary

Use these terms consistently in this repository.

- `scan`: a broad rough-cut grid search across wide parameter ranges to find a usable region. This is the first pass.
- `broad scan` or `broad grid`: the same thing as `scan`, used when the width of the search space matters.
- `sweep`: a systematic parameter search. In this repo, `sweep` is often used as a synonym for `scan`, but it should be qualified when possible.
- `neighborhood scan`: a local refinement search around one or more promising settings from a prior broad scan.
- `grid search`: the concrete method used to execute a `scan` or `sweep`, usually over a predefined set of parameter combinations.
- `refinement run`: a narrower follow-up search after a broad scan, usually centered on the best local hit.
- `full scan`: a broad scan that covers the intended global parameter space for the current question.

Rules of use:

- Use `scan` for the broad first pass.
- Use `neighborhood scan` for local follow-up runs around a known good setting.
- Avoid using `sweep` alone when the scope is ambiguous.
- Prefer `grid` when referring to the actual set of parameter combinations, not the scientific purpose of the run.

# Token Reduction

Use the installed token-reduction tools to keep agent runs compact.

- Use `rtk` for noisy shell output such as `git status`, diffs, logs, build
  output, and test output.
- Use `sqz` for repeated context, multi-file summaries, long logs, and
  MCP/tool-output compression.
- Do not use the legacy proxy tooling in this repository.
- Prefer narrow reads and summaries before asking for more file content.
