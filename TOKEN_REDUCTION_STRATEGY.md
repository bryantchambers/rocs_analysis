# TOKEN REDUCTION STRATEGY

## Purpose

This repository is configured for aggressive token reduction when used with coding agents such as Codex, Claude Code, Gemini CLI, Cursor, or similar tools.

## Integration Status In This Repo (2026-05-27)

- `scripts/start_codex.sh`: supports `codex-squeezr` and `codex-full-context`.
- `scripts/start_codex_v1_bridge.sh`: supports `codex-v1-squeezr`,
  `codex-v1-full-context`, and compatibility aliases `codex-squeezr`,
  `codex-full-context`.
- `scripts/start_codex_v2.sh`: supports `codex-v2-squeezr`,
  `codex-v2-full-context`, and compatibility aliases `codex-squeezr`,
  `codex-full-context`.
- `scripts/v2/check_tools.sh` reports host availability for `rtk`, `sqz`,
  `squeezr`, and `squeezr-mcp`.

Mode selection:

- Host default at container launch: `CODEX_TOKEN_MODE=balanced|proxy|full-context`
- In-container switching (no restart): exit Codex and relaunch with
  `codex`, `codex-squeezr`, or `codex-full-context`.

Use the installed token-reduction stack whenever available:

1. **RTK**: reduce shell command output.
2. **sqz**: reduce repeated context and MCP/tool output.
3. **Squeezr AI**: reduce model/API traffic through local proxying.

Correctness always comes before token savings. If compression hides information needed for debugging, inspect the raw source or raw command output narrowly.

---

## General Rules

Minimize unnecessary context.

Do not dump large outputs into the conversation.

Prefer:

- targeted search
- targeted file reads
- summarized command output
- small diffs
- narrow tests
- incremental investigation

Avoid:

- full recursive directory listings
- full logs
- full lockfiles
- full generated files
- full test-suite output when one test is enough
- repeated reads of the same file
- pasting minified, vendored, cached, or generated artifacts

Before reading more files, summarize what is already known and why the next read is necessary.

---

## Tool Roles

### RTK

RTK is the primary shell-output reducer.

Use RTK for noisy shell commands, especially:

- git status and diffs
- logs
- test output
- build output
- directory discovery
- grep/search output

Preferred forms:

```bash
rtk git status
rtk git diff
rtk git log
rtk ls
rtk find
rtk grep
rtk test
rtk npm test
rtk pytest
rtk cargo test
```

Use raw commands only when RTK hides necessary details or when the user explicitly asks for raw output.

For raw debugging through RTK:

```bash
rtk proxy <command>
```

RTK meta commands must be run directly:

```bash
rtk gain
rtk gain --history
rtk discover
rtk proxy <command>
```

---

### sqz

sqz is the context and MCP compression layer.

Use sqz when available for:

- repeated context
- large command outputs
- multi-file summaries
- long logs
- repo summaries
- MCP/tool-output compression

Useful checks:

```bash
sqz gain
sqz stats
sqz --version
which sqz
which sqz-mcp
```

If sqz MCP tools are visible inside the agent, prefer them for summarizing or compressing large intermediate context.

If sqz is unavailable, continue with RTK and targeted native commands.

---

### Squeezr AI

Squeezr AI is the model/API proxy compression layer.

Use it for supported clients when running Codex, Gemini CLI, Claude Code, or similar tools.

Check status:

```bash
squeezr status
squeezr gain
squeezr config
```

For Codex CLI, prefer launching through Squeezr with:

```bash
HTTPS_PROXY=http://localhost:8081 codex
```

If a local alias exists, use:

```bash
codex-squeezr
```

Do not globally export `HTTPS_PROXY` unless the user explicitly wants all HTTPS traffic in the shell routed through Squeezr.

For Gemini CLI, use the configured Squeezr base URL if present:

```bash
echo "$GEMINI_API_BASE_URL"
gemini
```

If Squeezr is not running, either start it or continue without it:

```bash
squeezr start
squeezr status
```

---

## Startup Verification

At the beginning of a session, verify tools only when relevant to the task.

Use this compact check:

```bash
which rtk || true
rtk --version || true
rtk gain || true

which sqz || true
sqz --version || true
sqz gain || sqz stats || true

which sqz-mcp || true

which squeezr || true
squeezr status || true
squeezr gain || true
```

For Codex CLI through Squeezr:

```bash
HTTPS_PROXY=http://localhost:8081 codex
```

For Gemini CLI through Squeezr:

```bash
echo "$GEMINI_API_BASE_URL"
gemini
```

---

## Shell Command Policy

Always prefer targeted commands.

Good first commands:

```bash
pwd
git status --short
git diff --stat
git diff --name-only
rg --files
rg -n "pattern"
```

Better when RTK is available:

```bash
rtk git status
rtk git diff
rtk git log
rtk find
rtk grep "pattern"
```

Avoid starting with:

```bash
find .
ls -R
cat large.log
cat package-lock.json
cat pnpm-lock.yaml
cat yarn.lock
cat Cargo.lock
git diff
npm test
pytest
cargo test
```

Use narrow alternatives:

```bash
rg --files | head -n 200
git diff --stat
git diff --name-only
git diff -- path/to/file
tail -n 120 path/to/log
sed -n '1,220p' path/to/file
rg -n "error|failed|exception|panic|traceback" path/to/log
```

---

## File Reading Policy

Read files only after locating the relevant target.

Preferred workflow:

1. Identify candidate files.
2. Search for relevant symbols or errors.
3. Read the smallest useful range.
4. Summarize findings.
5. Continue only if needed.

Use:

```bash
rg --files
rg -n "function_name|class_name|error_message"
sed -n '1,220p' path/to/file
sed -n '220,440p' path/to/file
```

Avoid:

```bash
cat path/to/large/file
cat path/to/generated/file
cat path/to/lockfile
cat path/to/minified/file
```

Do not read vendored dependencies, generated output, caches, or binary artifacts unless directly relevant.

Common directories to avoid unless needed:

```bash
node_modules
dist
build
target
.next
.nuxt
.cache
coverage
vendor
.venv
__pycache__
.git
```

---

## Git Policy

Use compact Git inspection first.

Preferred:

```bash
rtk git status
git status --short
git diff --stat
git diff --name-only
```

Then inspect specific files:

```bash
git diff -- path/to/file
```

For history:

```bash
rtk git log
git log --oneline -n 20
git log --stat -n 5
```

Do not paste large diffs unless the user explicitly asks.

---

## Test Policy

Run the narrowest useful test first.

Preferred order:

1. Static search.
2. Type check or compile check.
3. Single failing test.
4. Package-level test.
5. Full test suite only if needed.

Examples:

```bash
rtk npm test -- path/to/test
rtk pytest path/to/test.py -q
rtk cargo test test_name
rtk cargo check
```

If RTK is unavailable:

```bash
npm test -- path/to/test
pytest path/to/test.py -q
cargo test test_name
cargo check
```

When reporting test failures, include only:

- command run
- short failure summary
- first relevant error
- likely cause
- next action

Do not paste full test output unless requested.

---

## Log Policy

For logs, inspect the end and search errors first.

Preferred:

```bash
tail -n 120 path/to/log
rg -n "error|failed|exception|panic|traceback|fatal|timeout" path/to/log
```

For context around a match:

```bash
sed -n 'LINE_START,LINE_ENDp' path/to/log
```

Replace `LINE_START` and `LINE_END` with a narrow range around the relevant match.

Avoid:

```bash
cat path/to/log
```

---

## Large Output Policy

If a command may produce large output:

1. Add filters.
2. Limit line counts.
3. Prefer summaries.
4. Use RTK.
5. Redirect to a file only if necessary, then inspect targeted sections.

Examples:

```bash
command | head -n 120
command | tail -n 120
command 2>&1 | tee /tmp/agent-output.log | tail -n 120
rg -n "error|failed|exception" /tmp/agent-output.log
```

---

## Codex Usage

When using Codex CLI, prefer Squeezr:

```bash
HTTPS_PROXY=http://localhost:8081 codex
```

If configured:

```bash
codex-squeezr
```

Inside Codex, check available MCP tools if needed:

```text
/mcp
```

If sqz or Squeezr tools are visible, use them for compression/status checks instead of manually dumping long context.

---

## Gemini CLI Usage

For Gemini CLI, mirror this file into `GEMINI.md` when possible.

Check whether Squeezr is configured:

```bash
echo "$GEMINI_API_BASE_URL"
```

Then run:

```bash
gemini
```

Gemini-specific instructions should go in:

```bash
GEMINI.md
```

or globally in:

```bash
~/.gemini/GEMINI.md
```

---

## Claude Code Usage

RTK may be configured through Claude Code hooks.

If hooks are active, normal shell commands may be rewritten automatically.

Even with hooks, prefer explicit RTK commands when output may be large:

```bash
rtk git status
rtk git diff
rtk grep "pattern"
rtk find
```

If raw output is required:

```bash
rtk proxy <command>
```

If project-specific Claude instructions exist, keep token policy in:

```bash
.claude/CLAUDE.md
.claude/.RTK.md
```

---

## Fallback Behavior

If a tool is missing, broken, or not on PATH:

1. Continue with native targeted commands.
2. Avoid broad output.
3. Mention the missing tool only if it changes the result.
4. Do not stop the task just because token-reduction tooling is unavailable.

Example fallback checks:

```bash
command -v rtk || echo "rtk not found"
command -v sqz || echo "sqz not found"
command -v squeezr || echo "squeezr not found"
```

---

## Safety and Correctness

Do not let token reduction hide correctness issues.

Use raw or broader inspection when needed for:

- security-sensitive code
- data loss risks
- migration scripts
- deployment scripts
- authentication
- permissions
- billing/payment code
- destructive commands
- failing tests with unclear cause

Before destructive commands, explain intent and prefer dry runs.

Examples of destructive commands requiring care:

```bash
rm -rf
git reset --hard
git clean -fdx
dropdb
kubectl delete
terraform destroy
```

---

## Final Response Policy

When reporting back:

- summarize what was checked
- cite commands run when relevant
- include only important output
- explain failures concisely
- propose the next concrete step
- avoid dumping raw logs or full diffs

Prefer concise, actionable answers.

---

## Complete Local Setup Snippets

Use these snippets to install or repair local shell integration.

### Add Cargo binaries to PATH

```bash
if ! grep -q 'HOME/.cargo/bin' "$HOME/.bashrc" 2>/dev/null; then
  printf '\nexport PATH="$HOME/.cargo/bin:$PATH"\n' >> "$HOME/.bashrc"
fi

if [ -f "$HOME/.zshrc" ] && ! grep -q 'HOME/.cargo/bin' "$HOME/.zshrc" 2>/dev/null; then
  printf '\nexport PATH="$HOME/.cargo/bin:$PATH"\n' >> "$HOME/.zshrc"
fi

export PATH="$HOME/.cargo/bin:$PATH"
```

### Add common global npm binary locations to PATH

```bash
if [ -d "$HOME/.npm-global/bin" ]; then
  export PATH="$HOME/.npm-global/bin:$PATH"
fi

if [ -d "/opt/homebrew/bin" ]; then
  export PATH="/opt/homebrew/bin:$PATH"
fi

if [ -d "$HOME/.pixi/envs/nodejs/bin" ]; then
  export PATH="$HOME/.pixi/envs/nodejs/bin:$PATH"
fi
```

### Create a Codex alias for Squeezr

```bash
if ! grep -q 'alias codex-squeezr=' "$HOME/.bashrc" 2>/dev/null; then
  printf '\nalias codex-squeezr='\''HTTPS_PROXY=http://localhost:8081 codex'\''\n' >> "$HOME/.bashrc"
fi

if [ -f "$HOME/.zshrc" ] && ! grep -q 'alias codex-squeezr=' "$HOME/.zshrc" 2>/dev/null; then
  printf '\nalias codex-squeezr='\''HTTPS_PROXY=http://localhost:8081 codex'\''\n' >> "$HOME/.zshrc"
fi
```

### Verify token-reduction tools

```bash
echo "== RTK =="
command -v rtk || true
rtk --version || true
rtk gain || true

echo "== sqz =="
command -v sqz || true
sqz --version || true
sqz gain || sqz stats || true

echo "== sqz-mcp =="
command -v sqz-mcp || true

echo "== Squeezr =="
command -v squeezr || true
squeezr status || true
squeezr gain || true

echo "== Codex =="
command -v codex || true

echo "== Gemini =="
command -v gemini || true
echo "GEMINI_API_BASE_URL=${GEMINI_API_BASE_URL:-}"
```

### Start Squeezr and launch Codex

```bash
squeezr start || true
squeezr status || true
HTTPS_PROXY=http://localhost:8081 codex
```

### Launch Gemini

```bash
echo "GEMINI_API_BASE_URL=${GEMINI_API_BASE_URL:-}"
gemini
```


<!-- BEGIN sqz-agents-guidance (auto-installed by sqz init; remove this block to disable) -->

## sqz — Token-Optimized CLI Output

When running shell commands whose output may be long (directory listings,
git log/diff, test runners, build logs, `docker ps`, `kubectl get`, etc.),
pipe the output through `sqz compress` to reduce token consumption.

`sqz` is a stdin-to-stdout compressor, not a command wrapper. The correct
usage is to pipe the command's output into `sqz compress`:

```bash
# Instead of:     Use:
git status        git status 2>&1 | /opt/homebrew/bin/sqz compress
cargo test        cargo test 2>&1 | /opt/homebrew/bin/sqz compress
git log -10       git log -10 2>&1 | /opt/homebrew/bin/sqz compress
docker ps         docker ps 2>&1 | /opt/homebrew/bin/sqz compress
ls -la            ls -la 2>&1 | /opt/homebrew/bin/sqz compress
```

The `2>&1` captures stderr too, which is useful for commands like `cargo
test` where diagnostics go to stderr. `sqz compress` filters and compresses
the combined output while preserving filenames, paths, and identifiers.
It typically saves 60-90% tokens on verbose commands.

Do NOT pipe output for:
- Interactive commands (`vim`, `ssh`, `python`, REPLs)
- Compound commands with shell operators (`cmd && other`, `cmd > file.txt`,
  `cmd; other`) — run those directly
- Short commands whose output is already a few lines

If `sqz` is not on PATH, run commands normally.

The `sqz-mcp` MCP server is also available — Codex reads it from
`~/.codex/config.toml` under `[mcp_servers.sqz]`. It exposes three
tools: `compress` (the default pipeline), `passthrough` (return text
unchanged — the escape hatch below), and `expand` (resolve a
`§ref:HASH§` token back to the original bytes).

## Escape hatch — when sqz output confuses you

If you see a `§ref:HASH§` token and can't parse it, or compressed
output is leading you to make lots of small retries instead of one
big request, use one of these:

- **`/opt/homebrew/bin/sqz expand <prefix>`** — resolve a dedup ref back to the
  original bytes. Accepts bare hex (`sqz expand a1b2c3d4`) or the full
  token pasted verbatim (`sqz expand §ref:a1b2c3d4§`).
- **`SQZ_NO_DEDUP=1`** — set this env var for one command to disable
  dedup: `SQZ_NO_DEDUP=1 git status 2>&1 | sqz compress`. You'll get
  the full compressed output with no `§ref:…§` tokens.
- **`--no-cache`** — same opt-out as a CLI flag:
  `git status 2>&1 | sqz compress --no-cache`.

If you're using the MCP server, the `passthrough` tool returns raw
text and the `expand` tool resolves refs — call them when you need
data sqz hasn't touched.

<!-- END sqz-agents-guidance -->
