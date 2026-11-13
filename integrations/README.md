# Agent integration

ORFscape can be driven by an LLM agent instead of a human. This folder holds one
shell-agnostic runbook and thin adapters that point an agent at it.

## What is where

| Path | Purpose |
|---|---|
| `AGENT_RUNBOOK.md` | **Single source of truth.** Role, order of operations, stage gates, conventions, stop rules, report format. |
| `scripts/preflight.py` | Environment + input check. Any FAIL means stop. |
| `scripts/stage_check.py` | Deep per-stage gates (catches "finished but empty"). |
| `scripts/run_pipeline.py` | Stage driver: runs one stage, checks it, stops on failure. Delegates to `orfscape.orchestrator` when the package is importable, so the run also emits the JSONL event stream. |
| `zcode/` | ZCode skill (`SKILL.md`) and subagent (`agents/orfscape-runner.md`). |
| `claude-code/` | `CLAUDE.md` project snippet and a Claude Code subagent definition. |
| `generic/PROMPT.md` | Copy-paste system prompt for any chat model (DeepSeek, GPT, Qwen, local models…). |

## Wiring it up

**ZCode** — copy `zcode/skills/orfscape-proteogenomics/` into `.zcode/skills/` and
`zcode/agents/orfscape-runner.md` into `.zcode/agents/` of the workspace. Then
`/orfscape-proteogenomics` loads the skill, or a task can be delegated to the
`orfscape-runner` subagent.

**Claude Code** — copy `claude-code/CLAUDE.md` to the project root (or append it to
your existing one) and `claude-code/.claude/agents/orfscape-runner.md` into
`.claude/agents/`.

**Anything else** — give the model `generic/PROMPT.md` as its system prompt and let
it read `AGENT_RUNBOOK.md` for the details. No tool-specific syntax is required
beyond "run a shell command and read its output".

All adapters assume the same three things: the repository root is the working
directory, `python -m orfscape` is available, and the agent may run shell commands.
Nothing else in this folder is executable.

## The contract in one paragraph

The agent runs `preflight.py` first and stops on any FAIL. It then drives stages in
order via `run_pipeline.py`, runs `stage_check.py` after every stage, and treats a
failed gate as a hard stop. It reports concrete numbers (tasks, PSMs, peptides,
KEEP counts) rather than narration, quotes WARN lines verbatim, and never re-runs a
completed stage "to be sure". Progress is visible to humans through
`00_logs/progress.jsonl`, which the GUI can attach to while the agent is still running.
