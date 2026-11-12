#!/usr/bin/env bash
# Reject staged paths and content that look like LLM conversation exports or
# Cursor/agent prompt caches. Does not flag legitimate docs that mention AI.

set -euo pipefail

MSG="Remove the artifact from the commit or use 'git commit --no-verify' to override."

if ! git rev-parse --git-dir >/dev/null 2>&1; then
  exit 0
fi

# Path patterns (case-insensitive via grep).
bad_path_re='(^|/)(\.cursor/|agent-transcripts/|agent-tools/|\.git/cursor/)|\.jsonl$|(^|/)AGENTS\.md$|(^|/)\.cursorrules$|(^|/)prompt-history\.(json|jsonl|txt)$'

staged_paths="$(git diff --cached --name-only --diff-filter=ACMRT 2>/dev/null || true)"
if [[ -z "${staged_paths}" ]]; then
  exit 0
fi

path_hits=()
while IFS= read -r p; do
  [[ -z "$p" ]] && continue
  if echo "$p" | grep -Eiq "$bad_path_re"; then
    path_hits+=("$p")
  fi
done <<< "$staged_paths"

if [[ ${#path_hits[@]} -gt 0 ]]; then
  echo "Pre-commit: refusing paths that look like LLM/agent history:"
  printf '  %s\n' "${path_hits[@]}"
  echo "$MSG"
  exit 1
fi

# Content patterns in staged text-like files only.
content_re='user_query>|system_reminder>|available_subagent_types|"role"[[:space:]]*:[[:space:]]*"assistant".*"tool_calls"'

content_hits=()
while IFS= read -r p; do
  [[ -z "$p" ]] && continue
  if [[ ! -f "$p" ]]; then
    continue
  fi
  case "$p" in
    tools/check-no-prompt-artifacts.sh|tools/verify-no-prompt-history.sh)
      continue
      ;;
    *.png|*.jpg|*.jpeg|*.gif|*.webp|*.ico|*.pdf|*.rda|*.tar|*.tar.gz|*.zip|*.bin)
      continue
      ;;
  esac
  if grep -q -E "$content_re" "$p" 2>/dev/null; then
    content_hits+=("$p")
  fi
done <<< "$staged_paths"

if [[ ${#content_hits[@]} -gt 0 ]]; then
  echo "Pre-commit: refusing staged files with LLM conversation metadata:"
  printf '  %s\n' "${content_hits[@]}"
  echo "$MSG"
  exit 1
fi

exit 0
