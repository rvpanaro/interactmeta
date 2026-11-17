#!/usr/bin/env bash
# Scan all reachable commits for identifiable agent prompt / chat metadata in text paths.

set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

textglob='*.R *.Rmd *.md *.yml *.yaml *.TeX *.tex *.bib *.json *.jsonl *.sh *.css *.html *.txt'
# pattern_expr kept on one line so this script is not flagged by its own pre-commit grep.
pattern_expr='user_query>|system_reminder>|available_subagent_types|agent-transcripts/'

commits=$(git rev-list --all | wc -l | tr -d ' ')
hits=0
for c in $(git rev-list --all); do
  if git grep -q -i -E "$pattern_expr" "$c" -- $textglob 2>/dev/null; then
    echo "FAIL commit $c"
    git grep -i -E "$pattern_expr" "$c" -- $textglob 2>/dev/null || true
    hits=$((hits + 1))
  fi
done

echo "commits_scanned=$commits content_hits=$hits"
if git log --all --format='%s' | grep -Eiq 'user_query|chatgpt|cursor agent'; then
  echo "FAIL: suspicious commit message(s)"
  hits=$((hits + 1))
fi

if [[ $hits -eq 0 ]]; then
  echo "OK: no identifiable prompt metadata in text history"
  exit 0
fi
exit 1
