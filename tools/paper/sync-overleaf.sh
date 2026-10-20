#!/usr/bin/env bash
# Pushes paper/, as committed, to the Overleaf project behind the `overleaf` remote (Overleaf's Git
# integration), so that the project holds exactly the paper's files.
#
#     git remote add overleaf https://git@git.overleaf.com/6ad6c4d47e5cd74f01b8dde4    # once
#     make paper-sync              # after committing a change to paper/
#     make paper-sync FORCE=1      # the first time, to replace the new project's placeholder files
#
# Overleaf refuses forced pushes, so each sync is one new commit on top of the project's own history whose
# files are paper/'s. The sync is one way: paper/ is generated from tools/paper/, so an edit made in Overleaf
# would be lost at the next `make paper`. The script remembers the commit it last pushed (the local ref
# refs/overleaf/pushed) and stops when the project has moved past it, unless given --force; the edit then
# stays in the project's history, behind the new commit. (As in the dinara-align paper.)
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"

force=false
[ "${1:-}" = "--force" ] && force=true

if ! git remote get-url overleaf >/dev/null 2>&1; then
    echo "No overleaf remote: git remote add overleaf https://git@git.overleaf.com/<project id> (Overleaf: Menu > Git)" >&2
    exit 1
fi
if ! git diff --quiet HEAD -- paper || [ -n "$(git ls-files --others --exclude-standard paper)" ]; then
    echo "paper/ has changes that are not committed; commit them, then sync" >&2
    exit 1
fi

# The project's one branch, whatever Overleaf calls it (main for new projects, master for old ones).
branch=$(git ls-remote --symref overleaf HEAD | sed -n 's|^ref: refs/heads/\(.*\)\tHEAD$|\1|p')
branch=${branch:-main}
git fetch --quiet overleaf "$branch"
project=$(git rev-parse FETCH_HEAD)

pushed=$(git rev-parse --quiet --verify refs/overleaf/pushed || true)
if [ "$project" != "$pushed" ] && ! $force; then
    if [ -z "$pushed" ]; then
        echo "This repository has not pushed to the project yet: run with --force to replace its files." >&2
    else
        echo "The project has changes made in Overleaf since the last sync:" >&2
        git log --oneline "$pushed..$project" >&2 || true
        echo "Carry them into tools/paper/ if they should stay, then run with --force." >&2
    fi
    exit 1
fi

tree=$(git rev-parse HEAD:paper)
if [ "$(git rev-parse "$project^{tree}")" = "$tree" ]; then
    git update-ref refs/overleaf/pushed "$project"
    echo "Overleaf already has paper/ as at $(git rev-parse --short HEAD)."
    exit 0
fi
commit=$(git commit-tree "$tree" -p "$project" -m "pgphase $(git rev-parse --short HEAD): $(git log -1 --format=%s)")
git push --quiet overleaf "$commit:refs/heads/$branch"
git update-ref refs/overleaf/pushed "$commit"
echo "Overleaf now has paper/ as at $(git rev-parse --short HEAD)."
