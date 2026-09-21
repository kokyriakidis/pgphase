#!/usr/bin/env bash
# Regenerate src/test_gap_windows_expect.tsv from a measured run.
#
# Runs the window-test binary in emit mode, so the numbers are produced by the
# same code that asserts them. Do this only when a measured improvement is
# intended, and say in the commit which arm moved and why -- the file is the
# record of what the pipeline is expected to achieve on these windows, and
# refreshing it to make a red test green is how a regression gets committed.
set -euo pipefail
readonly REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
readonly OUT="${OUT:-${REPO}/src/test_gap_windows_expect.tsv}"
cd "${REPO}"
[[ -x ./test_gap_windows ]] || make test_gap_windows
# The binary writes the file itself; Catch2 captures stdout, so piping it would
# yield an indented, unparseable copy.
PGPHASE_EMIT_EXPECTATIONS="${OUT}.tmp" ./test_gap_windows "[windows]" > /dev/null
# The active graph arm must be present, or the file would silently install
# without its assertions. The name tracks the arms in src/test_gap_windows.cpp.
grep -q $'^graph\t' "${OUT}.tmp"
mv "${OUT}.tmp" "${OUT}"
echo "wrote ${OUT}: $(grep -vc '^#' "${OUT}") rows"
