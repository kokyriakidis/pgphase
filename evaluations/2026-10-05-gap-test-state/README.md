# Persistent gap test replay state

The native gap harness now shares successful pipeline output state across test
processes and shards. Every invocation still executes its assertions and parental
truth scoring; test verdicts are not cached. Failed or incomplete replays are not
published. File locks prevent duplicate concurrent production of the same replay.

Keys include the executable content hash, helper content hash, normalized command,
input/index and dynamic-library file identities (path, device, inode, size, mtime,
ctime), and relevant environment. Output identities are checked before reuse.
Changing production code requires a fresh replay after rebuilding. Set
`PGPHASE_TEST_CACHE` to an empty directory to request cold verification.

Four shards covered all 86 cases and 12,264 assertions in both runs. The first
run took 548.21 seconds; the warm run took 31.94 seconds (17.2 times faster).
The initial cache already contained the two focused measurements; it otherwise
needed the full panel's replay state. All 114 states and their timestamps remained
unchanged on the warm run, so it launched zero new pipeline replays.

The new-gap focused case ran 747 assertions in 15.945 seconds cold and 0.280
seconds warm. The phase-matrix diagnostic case took 16.296 seconds cold and
0.189 seconds warm, including restoration of its diagnostic matrix products.

Four helper tests cover cross-process reuse, input/index/binary/output
invalidation, failed replay retries, concurrent ownership, quoted paths, and
matrix restoration. The production build and all C++ unit tests passed. No
production phasing behavior or test expectations changed for this task.

See validation.json and the copied logs for measured shard results.

The normal unsharded `make window-tests` also passed all 86 cases (11,997
assertions) and all four helper tests using the populated default cache, whose
114 states remained intact. Assertion counts differ with process boundaries
because the existing harness shares previously measured outcomes within a
process. `make check` and all 47 predicate cases / 1,536 assertions passed.
