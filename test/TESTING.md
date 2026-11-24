# Test policy for Topology.jl

Tests are split by module and by tier. A change is checked by the unit tests of
the modules it touches, not by the whole suite. The whole suite runs once a day.

## Tiers

| Tier | What belongs here | Who runs it, when |
| --- | --- | --- |
| unit | Small fixed examples of one module. A file should finish within 60 s of wall time after `using Topology`. | Whoever changes the module, on every change, for the touched modules only. |
| integration | Chains across modules on small examples (HC to Koda to triangulation to invariant), Live Script ports, database round trips. | The lead, before judging a delegated turn or committing. |
| acceptance | Corpus sweeps, random walks, inverse walks, census runs, classification campaigns. | Once a day, unattended, with the log kept. |
| manual | One-off campaigns and files that no runner included before 2026-11-17. | Only when named by path. |

Random walks and exhaustive searches are experiments. They never belong in the
unit tier, and a unit test file must not `include` one.

## Running

From `lib/topology.jl`:

```sh
julia --project=. test/run.jl heegaard.koda                 # unit tests of one module
julia --project=. test/run.jl --tier integration heegaard   # unit + integration of a module family
julia --project=. test/run.jl --tier acceptance              # the daily full run
julia --project=. test/run.jl --list --tier unit             # what a selection would run
julia --project=. test/run.jl --check                        # every test file is in the manifest
```

Use the persistent runner for repeated tests and experiments. It keeps `Topology`
loaded under Revise and serves one request at a time per queue. Run these commands
from the repository root in a POSIX shell (Git Bash on Windows):

```sh
sh lib/topology.jl/test/client.sh heegaard.koda              # same arguments as run.jl
sh lib/topology.jl/test/client.sh --exec .research-scratch/my-experiment.jl 'argument with spaces'
```

With `--exec`, the script runs in a fresh Julia module from the caller's working
directory. `ARGS`, `PROGRAM_FILE` and relative `include` work as in a standalone
script. The client returns the script's output and an error status if it throws.
Use standalone Julia for scripts that call `exit`, change process-wide settings,
or require a different project or startup flags. Use local, trusted queues only.

The first request starts the background server and waits for it. If Revise
rejects a source edit before dispatch, the client retries once on a fresh
server. `TOPOLOGY_TEST_WAIT` bounds total client wait time (seconds, default
3600). A timed-out request may still run. Each server generation has its own
`server.out` under the queue (`test/.queue` by default).
Set `TOPOLOGY_JULIA` to the Julia executable path if `julia` on `PATH` is
unavailable in the launching shell.

Share the default queue for short requests. Use `TOPOLOGY_TEST_QUEUE` with a
dedicated directory when a lane needs independent execution, for example
`TOPOLOGY_TEST_QUEUE=.research-scratch/my-lane/julia-queue sh lib/topology.jl/test/client.sh --exec .research-scratch/my-experiment.jl`.
Each queue starts its own Julia process. If the shell supervising a server is
forcibly stopped, verify the server is gone using the generation's `server.pid`
before removing the stale `server.lock` directory or switching queues.

The acceptance tier runs every day at 03:00 as the Windows scheduled task
"Topology acceptance (daily)" (`test/acceptance.ps1`, started late if the machine
was asleep). Logs are in `test/acceptance-logs/` (not tracked): `<date>.log`,
`timings.tsv`, and `latest.txt` with the last result line.

`Pkg.test()` and `test/heegaard/runtests.jl` call the same runner at the tier named
by `TOPOLOGY_TEST_TIER` (default `unit`). Add `--pkgimages=no` where the README
requires it on this Windows machine.

## The manifest

`test/manifest.toml` lists every test file with its suite, its module tags and its
tier. A module argument matches a tag exactly or as a dotted prefix
(`heegaard` selects `heegaard.koda`). `--check` fails when a test file is missing
from the manifest, so a new test file must be added there with a tier.

Each file's wall time is appended to `test/timings.tsv` (not tracked). A file in
the unit tier that exceeds 60 s is moved to integration, and the reason is
written in its `note`. Tiers marked `provisional` were assigned by name and are
reclassified from measured timings.

## Known costs (2026-11-17)

- Loading `Topology` after any source edit recompiles the package: about 180 s.
  With a valid cache it takes about 9 s. Most of the first minute of each small
  test file is first-call compilation.
- `heegaard/koda_joining_prepare.jl` takes about 500 s for three small diagrams;
  the time is in `hc_prepare_joining_disk` itself, not in the test.
- `heegaard/lc_moves_inverse_walk.jl` (random inverse walks) takes about 95 min.
  It is acceptance only.
