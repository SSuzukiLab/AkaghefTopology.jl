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
sh lib/topology.jl/test/client.sh --instant --exec .research-scratch/one-call.jl
sh lib/topology.jl/test/client.sh --long --limit 45 --exec .research-scratch/acceptance-run.jl
```

### Tiers of the runner

A request goes to one of three tiers by how long it may run (Akaghef,
2026-11-25). The limit is enforced: a request that passes it is stopped, its
log ends with a `RUNNER: stopped after ...` line, and its client exits 124.

| Tier | Argument | Limit | Julia | Use |
| --- | --- | --- | --- | --- |
| instant | `--instant` | 15 s | one shared server (`test/.queue/instant`) | One call on a minimal example, one contract example. It is not delayed by the work tier. |
| work | none | 10 min | one shared server (`test/.queue/work`) | Unit tests of the touched modules, index generation, an experiment on a small case. |
| long | `--long --limit <minutes>` | the declared limit | its own Julia under `test/.queue/long/`, which exits when the request ends | A run of 10 minutes or more, that is, an acceptance check. `--batch` is the old name and means `--long --limit 60`. |

Each shared server runs one request at a time. A supervisor process beside each
server measures every request from the moment the server takes it. At the
limit it stops that Julia process (the exact pid in `server.pid`), writes the
request's exit status, and releases the queue; the next request starts a fresh
server, and requests that were still waiting are resubmitted to it by their
own clients. A request stopped at the limit of the instant tier belongs in the
work tier; one stopped in the work tier needs a smaller case or a split
question, and only an acceptance run goes to the long tier.

With `--exec`, the script runs in a fresh Julia module from the caller's working
directory. `ARGS`, `PROGRAM_FILE` and relative `include` work as in a standalone
script. The client returns the script's output and an error status if it throws.
Use standalone Julia for scripts that call `exit`, change process-wide settings,
or require a different project or startup flags. Use local, trusted queues only.

The first request to a queue starts its server and waits for it; a start after
source edits recompiles `Topology` and takes a few minutes. If a request was
not dispatched (Revise rejected a source edit, or the request ahead of it was
stopped at the limit), its client resubmits it to the next server, up to five
times. `TOPOLOGY_TEST_WAIT` bounds the client's total wait in seconds (default:
600 for the instant tier, 3600 for the work tier, the limit plus 900 for the
long tier). A client that times out leaves its request in the queue. For a
slow request the client reports whether it is waiting for startup, queue
service, execution setup or a result. Each server generation has its own
`server.out` under its queue. Set `TOPOLOGY_JULIA` to the Julia executable path
if `julia` on `PATH` is unavailable in the launching shell.

The supervisor touches `server.super` in the generation's directory every two
seconds, also while Julia is busy. A client that finds a queue locked by a
server whose `server.super` has been silent for 20 seconds (checked twice, five
seconds apart) treats the server as abandoned: it stops the Julia process
recorded in `server.pid` if that is still running, releases the lock and
starts a new server. No manual cleanup of `server.lock` is needed.

`TOPOLOGY_TEST_QUEUE` replaces the queue path when a lane needs independent
execution, for example
`TOPOLOGY_TEST_QUEUE=.research-scratch/my-lane/julia-queue sh lib/topology.jl/test/client.sh --exec .research-scratch/my-experiment.jl`.
Each queue starts its own Julia process with the limit of the tier named on
the command line that started it. Run a script that may not terminate (new
code, a diagnosis) on such a queue, so that stopping it restarts no shared
server. `TOPOLOGY_TEST_CAP` sets the limit in seconds of a server the call
starts; it exists for testing the runner.

The queues `test/.queue/interactive` and `test/.queue/batch` of the earlier
layout are no longer used.

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
