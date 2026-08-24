# Math Harness canonical boundary

ADR source: `20260813T054438Z-windows-codex-inko-message-math-harness-personal-work-package-adr-001`

Processed canonical packet:
`/Users/nisimoriyuuya/Akaghef-Bridge/done/windows_to_mac/20260813T054438Z-windows-codex-inko-message-math-harness-personal-work-package-adr-001.json`

Packet SHA-256:
`dba79202fcea05133ef91703b3720eb79077f37febd440bd810486f0d62166e0`

The Git repository rooted at `topology/` (origin
`https://github.com/SSuzukiLab/AkaghefTopology`) is the independent,
self-authored library. `julia/` is its Julia migration workbench; it is not a
separate repository and is not itself a portable Math Harness work package.
The Windows-local Math Harness at `C:\Users\Akaghef\data\21_研究` remains the
authority for mathematical execution evidence, claim status, and promotion
into canonical research results.

The A-sys shared session owns coordination and boundary enforcement for this
migration. That coordination ownership does not replace the Windows Harness as
the mathematical execution and evidence authority.

The boundary applied to this migration is:

1. Library source remains in the independent `AkaghefTopology` repository.
   A Windows Harness must consume a clean, accepted repository revision by
   pinned tag or commit; neither the repository nor the current worktree may be
   copied wholesale into the Harness.
2. `Project.toml` and `Manifest.toml` pin the Julia dependency environment.
   MATLAB baseline images and replay records in `artifacts/` are migration
   evidence, not Windows Harness verification or a promoted mathematical claim.
3. Any `.qmd`, rendered HTML/PDF, figures, small input data, or annotations sent
   to the personal/iPad surface are portable work-package material only.
4. iPad/personal edits return as an explicit packet or import candidate. They
   never overwrite this repository or the Windows Harness automatically.
5. Before a returned edit changes a canonical Harness file, ledger entry, claim
   status, or result, Windows must import it explicitly, run the pinned Julia
   environment, compare the required structural/numerical/plot evidence, and
   review the result.

## Current migration gate

- Mac-side MATLAB replay, Julia tests, and manual plot inspection may establish
  `mac-replayed` evidence only.
- The current migration worktree contains uncommitted source and generated
  artifacts. It is not an acceptable Harness dependency revision until the
  intended library contents are reviewed and recorded in a clean tag/commit.
- Completion of the Mac migration does not promote a mathematical result.
  Promotion requires the Windows Harness to consume that exact revision,
  recreate the pinned environment, rerun the required workflows, and record
  the Windows-side structural, numerical, and visual review.
- A personal/iPad edit is input to this process only through an explicit return
  packet or import candidate. No sync path may write it directly into either
  this repository or the Windows Harness.

## Required portable-package manifest

A package cut from this work must state at least:

- source repository URL and exact commit/tag;
- Julia and library versions (`Manifest.toml` or its digest);
- source Live Script(s) and baseline identifiers;
- exact run/test command;
- creation date;
- evidence status: `unverified`, `mac-replayed`, or `windows-verified`;
- for each plotted `VirtualLink`, source `REdgeTable.Position` polylines and
  crossing positions, so geometric Julia/MATLAB comparison is reproducible;
- explicit return/import packet identifier for edited material.

No artifact in this directory currently carries `windows-verified` status.
This ADR introduces no contradiction that requires stopping the Mac migration;
it adds a hard promotion boundary after Mac-side acceptance and before any
Windows canonical claim or ledger update.

## Finding logged 2026-08-20, corrected 2026-08-21: the visual clause needs a
## stateless `SageWrapper` on Windows

This section records a fact found after the ADR was written.  **It does not
amend the boundary.**  Changing the promotion conditions is the A-sys / ADR
owner's decision, not this repository's.  What follows is input to that
decision.

Windows-side logs supplied on 2026-08-20
(`sage-windows-matlab-logs-20260820.zip`; ChatGPT session
`260221_Git_submodule_commit_error`, Codex thread
`019d43a5-275b-7ed0-aff5-4050d18aba55`) establish that Windows MATLAB cannot
reach a WSL Sage through `pyenv`: `pyenv` embeds a Windows-native CPython into
the MATLAB process, while WSL's Sage is a Linux ELF build.  `SW/SageWrapper.m`
depends on that embedding and on a persistent `sage.all.__dict__` namespace, so
it has no Windows form.  An external-process fallback through `wsl.exe` was
attempted and did not land.  Details in
[`PLOT_PIPELINE.md`](PLOT_PIPELINE.md) section 7.

**Correction, 2026-08-21.**  The paragraph above stated the constraint too
strongly, and the consequences below are the corrected version.  Sage *is*
reachable from Windows MATLAB: `Execution/experiment/sage/run_sage4.m`
(2026-06-23) drives WSL Sage over an HTTP bridge and renders a plot back into
MATLAB.  What has no Windows form is `SageWrapper`'s persistent `sage.all`
namespace, because the bridge runs `sage -c` as a fresh process per call.

Consequences for the conditions stated above:

1. The requirement that Windows rerun and compare structural, numerical **and
   visual** evidence is **conditionally satisfiable**.  It requires
   `SageWrapper` to be made stateless — prepending the construction expression
   to each `exec`, which `setSageLink` already builds.  Thirteen methods in
   `VirtualLink.m` reach Sage.  Until that rewrite, no Windows `VirtualLink`
   workflow runs.
2. The manifest item "for each plotted `VirtualLink`, source
   `REdgeTable.Position` polylines and crossing positions" is likewise
   producible on Windows after that rewrite.  Mac-side export is possible now.
3. The Julia path has no Python, Sage or system dependency and is expected to
   run natively on Windows, which would satisfy a plot-evidence clause without
   the rewrite.  **This has not been verified on Windows.**

Nothing here relaxes the boundary.  Until the ADR owner rules, plot evidence
produced anywhere remains `mac-replayed` at best, and no promotion follows from
this repository's Mac-side work.
