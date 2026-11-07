# multiVencISMRM2026

## Branches

- **Working branch: `dev-inflow`.** Work and PRs target it.
- **`dev` is frozen** as the reproducible record: never edit it. It stays the GitHub default branch, so
  do not base work on the default; branch from `dev-inflow`.
- **`dev-multiVencISMRM2026` and `dev-postISMRM` are also frozen.** Do not edit them.
- **`errata` and `scratch-setup` are forward branches.** They are kept in step with `dev-inflow` by
  separate PRs (one PR per branch for a change that all three need).

## Tools

- A tool pinned to a non-default branch is cloned into the project's own `devTools/<tool>`
  (gitignored real clones, set up by the doIts), never into the shared `/scratch/bass/tools/<tool>`.
  Never switch the branch of a shared tools clone.
