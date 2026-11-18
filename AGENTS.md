# Agent instructions

Read `README.md` for the case inputs and build commands.

## Repository rules

- Keep Basilisk headers in `src-local/`, the three standalone impact drivers
  in `simulationCases/` and snapshot analysis in `postProcess/`.
- Preserve each driver's positional arguments, boundary conditions, event
  order and restart format when moving code. The direct `El Wi` build uses
  `bounce-ve.c`. The Kelvin–Voigt build uses the infinite-relaxation limit.
- Build against the project-local Basilisk release pinned in `basilisk.ref`.
  Keep the installation, build artefacts and executables out of Git.
- Write runs and analysis outputs outside the repository. Keep generated
  results, figures, run lists, data-root maps, run directory identities and
  machine-specific paths out of commits.
- Leave the force-scaling notebook and its paired input data together.
  A published figure identifies the solver commit and Basilisk ref used to
  generate it.
- Use the repository's Git author identity. Commit messages carry no
  tool-generated attribution. Keep a local `CLAUDE.md` untracked and point
  it to this file.
