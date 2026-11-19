# Soft Sphere Impacts

## Problem

An axisymmetric sphere impacts a solid wall. The four drivers describe a
Newtonian drop, a viscoelastic drop with finite relaxation time, a
Kelvin–Voigt sphere in the infinite-relaxation limit and a Bingham
(viscoplastic) drop.

## Layout

```
src-local/twoPhase.h - filtered two-phase flow and momentum-conserving advection
src-local/twoPhase-ve.h - log-conformation update and stress coupling
src-local/twoPhase-VP.h - two-phase flow with a regularised Bingham drop
simulationCases/bounce.c - Newtonian drop
simulationCases/bounce-ve.c - finite-relaxation viscoelastic drop
simulationCases/bounce-kv.c - Kelvin–Voigt sphere
simulationCases/bounce-bingham.c - Bingham drop
simulationCases/runSimulation.sh - compile and run a driver
postProcess/ - snapshot extraction, force, footprint and video tools
Force-data-scalings-elastic/ - existing force-scaling notebook and its input data
```

## Method

The drivers use Basilisk's centred Navier–Stokes solver on an adaptive
axisymmetric grid. The volume fraction tracks the interface. Mesh adaptation
uses volume fraction, velocity and curvature, with conformation components
added for the viscoelastic and Kelvin–Voigt cases. The Bingham drop
uses the inelastic limit of the Saramito yield criterion, with the
effective viscosity capped at $10^3$ times the plastic viscosity in
unyielded regions. Each driver can restore
from `dump` in its working directory.

Reproducing a published figure requires rerunning the simulations with the
solver commit and the release in `basilisk.ref` that produced it.

## Running

From the repository root, download and inspect the installer at the pinned
release, then install Basilisk into the checkout.

```sh
tag=$(cat basilisk.ref)
curl -fsSLO "https://raw.githubusercontent.com/comphy-lab/basilisk-C/$tag/reset_install_basilisk-ref-locked.sh"
# Inspect the downloaded script before running it.
bash reset_install_basilisk-ref-locked.sh --ref="$tag"
```

From a separate output directory for each run, call the runner using its
absolute path. It passes the positional arguments directly to the driver.

```sh
/path/to/Soft-Sphere-Impacts/simulationCases/runSimulation.sh bounce 20
/path/to/Soft-Sphere-Impacts/simulationCases/runSimulation.sh bounce-ve 20 0.1 1.0
/path/to/Soft-Sphere-Impacts/simulationCases/runSimulation.sh bounce-kv 20 0.1
/path/to/Soft-Sphere-Impacts/simulationCases/runSimulation.sh bounce-bingham 100 0.001
/path/to/Soft-Sphere-Impacts/simulationCases/runSimulation.sh bounce-ve-wi 100 0.1 0.1
```

The drivers take `We` (Weber number), followed by `Ec` (elastocapillary
number) for Kelvin–Voigt or `Ec De` (Deborah number last) for finite relaxation.
Both viscoelastic drivers use `Ec/We` for the modulus. The Bingham
driver takes `We J`, with `J = tau_y/(rho V^2)` the dimensionless yield
stress. `bounce-kv` and `bounce-bingham` accept an optional final
`MAXlevel` argument (default 10). The finite-relaxation
driver uses `sqrt(We)*De` for the relaxation time in solver units. The `bounce-ve-wi` runner builds the same finite-relaxation
source with `BOUNCE_VE_EL_WI` and takes the dimensionless modulus
`El = G/(rho*V^2)` and relaxation time `Wi = lambda*V/R` directly. The runner defaults to a serial build. Set `OPENMP=1` with a
compiler that supports `-fopenmp` to use the eight threads of the earlier
launchers. Override the count with `OMP_NUM_THREADS`.

Run the post-processing scripts from the simulation output directory. They
read `intermediate/` and `log`. Force extraction replaces the per-run
pressure and footprint tables and appends to case-specific summary files in
the parent directory (for example, `force_max_bounce-ve.dat`). The Python tools use NumPy, pandas, Matplotlib and
SciPy. Video generation also uses LaTeX and FFmpeg.

```sh
/path/to/Soft-Sphere-Impacts/postProcess/pforce.sh bounce 20
/path/to/Soft-Sphere-Impacts/postProcess/video.sh 4 512 20 0.01
```

For force extraction, `bounce-ve` takes `We`, `bounce-kv` takes `We Ec`,
`bounce-ve-wi` takes `We El` and `bounce-bingham` takes `We J`. A final video argument sets the output
filename, for example `elastic_We20.mp4` for the Kelvin–Voigt case.
