/**
# Bingham drop impacting a non-wetting wall

An axisymmetric Bingham (viscoplastic) sphere of radius $R_0$ hits a
rigid, non-wetting wall with velocity $V_0$. The geometry, gas phase,
Ohnesorge numbers, adaptation and force logging match `bounce-kv.c`; only
the drop rheology differs.

Lengths, velocities and densities are scaled by $R_0$, $V_0$ and the drop
density $\rho_l$, so time is in units of $R_0/V_0$ and stresses in units
of $\rho_l V_0^2$.

## Parameters

- `We` $= \rho_l V_0^2 R_0/\gamma$: Weber number.
- `J` $= \tau_y/(\rho_l V_0^2)$: dimensionless yield stress.
- `MAXlevel` (optional, default 10): maximum refinement level.

The plastic viscosity is $\mu_p = Oh_d/\sqrt{We}$ with $Oh_d = 10^{-2}$,
so the Bingham number is $Bi = \tau_y R_0/(\mu_p V_0) = J\sqrt{We}/Oh_d$.
The unyielded branch is regularised with
$\mu_{\max} = 10^3\mu_p$ (see `twoPhase-VP.h`).

## Usage

```sh
./bounce-bingham We J [MAXlevel]
```
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "twoPhase-VP.h"
#include "tension.h"

// Adaptation tolerances: volume fraction, velocity and curvature.
#define fErr (1e-3)
#define VelErr (1e-3)
#define KErr (1e-3)

// Gas-to-drop density ratio.
#define Rho21 (1e-3)
// Initial centre position: a gap of 0.02 separates the sphere from the wall.
#define Xdist (1.02)
#define R2Drop(x,y) (sq(x - Xdist) + sq(y))
// Viscosity cap relative to the plastic viscosity.
#define MUMAX_RATIO (1e3)

/**
The wall is the left boundary: no slip and $f = 0$, so the drop never
wets it. The axis is the bottom boundary. */

u.t[left] = dirichlet(0.0);
f[left] = dirichlet(0.0);

int MAXlevel;
double tmax, We, Ohd, Ohs, J, Ldomain;
#define tsnap (1e-2)
#define tlog (5e-4)

int main (int argc, char const *argv[]) {

  if (argc != 3 && argc != 4) {
    fprintf (stderr, "Usage: %s We J [MAXlevel]\n", argv[0]);
    return 1;
  }

  We = atof(argv[1]);
  J = atof(argv[2]);
  MAXlevel = (argc == 4) ? atoi(argv[3]) : 10;
  tmax = 1.0;

  Ohd = 1e-2; // drop Ohnesorge number (plastic viscosity)
  Ohs = 1e-4; // gas Ohnesorge number

  Ldomain = 4.0;

  fprintf (ferr, "Level %d tmax %g. We %g, Ohd %3.2e, Ohs %3.2e, J %g, Lo %g\n",
           MAXlevel, tmax, We, Ohd, Ohs, J, Ldomain);

  L0 = Ldomain;
  X0 = 0.; Y0 = 0.;
  init_grid (1 << (4));

  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system (comm);

  rho1 = 1.0; mu1 = Ohd/sqrt(We);
  rho2 = Rho21; mu2 = Ohs/sqrt(We);
  f.sigma = 1.0/We;

  tauy = J;
  mumax = MUMAX_RATIO*mu1;

  run();
}

event init (t = 0) {
  if (!restore (file = "dump")) {
    refine ((R2Drop(x,y) < 1.05) && (level < MAXlevel));
    fraction (f, 1. - R2Drop(x,y));
    foreach () {
      u.x[] = -1.0*f[];
      u.y[] = 0.0;
    }
  }
}

scalar KAPPA[];
event adapt (i++) {
  curvature (f, KAPPA);
  adapt_wavelet ((scalar *){f, u.x, u.y, KAPPA},
    (double[]){fErr, VelErr, VelErr, KErr},
    MAXlevel, MAXlevel-4);
}

/**
## Outputs

Snapshots every `tsnap` keep the pressure for force post-processing. The
log records the kinetic energy and the drop's centre-of-mass velocity
$V_{cm}$ every `tlog`; the wall force follows from
$F = (4\pi/3)\,\mathrm{d}V_{cm}/\mathrm{d}t$. */

event writingFiles (t = 0, t += tsnap; t <= tmax) {
  p.nodump = false;
  dump (file = "dump");
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

event logWriting (t = 0, t += tlog; t <= tmax) {
  double ke = 0., Vcm = 0., vol = 0.;
  foreach (reduction(+:ke) reduction(+:Vcm) reduction(+:vol)) {
    ke += 2*pi*y*(0.5*rho(f[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
    Vcm += 2*pi*y*((f[])*(u.x[]))*sq(Delta);
    vol += 2*pi*y*(f[])*sq(Delta);
  }
  Vcm /= vol;

  static FILE * fp;

  if (pid() == 0) {
    if (i == 0) {
      fprintf (ferr, "i dt t ke Vcm\n");
      fp = fopen ("log", "w");
      fprintf (fp, "Level %d tmax %g. We %g, Ohd %3.2e, Ohs %3.2e, J %g, Lo %g\n",
               MAXlevel, tmax, We, Ohd, Ohs, J, Ldomain);
      fprintf (fp, "i dt t ke Vcm\n");
      fprintf (fp, "%d %g %g %g %g\n", i, dt, t, ke, Vcm);
      fclose (fp);
    } else {
      fp = fopen ("log", "a");
      fprintf (fp, "%d %g %g %g %g\n", i, dt, t, ke, Vcm);
      fclose (fp);
    }
    fprintf (ferr, "%d %g %g %g %g\n", i, dt, t, ke, Vcm);
  }
}
