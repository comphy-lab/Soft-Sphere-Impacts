/**
# Kelvin–Voigt sphere impacting a wall

An axisymmetric sphere of radius $R_0$ hits a rigid, non-wetting wall at
speed $V_0$ through a surrounding gas. The sphere is an Oldroyd-B material in
the infinite-relaxation limit ($De, Wi \to \infty$): a neo-Hookean elastic
solid with modulus $G$ plus a small solvent viscosity, that is, a
Kelvin–Voigt solid. This is the elastic-memory limit of Jana et al.,
"Impacting spheres: from liquid drops to elastic beads", Soft Matter 22,
2226 (2026).

The solver integrates the incompressible two-phase Navier–Stokes equations
with surface tension. The extra stress in the sphere is
$\boldsymbol{\tau}_p = G(\mathbf{A} - \mathbf{I})$, with the conformation
tensor $\mathbf{A}$ advanced by the log-conformation method of
[twoPhase-ve.h](../src-local/twoPhase-ve.h).

Solver units: $\rho_l = R_0 = V_0 = 1$, so time is in $R_0/V_0$ and stress
in $\rho_l V_0^2$.

## Usage

```
bounce-kv We Ec
```

- `We = rho_l V0^2 R0/gamma`, the Weber number. Surface tension is `1/We`.
- `Ec = G R0/gamma`, the elastocapillary number. The modulus in solver units
  is the elasticity number `El = Ec/We = G/(rho_l V0^2)`.

Fixed in the source: sphere Ohnesorge number $Oh = 10^{-2}$, gas
$Oh = 10^{-4}$, gas-to-sphere density ratio $10^{-3}$, domain $4R_0$, end
time $t = 1$ and maximum level 10, which gives a smallest cell
$\Delta = R_0/256$.

## Outputs

Written to the working directory:

- `log`: a parameter line, the column header `i dt t ke Vcm`, then one row
  every `tlog`. `ke` is the total kinetic energy and `Vcm` the sphere's
  centre-of-mass velocity normal to the wall.
- `intermediate/snapshot-<t>`: full dumps, including pressure, every
  `tsnap`.
- `dump`: the latest snapshot. A run restores from it when present.

The impact force is estimated from the sphere's momentum balance,
$F = (4\pi/3)\,\mathrm{d}V_{cm}/\mathrm{d}t$ in units of
$\rho_l V_0^2 R_0^2$ (`postProcess/f_vcm_kv.py`). This is the net axial
force on the sphere, including gas traction. `postProcess/getBetaMax.c`
integrates the wall pressure instead.

Author: Vatsal Sanjay.
*/

// f = 1 is the sphere, f = 0 the gas.
#include "axi.h"
#include "navier-stokes/centered.h"
// No relaxation step: A is only stretched and advected (elastic memory).
#define INFINITE_RELAXATION
#include "twoPhase-ve.h"
#include "tension.h"

// Wavelet error tolerances for adaptation: volume fraction, velocity,
// curvature and conformation components.
#define fErr (1e-3)
#define VelErr (1e-3)
#define KErr (1e-3)
#define AErr (1e-3)

// Gas-to-sphere density ratio.
#define Rho21 (1e-3)
// The axis is the bottom boundary and the wall is the left one (x = 0).
// The sphere centre starts at x = 1.02, leaving a 0.02 gas gap.
#define Xdist (1.02)
#define R2Drop(x,y) (sq(x - Xdist) + sq(y))

// Wall: no slip, and f = 0 imposes a non-wetting wall. Top and right keep
// Basilisk's default free-slip boundaries.
u.t[left] = dirichlet(0.0);
f[left] = dirichlet(0.0);

int MAXlevel;
double tmax, We, Ohd, Ohs, Ec, De, Bo, Ldomain;
#define tsnap (1e-2) // snapshot interval
#define tlog (5e-4)  // log interval, fine enough to differentiate Vcm

scalar Gpd[]; // elastic modulus field

int main(int argc, char const *argv[]) {

  if (argc != 3) {
    fprintf (stderr, "Usage: %s We Ec\n", argv[0]);
    return 1;
  }

  MAXlevel = 10;
  We = atof(argv[1]);
  tmax = 1.0;

  Ohd = 1e-2; // sphere Ohnesorge number
  Ohs = 1e-4; // gas Ohnesorge number

  Ldomain = 4.0; // size of domain. must keep Ldomain \gg 1

  Ec = atof(argv[2]);

  fprintf(ferr, "Level %d tmax %g. We %g, Ohd %3.2e, Ohs %3.2e, Ec %g, De infty, Lo %g\n", MAXlevel, tmax, We, Ohd, Ohs, Ec, Ldomain);

  L0=Ldomain;
  X0=0.; Y0=0.;
  init_grid (1 << (4));

  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  // Inertial scaling: mu = Oh/sqrt(We) and gamma = 1/We.
  rho1 = 1.0; mu1 = Ohd/sqrt(We);
  rho2 = Rho21; mu2 = Ohs/sqrt(We);
  f.sigma = 1.0/We;

  Gp = Gpd;

  run();
}

/**
The modulus $G = Ec/We$ acts only inside the sphere; the gas carries no
elastic stress. */

 event properties (i++) {
   foreach (){
     Gpd[] = (f[] < 1e-6) ? 0.0: Ec/We;
   }
 }

/**
A fresh run refines around the sphere, sets its volume fraction and gives
it the impact velocity $-V_0$ towards the wall. A restart skips this. */

event init(t = 0){
  if(!restore (file = "dump")){
    refine((R2Drop(x,y) < 1.05) && (level < MAXlevel));
    fraction (f, 1. - R2Drop(x,y)); //initialise a drop with radius 1
    foreach () {
      u.x[] = -1.0*f[];
      u.y[] = 0.0;
    }
  }
}

/**
Adaptation also follows the conformation components, so the grid stays
refined where the solid is strained. */

scalar KAPPA[];
event adapt(i++) {
  curvature(f, KAPPA);
   adapt_wavelet ((scalar *){f, u.x, u.y, KAPPA, conform_p.x.x, conform_p.x.y, conform_p.y.y, conform_qq},
     (double[]){fErr, VelErr, VelErr, KErr, AErr, AErr, AErr, AErr},
     MAXlevel, MAXlevel-4);

}

// Outputs
event writingFiles (t = 0, t += tsnap; t <= tmax) {
   p.nodump = false; // keep pressure in dumps for wall-pressure forces
  dump (file = "dump");
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

/**
Volume integrals carry the axisymmetric weight $2\pi y$. */

event logWriting (t = 0, t += tlog; t <= tmax) {
  double ke = 0., Vcm = 0., vol = 0.;
  foreach (reduction(+:ke) reduction(+:Vcm) reduction(+:vol)) {
    ke += 2*pi*y*(0.5*rho(f[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
    Vcm += 2*pi*y*((f[])*(u.x[]))*sq(Delta);
    vol += 2*pi*y*(f[])*sq(Delta);
  }
  Vcm /= vol;

  static FILE * fp;

  if (pid() == 0){
    if (i == 0) {
      fprintf (ferr, "i dt t ke Vcm\n");
      fp = fopen ("log", "w");
      fprintf(fp, "Level %d tmax %g. We %g, Ohd %3.2e, Ohs %3.2e, Ec %g, De infty, Lo %g\n", MAXlevel, tmax, We, Ohd, Ohs, Ec, Ldomain);
      fprintf (fp, "i dt t ke Vcm\n");
      fprintf (fp, "%d %g %g %g %g\n", i, dt, t, ke, Vcm);
      fclose(fp);
    } else {
      fp = fopen ("log", "a");
      fprintf (fp, "%d %g %g %g %g\n", i, dt, t, ke, Vcm);
      fclose(fp);
    }
    fprintf (ferr, "%d %g %g %g %g\n", i, dt, t, ke, Vcm);
  }

}
