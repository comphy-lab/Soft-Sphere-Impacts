/**
# Kelvin–Voigt sphere impacting a non-wetting wall

An axisymmetric elastic sphere of radius $R_0$ hits a rigid wall with
velocity $V_0$ through a surrounding gas. The sphere is a Kelvin–Voigt
solid: an elastic network of modulus $G$ in parallel with a viscous
solvent. It is the infinite-relaxation (elastic-memory, $Wi \to \infty$)
limit of the Oldroyd-B drop in `bounce-ve.c`. The wall force follows
from the centre-of-mass velocity of the sphere.

## Scaling

Lengths are scaled by $R_0$, velocities by $V_0$ and densities by the
sphere density $\rho_l$. Time is therefore in units of $R_0/V_0$,
stresses in units of $\rho_l V_0^2$ and forces in units of
$\rho_l V_0^2 R_0^2$.

- `We` $= \rho_l V_0^2 R_0/\gamma$: Weber number; surface tension is $1/We$.
- `Ec` $= G R_0/\gamma$: elastocapillary number.
- $El = Ec/We = G/(\rho_l V_0^2)$: elasticity number, the dimensionless
  modulus applied inside the sphere.
- $Oh_d = 10^{-2}$ and $Oh_s = 10^{-4}$: Ohnesorge numbers of the sphere
  solvent and the gas; viscosities are $Oh/\sqrt{We}$.
- Gas-to-sphere density ratio $10^{-3}$.

## Constitutive model

`twoPhase-ve.h` evolves the conformation tensor $\mathbf{A}$ with the
log-conformation method and adds the elastic stress
$\boldsymbol{\tau}_p = G(\mathbf{A} - \mathbf{I})$ to the momentum
equation. Defining `INFINITE_RELAXATION` removes the relaxation term, so
$\mathbf{A}$ records the full deformation history and the sphere recovers
its shape elastically.

## Usage

```sh
./bounce-kv We Ec [MAXlevel]
```

`MAXlevel` is the maximum refinement level (default 10); the finest cell
size is $L_0/2^{MAXlevel}$ with $L_0 = 4$. The run restarts from `dump` if
the file exists in the working directory.

## Author
Vatsal Sanjay
Email: vatsal.sanjay@comphy-lab.org
CoMPhy Lab
*/

// Phase 1 (f = 1) is the sphere.
#include "axi.h"
#include "navier-stokes/centered.h"
#define INFINITE_RELAXATION
#include "twoPhase-ve.h"
#include "tension.h"

/**
## Numerical parameters

Wavelet-adaptation tolerances for the volume fraction, velocity,
curvature and conformation components. */

#define fErr (1e-3)
#define VelErr (1e-3)
#define KErr (1e-3)
#define AErr (1e-3)

// Gas-to-sphere density ratio.
#define Rho21 (1e-3)
// Initial centre position: a gap of 0.02 separates the sphere from the wall.
#define Xdist (1.02)
#define R2Drop(x,y) (sq(x - Xdist) + sq(y))

/**
## Boundary conditions

The axis of symmetry is the bottom boundary and the wall is the left
boundary. The wall is impermeable (default) and no-slip. Setting $f = 0$
there keeps the sphere from wetting it, so a thin gas film always
separates the two. The top and right boundaries keep the default
symmetry conditions. */

u.t[left] = dirichlet(0.0);
f[left] = dirichlet(0.0);

int MAXlevel;
double tmax, We, Ohd, Ohs, Ec, De, Bo, Ldomain;
#define tsnap (1e-2) // snapshot interval
#define tlog (5e-4)  // log interval for the force signal

scalar Gpd[]; // elastic modulus field

int main(int argc, char const *argv[]) {

  if (argc != 3 && argc != 4) {
    fprintf (stderr, "Usage: %s We Ec [MAXlevel]\n", argv[0]);
    return 1;
  }

  MAXlevel = (argc == 4) ? atoi(argv[3]) : 10;
  We = atof(argv[1]);
  tmax = 1.0;

  Ohd = 1e-2; // sphere solvent
  Ohs = 1e-4; // gas

  Ldomain = 4.0; // domain size; must be much larger than the sphere

  Ec = atof(argv[2]);

  fprintf(ferr, "Level %d tmax %g. We %g, Ohd %3.2e, Ohs %3.2e, Ec %g, De infty, Lo %g\n", MAXlevel, tmax, We, Ohd, Ohs, Ec, Ldomain);

  L0=Ldomain;
  X0=0.; Y0=0.;
  init_grid (1 << (4));

  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  rho1 = 1.0; mu1 = Ohd/sqrt(We);
  rho2 = Rho21; mu2 = Ohs/sqrt(We);
  f.sigma = 1.0/We;

  Gp = Gpd;

  run();
}

/**
## Elastic modulus

The modulus is $El = Ec/We$ wherever the sphere is present and zero in
the gas. */

 event properties (i++) {
   foreach (){
     Gpd[] = (f[] < 1e-6) ? 0.0: Ec/We;
   }
 }

/**
## Initial condition

A unit sphere centred at $x = 1.02$ moves towards the wall with unit
velocity; the gas starts at rest. The grid is refined to `MAXlevel`
around the sphere. */

event init(t = 0){
  if(!restore (file = "dump")){
    refine((R2Drop(x,y) < 1.05) && (level < MAXlevel));
    fraction (f, 1. - R2Drop(x,y));
    foreach () {
      u.x[] = -1.0*f[];
      u.y[] = 0.0;
    }
  }
}

/**
## Adaptation

The grid adapts between levels `MAXlevel - 4` and `MAXlevel` on the
interface, velocity, curvature and conformation fields. */

scalar KAPPA[];
event adapt(i++) {
  curvature(f, KAPPA);
   adapt_wavelet ((scalar *){f, u.x, u.y, KAPPA, conform_p.x.x, conform_p.x.y, conform_p.y.y, conform_qq},
     (double[]){fErr, VelErr, VelErr, KErr, AErr, AErr, AErr, AErr},
     MAXlevel, MAXlevel-4);

}

/**
## Outputs

Every `tsnap`, the run writes the restart file `dump` and a snapshot in
`intermediate/`, both including pressure. */

event writingFiles (t = 0, t += tsnap; t <= tmax) {
   p.nodump = false;
  dump (file = "dump");
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

/**
Every `tlog`, `log` records the kinetic energy of both phases and the
centre-of-mass velocity of the sphere,
$$
V_{cm} = \frac{\int f\,u_x\,\mathrm{d}V}{\int f\,\mathrm{d}V},
\qquad \mathrm{d}V = 2\pi y\,\mathrm{d}x\,\mathrm{d}y.
$$
The wall force follows from the momentum balance of the sphere,
$F = (4\pi/3)\,\mathrm{d}V_{cm}/\mathrm{d}t$ in units of
$\rho_l V_0^2 R_0^2$ (`postProcess/f_vcm_kv.py`). */

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
