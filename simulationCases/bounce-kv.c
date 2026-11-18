/**
# Kelvin–Voigt sphere impact

Simulates an axisymmetric sphere impacting a rigid, non-wetting wall through
air. Viscous stresses and an elastic stress evolved without relaxation
give the Kelvin–Voigt limit of the two-phase viscoelastic model.

Run `bounce-kv We Ec` with the Weber number `We = rho V^2 R/sigma`
and the elastocapillary number `Ec = G R/sigma`. The solver uses `R`, `V`
and `rho` of the sphere as units. The sphere Ohnesorge number is fixed at
`0.01`, and the output time horizon is `t V/R = 1`.

- `dump` is the latest checkpoint and is restored when present on launch.
- `intermediate/snapshot-*` stores time-stamped states, including pressure.
- `log` records iteration, time step, time, total two-phase kinetic energy
  and the sphere's axial centre-of-mass velocity. Its time derivative
  gives the force after multiplication by the sphere mass.

Original driver by Vatsal Sanjay, 30 May 2023.
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#define INFINITE_RELAXATION
#include "twoPhase-ve.h"
#include "tension.h"

// Wavelet tolerances for interface, velocity, curvature and conformation.
#define fErr (1e-3)
#define VelErr (1e-3)
#define KErr (1e-3)
#define AErr (1e-3)

// Ambient-to-sphere density ratio and initial sphere centre.
#define Rho21 (1e-3)
#define Xdist (1.02)
#define R2Drop(x,y) (sq(x - Xdist) + sq(y))

// The left wall is no-slip. Zero wall volume fraction maintains an air film.
u.t[left] = dirichlet(0.0);
f[left] = dirichlet(0.0);

int MAXlevel;
double tmax, We, Ohd, Ohs, Ec, De, Bo, Ldomain;
#define tsnap (1e-2)
#define tlog (5e-4)

scalar Gpd[];

int main(int argc, char const *argv[]) {

  if (argc != 3) {
    fprintf (stderr, "Usage: %s We Ec\n", argv[0]);
    return 1;
  }

  MAXlevel = 10;
  We = atof(argv[1]);
  tmax = 1.0;

  Ohd = 1e-2;
  Ohs = 1e-4; // Uses the sphere's density and surface tension in the scaling.

  Ldomain = 4.0;

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

 // Apply the modulus in sphere and interfacial cells, excluding ambient gas.
 event properties (i++) {
   foreach (){
     Gpd[] = (f[] < 1e-6) ? 0.0: Ec/We;
   }
 }

event init(t = 0){
  if(!restore (file = "dump")){
    // Resolve the initially spherical interface before assigning its velocity.
    refine((R2Drop(x,y) < 1.05) && (level < MAXlevel));
    fraction (f, 1. - R2Drop(x,y));
    foreach () {
      u.x[] = -1.0*f[];
      u.y[] = 0.0;
    }
  }
}

scalar KAPPA[];
event adapt(i++) {
  curvature(f, KAPPA);
   // Include the azimuthal conformation so hoop stresses remain resolved.
   adapt_wavelet ((scalar *){f, u.x, u.y, KAPPA, conform_p.x.x, conform_p.x.y, conform_p.y.y, conform_qq},
     (double[]){fErr, VelErr, VelErr, KErr, AErr, AErr, AErr, AErr},
     MAXlevel, MAXlevel-4);

}

event writingFiles (t = 0, t += tsnap; t <= tmax) {
   p.nodump = false; // Retain wall pressure for post-processing.
  dump (file = "dump");
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

event logWriting (t = 0, t += tlog; t <= tmax) {
  double ke = 0., Vcm = 0., vol = 0.;
  // Integrate with the ring volume 2*pi*y*Delta^2 over the meridional grid.
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
