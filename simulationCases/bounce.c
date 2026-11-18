/**
 * @file bounce.c
 * @author Vatsal Sanjay (vatsalsanjay@gmail.com)
 * vatsalsanjay.com
 * Physics of Fluids
 * @date 2023-05-30
 *
*/

// 1 is drop
#include "axi.h"
#include "navier-stokes/centered.h"
#include "twoPhase.h"
#include "tension.h"

// error tolerances
#define fErr (1e-3)
#define VelErr (1e-3)
#define KErr (1e-3)

// air-water
#define Rho21 (1e-3)
// Calculations!
#define Xdist (1.02)
#define R2Drop(x,y) (sq(x - Xdist) + sq(y))

// boundary conditions
u.t[left] = dirichlet(0.0);
f[left] = dirichlet(0.0);

int MAXlevel;
double tmax, We, Ohd, Ohs, Ec, De, Bo, Ldomain;
#define tsnap (0.01)

int main(int argc, char const *argv[]) {

  if (argc != 2) {
    fprintf (stderr, "Usage: %s We\n", argv[0]);
    return 1;
  }

  MAXlevel = 10;
  We = atof(argv[1]);
  tmax = 1.0;

  Ohd = 1e-2; // <\mu/sqrt(1250*0.060*0.001)>
  Ohs = 1e-4; //\mu_r * Ohd

  Ldomain = 4.0; // size of domain. must keep Ldomain \gg 1

  fprintf(ferr, "Level %d tmax %g. We %g, Ohd %3.2e, Ohs %3.2e, Ec %g, De %g, Lo %g\n", MAXlevel, tmax, We, Ohd, Ohs, Ec, De, Ldomain);

  L0=Ldomain;
  X0=0.; Y0=0.;
  init_grid (1 << (4));

  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  rho1 = 1.0; mu1 = Ohd/sqrt(We);
  rho2 = Rho21; mu2 = Ohs/sqrt(We);
  f.sigma = 1.0/We;

  run();
}

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

scalar KAPPA[];
event adapt(i++) {
  curvature(f, KAPPA);
  adapt_wavelet ((scalar *){f, u.x, u.y, KAPPA},
    (double[]){fErr, VelErr, VelErr, KErr},
    MAXlevel, MAXlevel-4);
}

// Outputs
// static
event writingFiles (t = 0, t += tsnap; t <= tmax) {
   p.nodump = false; // dump pressure to calculate force in post-processing
  dump (file = "dump");
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

event logWriting (i++) {
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
      fprintf(fp, "Level %d tmax %g. We %g, Ohd %3.2e, Ohs %3.2e, Ec %g, De %g, Lo %g\n", MAXlevel, tmax, We, Ohd, Ohs, Ec, De, Ldomain);
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
