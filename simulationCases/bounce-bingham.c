/**
# Axisymmetric impact of a Bingham drop

A unit-radius drop strikes a no-slip wall with unit speed. Length, velocity
and stress are scaled by $R$, $V$ and $\rho V^2$. Pass `We` and `J` as
positional arguments, where $We=\rho V^2R/\sigma$ and $J=\tau_y/(\rho V^2)$.
The plastic Ohnesorge number is $0.01$ and the drop viscosity in solver units
is $Oh/\sqrt{We}$ before adding the regularised yield contribution. `J=0`
recovers the Newtonian constitutive law.

The driver restores `dump` when present. It writes pressure-bearing snapshots
to `intermediate/`, updates `dump` and logs kinetic energy and drop
centre-of-mass velocity to `log`.
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "twoPhase-VP.h"
#include "tension.h"
#include <errno.h>

#define fErr (1e-3)
#define VelErr (1e-3)
#define KErr (1e-3)
#define Rho21 (1e-3)
#define Xdist (1.02)
#define R2Drop(x,y) (sq(x - Xdist) + sq(y))
#define tsnap (1e-2)
#define tlog (5e-4)

u.t[left] = dirichlet(0.);
f[left] = dirichlet(0.);

static const int MAXlevel = 10;
static const double tmax = 1., Ohd = 1e-2, Ohs = 1e-4, Ldomain = 4.;
static double We;

static bool parse_number (const char *arg, double *value) {
  char *end;
  errno = 0;
  *value = strtod (arg, &end);
  return end != arg && *end == '\0' && errno == 0 && isfinite(*value);
}

static void write_parameters (FILE *fp) {
  fprintf (fp, "Level %d tmax %g. We %g, Ohd %.2e, Ohs %.2e, J %g, m %g, Lo %g\n",
           MAXlevel, tmax, We, Ohd, Ohs, bingham_yield_stress,
           bingham_regularisation, Ldomain);
}

int main (int argc, char const *argv[]) {
  if (argc != 3 || !parse_number (argv[1], &We) || We <= 0. ||
      !parse_number (argv[2], &bingham_yield_stress) ||
      bingham_yield_stress < 0.) {
    fprintf (stderr, "Usage: %s We J (finite We > 0 and J >= 0)\n", argv[0]);
    return 1;
  }

  if (system ("mkdir -p intermediate") != 0)
    return 1;

  write_parameters (ferr);
  size (Ldomain);
  origin (0., 0.);
  init_grid (1 << 4);

  rho1 = 1.; mu1 = Ohd/sqrt(We);
  rho2 = Rho21; mu2 = Ohs/sqrt(We);
  f.sigma = 1./We;

  run();
}

event init (t = 0) {
  if (!restore (file = "dump")) {
    refine ((R2Drop(x,y) < 1.05) && (level < MAXlevel));
    fraction (f, 1. - R2Drop(x,y));
    foreach() {
      u.x[] = -f[];
      u.y[] = 0.;
    }
  }
}

scalar KAPPA[];
event adapt (i++) {
  curvature (f, KAPPA);
  adapt_wavelet ((scalar *){f, u.x, u.y, KAPPA},
                 (double[]){fErr, VelErr, VelErr, KErr}, MAXlevel, MAXlevel - 4);
}

event writingFiles (t = 0; t += tsnap; t <= tmax) {
  p.nodump = false;
  dump (file = "dump");
  char nameOut[80];
  snprintf (nameOut, sizeof nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

event logWriting (t = 0; t += tlog; t <= tmax) {
  double ke = 0., momentum = 0., vol = 0.;
  foreach (reduction(+:ke) reduction(+:momentum) reduction(+:vol)) {
    double dv_axi = 2.*pi*y*sq(Delta);
    ke += 0.5*rho(f[])*(sq(u.x[]) + sq(u.y[]))*dv_axi;
    momentum += f[]*u.x[]*dv_axi;
    vol += f[]*dv_axi;
  }

  if (pid() == 0) {
    FILE *fp = fopen ("log", i == 0 ? "w" : "a");
    if (!fp) {
      perror ("log");
      exit (EXIT_FAILURE);
    }
    if (i == 0) {
      write_parameters (fp);
      fputs ("i dt t ke Vcm\n", fp);
      fputs ("i dt t ke Vcm\n", ferr);
    }
    fprintf (fp, "%d %.12g %.12g %.12g %.12g\n", i, dt, t, ke, momentum/vol);
    fclose (fp);
    fprintf (ferr, "%d %.12g %.12g %.12g %.12g\n", i, dt, t, ke, momentum/vol);
  }
}
