/**
# Bingham (Saramito) elasto-viscoplastic sphere impacting a wall

An axisymmetric sphere of radius $R_0$ hits a rigid, non-wetting wall at
speed $V_0$. Below its yield stress the sphere is a Kelvin–Voigt solid. Above
it, the sphere flows with the Saramito model in
[twoPhase-VP.h](../src-local/twoPhase-VP.h). The geometry, solver and
outputs match [bounce-ve.c](bounce-ve.c), at maximum level 10 ($\Delta =
R_0/256$) like [bounce-kv.c](bounce-kv.c).

Solver units: $\rho_l = R_0 = V_0 = 1$. The drop and gas viscosities are
$Oh/\sqrt{We}$ with $Oh = 10^{-2}$ and $10^{-4}$, and $\gamma = 1/We$.

## Usage

```
bounce-bingham We Ec De J
```

- `We = rho V0^2 R0/gamma`, the Weber number.
- `Ec = G R0/gamma`, the elastocapillary number. The modulus is `Ec/We`.
- `De = lambda/sqrt(rho R0^3/gamma)`, the Deborah number. The relaxation time
  is `sqrt(We)*De`, so the Weissenberg number is `Wi = De*sqrt(We)`.
- `J = tau_y R0/gamma`, the plasticity number. The yield stress is `J/We`.

All four arguments must be finite numbers; `We`, `Ec` and `De` positive and
`J` non-negative.

The plastic viscosity is $G\lambda$, and the Bingham number is
$Bi = \tau_y R_0/(G\lambda V_0) = J/(Ec\,De\sqrt{We})$.

## Outputs

In the working directory: `log` (a parameter line, the column header
`i dt t ke Vcm`, then one row every `tlog`), `dump` (the restart file), and
`intermediate/snapshot-<t>` every `tsnap`. A run restores from `dump` when it
exists and drops log rows written after the restored time. `postProcess/pforce.sh`
has no Bingham case; the impact force follows from `Vcm` as in
[bounce-kv.c](bounce-kv.c).
*/

// 1 is drop
#include "axi.h"
#include "navier-stokes/centered.h"
#include "twoPhase-VP.h"
#include "tension.h"
#include <errno.h>

// error tolerances
#define fErr (1e-3)
#define VelErr (1e-3)
#define KErr (1e-3)
#define AErr (1e-3)

// air-water
#define Rho21 (1e-3)
// Calculations!
#define Xdist (1.02)
#define R2Drop(x,y) (sq(x - Xdist) + sq(y))

// boundary conditions
u.t[left] = dirichlet(0.0);
f[left] = dirichlet(0.0);

int MAXlevel;
double tmax, We, Ohd, Ohs, Ec, De, J, Ldomain;
#define tsnap (1e-2)
#define tlog (5e-4)

scalar Gpd[], lambdav[], tau_yv[];
static double Gdrop, lambdaDrop, tauyDrop;

static double parse_arg (const char * s, const char * name)
{
  char * end;
  errno = 0;
  double v = strtod (s, &end);
  if (end == s || *end != '\0' || errno == ERANGE || !isfinite (v)) {
    fprintf (stderr, "Invalid %s: %s\n", name, s);
    exit (1);
  }
  return v;
}

static FILE * open_log (const char * mode)
{
  FILE * fp = fopen ("log", mode);
  if (!fp) {
    perror ("log");
    exit (1);
  }
  return fp;
}

static void close_log (FILE * fp)
{
  int failed = ferror (fp);
  failed |= fclose (fp);
  if (failed) {
    fprintf (stderr, "log: write failed\n");
    exit (1);
  }
}

/**
On restart, keep the two header lines and the complete rows up to the
restored time; `restore()` does not re-fire events due at that time, so
the run appends from the next row. It reads the global `t`: inside an event
`t` is the event's own argument, which `restore()` leaves unchanged. */

static void trim_log (void)
{
  double t0 = t;
  FILE * in = fopen ("log", "r");
  if (!in)
    return;
  FILE * out = fopen ("log.tmp", "w");
  if (!out) {
    perror ("log.tmp");
    exit (1);
  }
  char line[512];
  int n = 0;
  while (fgets (line, sizeof (line), in)) {
    if (!strchr (line, '\n'))
      break; // partial last line
    int it;
    double dti, ti;
    if (n++ < 2 || (sscanf (line, "%d %lf %lf", &it, &dti, &ti) == 3 &&
                    ti <= t0 + 1e-9))
      fputs (line, out);
  }
  int failed = ferror (in);
  fclose (in);
  close_log (out);
  if (failed || rename ("log.tmp", "log")) {
    fprintf (stderr, "log: restart trim failed\n");
    exit (1);
  }
}

static void require_finite (double v, const char * name, int positive)
{
  if (!isfinite (v) || (positive && v <= 0.)) {
    fprintf (stderr, "Derived %s out of range: %g\n", name, v);
    exit (1);
  }
}

int main(int argc, char const *argv[]) {

  if (argc != 5) {
    fprintf (stderr, "Usage: %s We Ec De J\n", argv[0]);
    return 1;
  }

  MAXlevel = 10;
  We = parse_arg (argv[1], "We");
  tmax = 1.0;

  Ohd = 1e-2;
  Ohs = 1e-4;

  Ldomain = 4.0; // size of domain. must keep Ldomain \gg 1

  Ec = parse_arg (argv[2], "Ec");
  De = parse_arg (argv[3], "De");
  J = parse_arg (argv[4], "J");
  if (We <= 0. || Ec <= 0. || De <= 0. || J < 0.) {
    fprintf (stderr, "Need We, Ec, De > 0 and J >= 0.\n");
    return 1;
  }
  Gdrop = Ec/We;
  lambdaDrop = sqrt(We)*De;
  tauyDrop = J/We;
  require_finite (Gdrop, "modulus Ec/We", 1);
  require_finite (lambdaDrop, "relaxation time sqrt(We)*De", 1);
  require_finite (tauyDrop, "yield stress J/We", J > 0.);
  require_finite (1./We, "surface tension 1/We", 1);
  require_finite (1./sqrt(We), "viscosity scale 1/sqrt(We)", 1);

  fprintf(ferr, "Level %d tmax %g. We %g, Ohd %3.2e, Ohs %3.2e, Ec %g, De %g, J %g, Lo %g\n", MAXlevel, tmax, We, Ohd, Ohs, Ec, De, J, Ldomain);

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
  lambda = lambdav;
  tau_y = tau_yv;

  run();
}

/**
Material properties act only inside the drop. They are refreshed from the
volume fraction after the VOF step and before the constitutive update: this
driver's `tracer_advection` runs before the one in
[twoPhase-ve.h](../src-local/twoPhase-ve.h). A cell the drop has just
entered therefore keeps its strain and gets the drop modulus at once. */

event tracer_advection (i++) {
  foreach () {
    double inside = (f[] < 1e-6) ? 0. : 1.;
    Gpd[] = inside*Gdrop;
    lambdav[] = inside*lambdaDrop;
    tau_yv[] = inside*tauyDrop;
  }
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
  else if (pid() == 0)
    trim_log ();
}

scalar KAPPA[];
event adapt(i++) {
  curvature(f, KAPPA);
   adapt_wavelet ((scalar *){f, u.x, u.y, KAPPA, conform_p.x.x, conform_p.x.y, conform_p.y.y, conform_qq},
     (double[]){fErr, VelErr, VelErr, KErr, AErr, AErr, AErr, AErr},
     MAXlevel, MAXlevel-4);

}

// Outputs
event writingFiles (t = 0, t += tsnap; t <= tmax) {
  p.nodump = false; // dump pressure to calculate force in post-processing
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

  if (pid() == 0){
    if (i == 0) {
      fprintf (ferr, "i dt t ke Vcm\n");
      fp = open_log ("w");
      fprintf(fp, "Level %d tmax %g. We %g, Ohd %3.2e, Ohs %3.2e, Ec %g, De %g, J %g, Lo %g\n", MAXlevel, tmax, We, Ohd, Ohs, Ec, De, J, Ldomain);
      fprintf (fp, "i dt t ke Vcm\n");
      fprintf (fp, "%d %g %g %g %g\n", i, dt, t, ke, Vcm);
      close_log (fp);
    } else {
      fp = open_log ("a");
      fprintf (fp, "%d %g %g %g %g\n", i, dt, t, ke, Vcm);
      close_log (fp);
    }
    fprintf (ferr, "%d %g %g %g %g\n", i, dt, t, ke, Vcm);
  }

}
