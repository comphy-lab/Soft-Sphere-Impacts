/**
# Two-phase viscoplastic flow for drop impact

The drop follows the $n=1$ limit of the Saramito–Herschel–Bulkley yield law
in [MultiRheoFlow](https://github.com/comphy-lab/MultiRheoFlow/blob/68222b6ff926971b623cf4e56d29e075b57ea086/src-local/log-conform-elastoviscoplastic-HB-scalar-2D.h)
as $\lambda \to 0$ at fixed $G_p\lambda$. Its ideal yielded stress
magnitude is $\tau_y + \mu_1\dot\gamma$, where $\mu_1$ is the total plastic
viscosity.
Papanastasiou regularisation replaces the singular yield contribution to the
apparent viscosity with $\tau_y(1-e^{-m\dot\gamma})/\dot\gamma$.
Here $\dot\gamma = \sqrt{2\mathbf{D}'{:}\mathbf{D}'}$, with $\mathbf{D}'$
the deviatoric strain-rate tensor. The ambient phase remains Newtonian.
*/

#ifndef SOFT_SPHERE_TWO_PHASE_VP_H
#define SOFT_SPHERE_TWO_PHASE_VP_H

#if !AXI
#error "twoPhase-VP.h requires axisymmetric geometry"
#endif

#include "twoPhase.h"

double bingham_yield_stress = 0.;
/* m is measured in the solver's time units. */
static const double bingham_regularisation = 100.;

event defaults (i = 0) {
  /* The base two-phase module allocates mu when either phase has viscosity. */
  if (bingham_yield_stress > 0. && !mu1 && !mu2)
    mu = new face vector;
}

static void update_bingham_viscosity (void) {
  scalar shear_rate[];
  foreach() {
    double inv_2delta = 0.5/Delta;
    double dxx = (u.x[1,0] - u.x[-1,0])*inv_2delta;
    double drr = (u.y[0,1] - u.y[0,-1])*inv_2delta;
    double dtt = u.y[]/y;
    double dxr = 0.5*(u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0])*
      inv_2delta;
    shear_rate[] = sqrt((2./3.)*(sq(dxx - drr) + sq(drr - dtt) +
                                sq(dtt - dxx)) + 4.*sq(dxr));
  }

#if TREE
  /* Match the two-phase properties event's interpolation at coarse–fine faces. */
  set_prolongation (sf, refine_bilinear);
#endif
  face vector muv = mu;
  foreach_face() {
    double ff = clamp(0.5*(sf[] + sf[-1]), 0., 1.);
    double rate = 0.5*(shear_rate[] + shear_rate[-1]);
    double factor = rate > 0. ? -expm1(-bingham_regularisation*rate)/rate :
      bingham_regularisation;
    muv.x[] = fm.x[]*(mu(ff) + ff*bingham_yield_stress*factor);
  }
#if TREE
  set_prolongation (sf, fraction_refine);
#endif
}

event viscous_term (i++) {
  /* Freeze the apparent viscosity for the centred solver's implicit step. */
  if (bingham_yield_stress > 0.)
    update_bingham_viscosity();
}

#endif
