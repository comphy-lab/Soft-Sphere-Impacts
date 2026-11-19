/**
# Two-phase Bingham (viscoplastic) flow for drop impact

Phase one (`f = 1`) is a Bingham material with plastic viscosity `mu1` and
yield stress `tauy`; phase two is Newtonian with viscosity `mu2`. The header
replaces `twoPhase.h`: it keeps filtered properties and momentum-conserving
VOF advection, and changes only the face viscosity of phase one.

## Constitutive law

The yield criterion is the one used by the Saramito
elasto-viscoplastic model in
[MultiRheoFlow](https://github.com/comphy-lab/MultiRheoFlow):
$\mathcal{F} = \max\left(0, (\|\boldsymbol{\tau}_d\| - \tau_y)/\|\boldsymbol{\tau}_d\|\right)$
on the deviatoric stress $\boldsymbol{\tau}_d$, with
$\|\mathbf{X}\| = \sqrt{\mathbf{X}:\mathbf{X}/2}$. In the inelastic
limit (relaxation time $\lambda \to 0$ at fixed plastic viscosity
$\mu_p = G\lambda$), $\mathcal{F}\boldsymbol{\tau} = 2\mu_p\mathbf{D}$
gives the Bingham law
$$
\boldsymbol{\tau} = 2\left(\mu_p + \frac{\tau_y}{2\|\mathbf{D}\|}\right)\mathbf{D}
\quad\text{for}\quad \|\boldsymbol{\tau}\| > \tau_y,
\qquad \mathbf{D} = 0 \quad\text{otherwise}.
$$
The unyielded branch is regularised by capping the effective viscosity:
$$
\mu_{\text{eff}} = \min\left(\mu_p + \frac{\tau_y}{2\|\mathbf{D}\|},\ \mu_{\max}\right).
$$
Unyielded regions therefore move as a highly viscous fluid rather than a
rigid body. Results should be checked for insensitivity to `mumax`.

The axisymmetric strain-rate norm uses $x$ as the axial and $y$ as the
radial coordinate:
$$
\|\mathbf{D}\| = \sqrt{\tfrac{1}{2}\left(D_{xx}^2 + D_{yy}^2 +
D_{\theta\theta}^2 + 2D_{xy}^2\right)}, \qquad D_{\theta\theta} = u_y/y.
$$
It is evaluated on each face from the neighbouring cell-centred
velocities. Phase averaging remains arithmetic, as in `two-phase.h`.

## Parameters

- `tauy`: yield stress of phase one. `tauy = 0` recovers `twoPhase.h`.
- `mumax`: viscosity cap, which must exceed `mu1` when `tauy > 0`.
*/

#ifndef SOFT_SPHERE_TWO_PHASE_VP_H
#define SOFT_SPHERE_TWO_PHASE_VP_H

#if !AXI
# error "twoPhase-VP.h evaluates the axisymmetric strain rate; include axi.h first."
#endif

#include "vof.h"

scalar f[], * interfaces = {f};

double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0.;
double tauy = 0., mumax = 0.;

face vector alphav[];
scalar rhov[];
scalar sf[];

event defaults (i = 0)
{
  if (tauy > 0. && mumax <= mu1) {
    fprintf (ferr, "twoPhase-VP.h: set mumax > mu1 when tauy > 0.\n");
    exit (1);
  }
  alpha = alphav;
  rho = rhov;
  mu = new face vector;
  reset ((scalar *){mu}, 0);
}

#define rho(f) (clamp(f,0.,1.)*(rho1 - rho2) + rho2)

/**
Effective viscosity of phase one for a given strain-rate norm. */

static inline double bingham_viscosity (double Dnorm)
{
  if (tauy <= 0.)
    return mu1;
  if (Dnorm <= 0.)
    return mumax;
  return min (mu1 + tauy/(2.*Dnorm), mumax);
}

static inline double strain_rate_norm (double Dxx, double Dyy,
                                       double Dqq, double Dxy)
{
  return sqrt ((sq(Dxx) + sq(Dyy) + sq(Dqq) + 2.*sq(Dxy))/2.);
}

/**
The filtered volume fraction smooths the density and viscosity jumps. */

event tracer_advection (i++)
{
  foreach()
    sf[] = (4.*f[] +
            2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
            f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
#if TREE
  set_prolongation (sf, refine_bilinear);
#endif
}

#include "fractions.h"

event properties (i++)
{
  face vector muv = mu;

  /**
  Faces normal to the axial direction $x$. */

  foreach_face (x) {
    double ff = (sf[] + sf[-1])/2.;
    alphav.x[] = fm.x[]/rho(ff);
    double Dxx = (u.x[] - u.x[-1])/Delta;
    double Dyy = (u.y[0,1] - u.y[0,-1] + u.y[-1,1] - u.y[-1,-1])/(4.*Delta);
    double Dqq = (u.y[] + u.y[-1])/(2.*y);
    double Dxy = 0.5*((u.x[0,1] - u.x[0,-1] + u.x[-1,1] - u.x[-1,-1])/(4.*Delta)
                      + (u.y[] - u.y[-1])/Delta);
    double muVP = bingham_viscosity (strain_rate_norm (Dxx, Dyy, Dqq, Dxy));
    muv.x[] = fm.x[]*(clamp(ff,0.,1.)*(muVP - mu2) + mu2);
  }

  /**
  Faces normal to the radial direction $y$. On the axis, $u_y/y$ tends
  to $\partial_y u_y$. */

  foreach_face (y) {
    double ff = (sf[] + sf[0,-1])/2.;
    alphav.y[] = fm.y[]/rho(ff);
    double Dyy = (u.y[] - u.y[0,-1])/Delta;
    double Dxx = (u.x[1,0] - u.x[-1,0] + u.x[1,-1] - u.x[-1,-1])/(4.*Delta);
    double Dqq = y > 0. ? (u.y[] + u.y[0,-1])/(2.*y) : Dyy;
    double Dxy = 0.5*((u.x[] - u.x[0,-1])/Delta
                      + (u.y[1,0] - u.y[-1,0] + u.y[1,-1] - u.y[-1,-1])/(4.*Delta));
    double muVP = bingham_viscosity (strain_rate_norm (Dxx, Dyy, Dqq, Dxy));
    muv.y[] = fm.y[]*(clamp(ff,0.,1.)*(muVP - mu2) + mu2);
  }

  foreach()
    rhov[] = cm[]*rho(sf[]);

#if TREE
  set_prolongation (sf, fraction_refine);
#endif
}

#include "navier-stokes/conserving.h"

#endif // SOFT_SPHERE_TWO_PHASE_VP_H
