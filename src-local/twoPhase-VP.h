/**
# Two-phase Bingham (Saramito) elasto-viscoplastic flow for drop impact

Adds a von Mises yield criterion to the log-conformation update of
[twoPhase-ve.h](twoPhase-ve.h). The extra stress obeys the Saramito (2007)
model, which follows the Bingham yield criterion. In continuum form,
$$
\lambda \overset{\nabla}{\boldsymbol{\tau}}_p
+ \mathcal{F}\,\boldsymbol{\tau}_p = 2 G_p \lambda \mathbf{D},
\qquad
\mathcal{F} = \max\left(0, \frac{\|\boldsymbol{\tau}_d\| - \tau_y}
{\|\boldsymbol{\tau}_d\|}\right),
$$
where $\|\boldsymbol{\tau}_d\| = \sqrt{\boldsymbol{\tau}_d :
\boldsymbol{\tau}_d / 2}$ is the von Mises norm of the deviatoric stress.
Below yield ($\mathcal{F} = 0$) the material is a Kelvin–Voigt solid with
modulus $G_p$. Above yield it relaxes like an Oldroyd-B fluid with plastic
viscosity $G_p\lambda$. When $G_p \to \infty$ at fixed $G_p\lambda$, the
model reduces to a Bingham fluid. $\tau_y = 0$ gives $\mathcal{F} = 1$,
which recovers Oldroyd-B exactly.

The yield factor is the `n = 1` case of the elasto-viscoplastic headers in
[MultiRheoFlow](https://github.com/comphy-lab/MultiRheoFlow). As there, it
uses the stress from the previous time step and is applied through the
analytic relaxation factor $\exp(-\Delta t\,\mathcal{F}/\lambda)$.

The driver sets the fields `Gp`, `lambda` and `tau_y`, and includes this
header in place of `twoPhase-ve.h`. This header cannot be
combined with `INFINITE_RELAXATION`, because the yield factor acts through
$\lambda$.
*/

#ifndef SOFT_SPHERE_TWO_PHASE_VP_H
#define SOFT_SPHERE_TWO_PHASE_VP_H

#ifdef INFINITE_RELAXATION
# error "twoPhase-VP.h requires a finite relaxation time"
#endif
#ifdef SOFT_SPHERE_TWO_PHASE_VE_H
# error "include twoPhase-VP.h instead of, not after, twoPhase-ve.h"
#endif

(const) scalar tau_y = zeroc;

/**
The implementation regularises the denominator,
$\mathcal{F} = \max[0, (\|\boldsymbol{\tau}_d\| - \tau_y)/
(\|\boldsymbol{\tau}_d\| + \epsilon)]$, with $\epsilon$ = `yield_eps`.
This matters only where $\|\boldsymbol{\tau}_d\|$ is comparable to
$\epsilon$. */

double yield_eps = 1e-6;

/**
`saramito_factor()` returns $\mathcal{F}$ from the stress components. In
the axisymmetric case the azimuthal component enters the deviatoric norm. */

static inline double saramito_factor (double txx, double tyy, double tqq,
                                      double txy, double ty)
{
  if (ty <= 0.)
    return 1.;
#if AXI
  double tau_d = sqrt ((sq(txx - tyy) + sq(tyy - tqq) + sq(tqq - txx))/6.
                       + sq(txy));
#else
  double tau_d = sqrt (sq(txx - tyy)/4. + sq(txy));
#endif
  return max (0., (tau_d - ty)/(tau_d + yield_eps));
}

#if AXI
# define RELAXATION_RATE_FACTOR(point)                                  \
  saramito_factor (tau_p.x.x[], tau_p.y.y[], tau_qq[], tau_p.x.y[],     \
                   tau_y[])
#else
# define RELAXATION_RATE_FACTOR(point)                                  \
  saramito_factor (tau_p.x.x[], tau_p.y.y[], 0., tau_p.x.y[], tau_y[])
#endif

#include "twoPhase-ve.h"

#endif // SOFT_SPHERE_TWO_PHASE_VP_H
