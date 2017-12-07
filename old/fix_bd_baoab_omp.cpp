/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Changed from fix_nve_omp.cpp May 14th 2016 version. Guang Shi, July 2016
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Brownian Dynamics integrator. Euler Algorithm.
------------------------------------------------------------------------- */

#include <math.h>
#include "math_extra.h"
#include "fix_bd_baoab_omp.h"
#include "atom.h"
#include "force.h"
#include "random_mars.h"

using namespace LAMMPS_NS;
using namespace FixConst;

typedef struct { double x,y,z; } dbl3_t;

/* ---------------------------------------------------------------------- */

FixBDBAOABOMP::FixBDBAOABOMP(LAMMPS *lmp, int narg, char **arg) :
  FixBDBAOAB(lmp, narg, arg) { }

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixBDBAOABOMP::initial_integrate(int vflag)
{
  // update v and x of atoms in group

  dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  const dbl3_t * _noalias const f = (dbl3_t *) atom->f[0];
  dbl3_t * _noalias const rn = (dbl3_t *) rnum[0];
  const int * const mask = atom->mask;
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst : atom->nlocal;
  int i;

  if (atom->rmass) {
    const double * const rmass = atom->rmass;
#if defined (_OPENMP)
#pragma omp parallel for private(i) default(none) schedule(static)
#endif
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        const double dtfm = dtf / rmass[i];
        const double randf = sqrt(rmass[i]) * gfactor;
        const double rx = random->gaussian();
        const double ry = random->gaussian();
        const double rz = random->gaussian();
        x[i].x += dtv * dtfm * (f[i].x+randf*(rx + rn[i].x));
        x[i].y += dtv * dtfm * (f[i].y+randf*(ry + rn[i].y));
        x[i].z += dtv * dtfm * (f[i].z+randf*(rz + rn[i].z));
        rn[i].x = rx;
        rn[i].y = ry;
        rn[i].z = rz;
      }

  } else {
    const double * const mass = atom->mass;
    const int * const type = atom->type;
#if defined (_OPENMP)
#pragma omp parallel for private(i) default(none) schedule(static)
#endif
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        const double dtfm = dtf / mass[type[i]];
        const double randf = sqrt(mass[type[i]]) * gfactor;
        const double rx = random->gaussian();
        const double ry = random->gaussian();
        const double rz = random->gaussian();
        x[i].x += dtv * dtfm * (f[i].x+randf*(rx + rn[i].x));
        x[i].y += dtv * dtfm * (f[i].y+randf*(ry + rn[i].y));
        x[i].z += dtv * dtfm * (f[i].z+randf*(rz + rn[i].z));
        rn[i].x = rx;
        rn[i].y = ry;
        rn[i].z = rz;
      }
  }
}
