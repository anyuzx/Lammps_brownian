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
   Changed from fix_nve.cpp May 14th 2016 version. Guang Shi, July 2016
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Brownian Dynamics integrator. Euler Algorithm. 
------------------------------------------------------------------------- */

#include <math.h>
#include "math_extra.h"
#include <stdio.h>
#include <string.h>
#include "fix_bd.h"    
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "update.h"
#include "error.h"
#include "random_mars.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBD::FixBD(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(style,"nve/sphere") != 0 && narg <= 5)
    error->all(FLERR,"Illegal fix nve command");

  t_target = force->numeric(FLERR,arg[3]); // set temperature
  t_period = force->numeric(FLERR,arg[4]); // same as t_period in fix_langevin_overdamp.cpp
  seed = force->inumeric(FLERR,arg[5]); //seed for random number generator. integer

  if (t_target <= 0.0) error->all(FLERR,"Fix bd temperature must be > 0.0");
  if (t_period <= 0.0) error->all(FLERR,"Fix bd period must be > 0.0");
  if (seed <= 0) error->all(FLERR,"Illegal fix bd command");

  // initialize Marsaglia RNG with processor-unique seed

  random = new RanMars(lmp,seed + comm->me);

  dynamic_group_allow = 1;
  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixBD::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBD::init()
{
  dtv = update->dt;  // timestep
  dtf = t_period * force->ftm2v;
  gfactor = sqrt(2*force->boltz*t_target/t_period/dtv/force->mvv2e) / force->ftm2v;
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixBD::initial_integrate(int vflag)
{
  double dtfm;
  double randf;

  // update v and x of atoms in group

  double **x = atom->x;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        randf = sqrt(rmass[i]) * gfactor;
        x[i][0] += dtv * dtfm * (f[i][0]+randf*random->gaussian());
        x[i][1] += dtv * dtfm * (f[i][1]+randf*random->gaussian());
        x[i][2] += dtv * dtfm * (f[i][2]+randf*random->gaussian());
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        randf = sqrt(mass[type[i]]) * gfactor;
        x[i][0] += dtv * dtfm * (f[i][0]+randf*random->gaussian());
        x[i][1] += dtv * dtfm * (f[i][1]+randf*random->gaussian());
        x[i][2] += dtv * dtfm * (f[i][2]+randf*random->gaussian());
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixBD::reset_dt()
{
  dtv = update->dt;
  dtf = t_period * force->ftm2v;
  gfactor = sqrt(2*force->boltz*t_target/t_period/dtv/force->mvv2e) / force->ftm2v;
}
