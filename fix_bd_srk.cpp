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
#include "fix_bd_srk.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "update.h"
#include "error.h"
#include "memory.h"
#include "random_mars.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBDSRK::FixBDSRK(LAMMPS *lmp, int narg, char **arg) :
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

  fprev = NULL;

  nvalues = 3;
  grow_arrays(atom->nmax);
  atom->add_callback(0);

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++){
    fprev[i][0] = 0.0;
    fprev[i][1] = 0.0;
    fprev[i][2] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */
/*
FixBDSRK::~FixBDSRK()
{
  delete random;
  memory->destroy(fprev);
  atom->delete_callback(id,0);
}
*/
/* ---------------------------------------------------------------------- */

int FixBDSRK::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBDSRK::init()
{
  dtv = update->dt;  // timestep
  dtf = t_period * force->ftm2v;
  gfactor = sqrt(2*force->boltz*t_target/t_period/dtv/force->mvv2e) / force->ftm2v;
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixBDSRK::initial_integrate(int vflag)
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
        fprev[i][0] = f[i][0];
        fprev[i][1] = f[i][1];
        fprev[i][2] = f[i][2];
        x[i][0] += dtv * dtfm * (f[i][0] + randf*random->gaussian());
        x[i][1] += dtv * dtfm * (f[i][1] + randf*random->gaussian());
        x[i][2] += dtv * dtfm * (f[i][2] + randf*random->gaussian());
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        randf = sqrt(mass[type[i]]) * gfactor;
        fprev[i][0] = f[i][0];
        fprev[i][1] = f[i][1];
        fprev[i][2] = f[i][2];
        x[i][0] += dtv * dtfm * (f[i][0] + randf*random->gaussian());
        x[i][1] += dtv * dtfm * (f[i][1] + randf*random->gaussian());
        x[i][2] += dtv * dtfm * (f[i][2] + randf*random->gaussian());
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixBDSRK::final_integrate()
{
  double dtfm;

  // update x of atoms in group

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
        x[i][0] += dtv * 0.5 * dtfm * (f[i][0] - fprev[i][0]);
        x[i][1] += dtv * 0.5 * dtfm * (f[i][1] - fprev[i][1]);
        x[i][2] += dtv * 0.5 * dtfm * (f[i][2] - fprev[i][2]);
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        x[i][0] += dtv * 0.5 * dtfm * (f[i][0] - fprev[i][0]);
        x[i][1] += dtv * 0.5 * dtfm * (f[i][1] - fprev[i][1]);
        x[i][2] += dtv * 0.5 * dtfm * (f[i][2] - fprev[i][2]);
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixBDSRK::reset_dt()
{
  dtv = update->dt;
  dtf = t_period * force->ftm2v;
  gfactor = sqrt(2*force->boltz*t_target/t_period/dtv/force->mvv2e) / force->ftm2v;
}

/* ----------------------------------------------------------------------
   allocate atom-based array for fprev
------------------------------------------------------------------------- */

void FixBDSRK::grow_arrays(int nmax)
{
  memory->grow(fprev,nmax,3,"fix_bd_srk:fprev");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixBDSRK::copy_arrays(int i, int j, int delflag)
{
  for (int m = 0; m < nvalues; m++)
    fprev[j][m] = fprev[i][m];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixBDSRK::pack_exchange(int i, double *buf)
{
  for (int m = 0; m < nvalues; m++) buf[m] = fprev[i][m];
  return nvalues;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixBDSRK::unpack_exchange(int nlocal, double *buf)
{
  for (int m = 0; m < nvalues; m++) fprev[nlocal][m] = buf[m];
  return nvalues;
}
