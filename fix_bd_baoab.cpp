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
#include "fix_bd_baoab.h"
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

FixBDBAOAB::FixBDBAOAB(LAMMPS *lmp, int narg, char **arg) :
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

  rnum = NULL;

  nvalues = 3;
  grow_arrays(atom->nmax);
  atom->add_callback(0);

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++){
    rnum[i][0] = random->gaussian();
    rnum[i][1] = random->gaussian();
    rnum[i][2] = random->gaussian();
  }
}

/* ---------------------------------------------------------------------- */

int FixBDBAOAB::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBDBAOAB::init()
{
  dtv = update->dt;  // timestep
  dtf = t_period * force->ftm2v;
  gfactor = sqrt(0.5*force->boltz*t_target/t_period/dtv/force->mvv2e) / force->ftm2v;
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixBDBAOAB::initial_integrate(int vflag)
{
  double dtfm;
  double randf;
  double rx, ry, rz;

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
        rx = random->gaussian();
        ry = random->gaussian();
        rz = random->gaussian();
        x[i][0] += dtv * dtfm * (f[i][0]+randf*(rx + rnum[i][0]));
        x[i][1] += dtv * dtfm * (f[i][1]+randf*(ry + rnum[i][1]));
        x[i][2] += dtv * dtfm * (f[i][2]+randf*(rz + rnum[i][2]));
        rnum[i][0] = rx;
        rnum[i][1] = ry;
        rnum[i][2] = rz;
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        randf = sqrt(mass[type[i]]) * gfactor;
        rx = random->gaussian();
        ry = random->gaussian();
        rz = random->gaussian();
        x[i][0] += dtv * dtfm * (f[i][0]+randf*(rx + rnum[i][0]));
        x[i][1] += dtv * dtfm * (f[i][1]+randf*(ry + rnum[i][1]));
        x[i][2] += dtv * dtfm * (f[i][2]+randf*(rz + rnum[i][2]));
        rnum[i][0] = rx;
        rnum[i][1] = ry;
        rnum[i][2] = rz;
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixBDBAOAB::reset_dt()
{
  dtv = update->dt;
  dtf = t_period * force->ftm2v;
  gfactor = sqrt(0.5*force->boltz*t_target/t_period/dtv/force->mvv2e) / force->ftm2v;
}

/* ----------------------------------------------------------------------
   allocate atom-based array for rnum
------------------------------------------------------------------------- */

void FixBDBAOAB::grow_arrays(int nmax)
{
  memory->grow(rnum,nmax,3,"fix_bd_baoab:rnum");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixBDBAOAB::copy_arrays(int i, int j, int delflag)
{
  for (int m = 0; m < nvalues; m++)
    rnum[j][m] = rnum[i][m];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixBDBAOAB::pack_exchange(int i, double *buf)
{
  for (int m = 0; m < nvalues; m++) buf[m] = rnum[i][m];
  return nvalues;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixBDBAOAB::unpack_exchange(int nlocal, double *buf)
{
  for (int m = 0; m < nvalues; m++) rnum[nlocal][m] = buf[m];
  return nvalues;
}
