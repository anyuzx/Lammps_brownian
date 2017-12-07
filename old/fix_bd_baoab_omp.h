/* -*- c++ -*- ----------------------------------------------------------
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
   Changed from fix_nve_omp.h May 14th 2016 version. Guang Shi, July 2016
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Brownian Dynamics integrator. Euler Algorithm.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(bd/baoab/omp,FixBDBAOABOMP)

#else

#ifndef LMP_FIX_BD_BAOAB_OMP_H
#define LMP_FIX_BD_BAOAB_OMP_H

#include "fix_bd_baoab.h"

namespace LAMMPS_NS {

class FixBDBAOABOMP : public FixBDBAOAB {
 public:
  FixBDBAOABOMP(class LAMMPS *, int, char **);
  virtual ~FixBDBAOABOMP() {}
  virtual void initial_integrate(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
