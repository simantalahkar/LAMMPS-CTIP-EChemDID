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

#ifdef FIX_CLASS

FixStyle(qeq/ctip,FixQEqCtip)

#else

#ifndef LMP_FIX_QEQ_CTIP_H
#define LMP_FIX_QEQ_CTIP_H

#include "fix_qeq.h"

namespace LAMMPS_NS {

class FixQEqCtip : public FixQEq {
 public:
  FixQEqCtip(class LAMMPS *, int, char **);
  ~FixQEqCtip();
  void init();
  void pre_force(int);

  double memory_usage();
  virtual int pack_forward_comm(int, int *, double *, int, int *);
  virtual void unpack_forward_comm(int, int, double *);

 private:
  void conjugate_gradient();
  double brent_solver(double *);
  double calculate_qf(double,double *);
  double dUdq_self(double,int,int);
  double dUdq_two_body(int,int,double);
  double dUdq_two_body_ewald(double,double,double,double,double);
  double dUdq_two_body_wolf(double,double,double,double,double);
  void extract_ctip();
  double *qmin, *qmax;  // extract from pair_coul_ctip
  double ctip_omega, wolfewald, cut_coul; // extract from pair_coul_ctip
  int kspacetype; // extract from pair_coul_ctip
  
  int get_names_ctip(char *,double *&); // EChemDID					   

  double brent_maxiter, brent_critical, brent_tolerance; // for brent solver

  class Pair *pair;
};
}
#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix qeq/ctip requires atom attribute q

Self-explanatory.

E: Fix qeq/ctip group has no atoms

Self-explanatory.

E: Invalid param file for fix qeq/ctip

Zeta value is 0.0.

E: No pair coul/ctip for fix qeq/ctip

These commands must be used together.

E: Fix qeq/ctip could not extract params from pair coul/ctip

This should not happen unless pair coul/ctip has been altered.

W: Charges did not converge at step ...

Occurs when maximum number of conjugate gradient iterations are reached.
Increase maxiter value in the fix_qeq command.

E:  Brent solver limit reached

Occurs when maximum number of brent solver iterations are reached.
Increase brent_maxiter keyword.

*/
