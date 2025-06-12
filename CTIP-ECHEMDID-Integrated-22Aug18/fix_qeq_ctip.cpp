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
   Contributing author: Ray Shan (Sandia)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fix_qeq_ctip.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "force.h"
#include "group.h"
#include "pair.h"
#include "kspace.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixQEqCtip::FixQEqCtip(LAMMPS *lmp, int narg, char **arg) :
  FixQEq(lmp, narg, arg)
{

  if (narg < 8) error->all(FLERR,"Illegal fix qeq/ctip command");

  if (nevery <= 0 || cutoff <= 0.0 || tolerance <= 0.0 || maxiter <=0.0)
    error->all(FLERR,"Illegal fix qeq/ctip command");

  // optional arg
  brent_maxiter = 100;
  brent_critical = 1.0e-8;
  brent_tolerance = 1.0e-11;
  int iarg = 8;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"brent_maxiter") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix qeq/ctip command");
      brent_maxiter = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"brent_critical") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix qeq/ctip command");
      brent_critical = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"brent_tolerance") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix qeq/ctip command");
      brent_tolerance = atof(arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix qeq/ctip command");
  }
  comm_forward = 1;

}

/* ---------------------------------------------------------------------- */

FixQEqCtip::~FixQEqCtip()
{
  memory->destroy(chi);
  memory->destroy(eta);
  memory->destroy(gamma);
  memory->destroy(zeta);
  memory->destroy(zcore);
  memory->destroy(qmin);
  memory->destroy(qmax);
  memory->destroy(qf);
}

/* ---------------------------------------------------------------------- */

void FixQEqCtip::init()
{


  if (!atom->q_flag)
    error->all(FLERR,"Fix qeq/ctip requires atom attribute q");

  ngroup = group->count(igroup);
  if (ngroup == 0) error->all(FLERR,"Fix qeq/ctip group has no atoms");

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix  = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  extract_ctip();
  for (int i = 1; i <= atom->ntypes; i++) {
    if (zeta[i] == 0.0)
      error->all(FLERR,"Invalid parameters for fix qeq/ctip");
  }

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

}

/* ---------------------------------------------------------------------- */

void FixQEqCtip::extract_ctip()
{

  memory->create(chi,atom->ntypes+1,"qeq:chi");
  memory->create(eta,atom->ntypes+1,"qeq:eta");
  memory->create(gamma,atom->ntypes+1,"qeq:gamma");
  memory->create(zeta,atom->ntypes+1,"qeq:zeta");
  memory->create(zcore,atom->ntypes+1,"qeq:zcore");
  memory->create(qmin,atom->ntypes+1,"qeq:qmin");
  memory->create(qmax,atom->ntypes+1,"qeq:qmax");
  memory->create(qf,atom->natoms,"qeq:qf");

  pair = (Pair *) force->pair_match("coul/ctip",1);
  if (pair == NULL) error->all(FLERR,"No pair coul/ctip for fix qeq/ctip");
  int tmp;
  chi = (double *) pair->extract("chi",tmp);
  eta = (double *) pair->extract("eta",tmp);
  gamma = (double *) pair->extract("gamma",tmp);
  zeta = (double *) pair->extract("zeta",tmp);
  zcore = (double *) pair->extract("zcore",tmp);
  qmin = (double *) pair->extract("qmin",tmp);
  qmax = (double *) pair->extract("qmax",tmp);
  kspacetype = *(int *) pair->extract("kspacetype",tmp);
  ctip_omega = *(double *) pair->extract("omega",tmp);
  wolfewald = *(double *) pair->extract("wolfewald",tmp);
  cut_coul = *(double *) pair->extract("cut_coul",tmp);
  if (chi == NULL || eta == NULL || gamma == NULL
                  || zeta == NULL || zcore == NULL
                  || qmin == NULL || qmax == NULL
                  || &ctip_omega == NULL || &kspacetype == NULL
                  || &cut_coul == NULL || &wolfewald == NULL)
    error->all(FLERR,
	"Fix qeq/ctip could not extract params from pair coul/ctip");
}

/* ---------------------------------------------------------------------- */

void FixQEqCtip::conjugate_gradient()
{

  int inum = list->inum, *ilist = list->ilist;
  int *numneigh = list->numneigh, **firstneigh = list->firstneigh;
  int iloop, i, j, ii, jj;
  double g[inum], g_sum, g_sum_all; // g = -dUdq
  double d[inum], d_sum, d_sum_all;
  double p[inum], max_p, max_p_all, sigma_k;
  double Error = 1.0, g_sum0 = 1.0;
  for (ii = 0; ii < inum; ii++) p[ii] = 0;

  //printf("%d| inum = %d\n",comm->me,inum);

  for (iloop = 0; iloop < maxiter; iloop++) {

    // calculate g = -dUdq
    for (ii = 0; ii < inum; ii++) {
        g[ii] = 0.0;
        i = ilist[ii];
        if (atom->mask[i] & groupbit)
             g[ii] = -dUdq_self(atom->q[i],atom->type[i]);
    }
    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
        if (atom->mask[i] & groupbit) {
         for (jj = 0; jj < numneigh[i]; jj++) {
          j = firstneigh[i][jj];
          j &= NEIGHMASK;
          g[ii] -= dUdq_two_body(i,j,atom->q[j]);
         }
      }
    }

    g_sum = 0.0; g_sum_all = 0.0;
    for (ii = 0; ii < inum; ii++) g_sum += g[ii];
    MPI_Allreduce(&g_sum,&g_sum_all,1,MPI_DOUBLE,MPI_SUM,world);
    g_sum = g_sum_all/g_sum0; g_sum0 = g_sum_all;
//    if (comm->me == 0) printf("g_sum = %f\n",g_sum_all);

    // calculate d
    for (ii = 0; ii < inum; ii++) d[ii] = g[ii]+p[ii]*g_sum;
    d_sum = 0.0; d_sum_all = 0.0;
    for (ii = 0; ii < inum; ii++) d_sum += d[ii];
    MPI_Allreduce(&d_sum,&d_sum_all,1,MPI_DOUBLE,MPI_SUM,world);
    d_sum = d_sum_all/ngroup;
//    if (comm->me == 0) printf("d_sum = %f\n",d_sum);

    // calculate p, solve for sigma_k, update charges
    for (ii = 0; ii < inum; ii++) p[ii] = d[ii]-d_sum; // lead to charge neutral
    sigma_k = brent_solver(p); // solve sigma_k in qf = qi + sigma_k*p
    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      if (atom->mask[i] & groupbit)
          atom->q[i] += sigma_k*p[ii]; // update charges
    }
    pack_flag = 2;
    comm->forward_comm_fix(this);
//    if (comm->me == 0) printf("sigma_k = %f\n",sigma_k);

    // check if Error is within tolerance
    max_p = -1.0; max_p_all = -1.0;
    for (ii = 0; ii < inum; ii++) {
        p[ii] = fabs(p[ii]);
        if (p[ii] > max_p) max_p = p[ii];
    }
    MPI_Allreduce(&max_p,&max_p_all,1,MPI_DOUBLE,MPI_MAX,world);
    Error = max_p_all*fabs(sigma_k); // max(|p|)*|sigma_k|
//    if (comm->me == 0) printf("max_p_all = %f\n",max_p_all);

//    if (comm->me == 0) printf("error = %lg, tolerance = %lg\n",Error,tolerance);
    if (Error <= tolerance) return;
  }
  if (comm->me == 0) {
      char str[128];
      sprintf(str,"Charges did not converge at step " BIGINT_FORMAT
          ": maxiter = %d, CG error = %lg",update->ntimestep,maxiter,Error);
      error->warning(FLERR,str);
  }
  return;
}

/* ---------------------------------------------------------------------- */

double FixQEqCtip::dUdq_self(double qi, int itype)
{

/*
  X_i + qi*J_i  (electronegativity + self-Coulomb repulsion)
  + 2w(1-(qi-qmin)/|qi-qmin|)(qi-qmin)  correction for cation
  + 2w(1-(qmax-qi)/|qi-qmax|)(qi-qmax)  correction for anion
  correction terms keep charges in bound
*/

 double correction = 0.0;
 if (qi < qmin[itype]) 
    correction = 4.0*ctip_omega * (qi-qmin[itype]);
 else if (qi > qmax[itype]) 
    correction = 4.0*ctip_omega * (qi-qmax[itype]);

 double qqrd2e = force->qqrd2e;
 if (kspacetype == 1) { // wolf
   double woself = 0.50*erfc(wolfewald*cut_coul)/cut_coul + wolfewald/MY_PIS; // kc constant not yet multiplied.
   return chi[itype] + qi*(eta[itype]-2*qqrd2e*woself) + correction;
  }
 if (kspacetype == 2) // ewald
   return chi[itype] + qi*(eta[itype]-2*qqrd2e*wolfewald/MY_PIS) + correction;

  return 0.0;
}

/* ---------------------------------------------------------------------- */

double FixQEqCtip::dUdq_two_body(int i, int j, double qj)
{

 double delr[3], rsq;
 for (int k = 0; k <= 2; k++) delr[k] = atom->x[i][k] - atom->x[j][k];
 rsq = delr[0]*delr[0] + delr[1]*delr[1] + delr[2]*delr[2];
 if (rsq > cutoff_sq) return 0.0;
 
 int itype = atom->type[i], jtype = atom->type[j];
 if (kspacetype == 1) // wolf
   return dUdq_two_body_wolf(zeta[itype],zeta[jtype],
            zcore[jtype],sqrt(rsq),qj);

 if (kspacetype == 2) // ewald
   return dUdq_two_body_ewald(zeta[itype],zeta[jtype],
            zcore[jtype],sqrt(rsq),qj);

 return 0.0;
}

/* ---------------------------------------------------------------------- */

double FixQEqCtip::brent_solver(double *p)
{ 
 // solves for sigma_k in  qf = qi + sigma_k * p
 // sigma_k is the distance, p is the direction

 double a = -1.0, b = 1.0;
 double fa = calculate_qf(a,p);
 double fb = calculate_qf(b,p);

 double c = b, fc = fb;
 double xm,s,d,e,x,y,z;
 
 for (int iloop = 0; iloop < brent_maxiter; iloop++) {

     if (fb*fc > 0.0) {
        c = a; d = b-a; e = d; fc = fa;
     }
     if (fabs(fc) < fabs(fb)) {
        a = b; b = c; c = a; fa = fb; fb = fc; fc = fa;
     }
     xm = 0.5*(c-b);

     if (fabs(xm) < brent_critical || 
         fabs(fb) < brent_tolerance) {
            return b;
         }

     if (fabs(e) >= brent_critical && 
         fabs(fa) > fabs(fb)) {
         s = fb/fa;
         if (a == c) {
             x = 2.0*xm*s; y = 1.0-s;
         } else {
             y = fa/fc; z = fb/fc;
             x = s * (2.0*xm*y*(y-z) - (b-a)*(z-1.0));
             y = (y-1.0)*(z-1.0)*(s-1.0);
         }
         if (x > 0.0) y *= -1;
         x = fabs(x);
         if (2.0*x < fmin(3.0*xm*y - 
                fabs(brent_critical*y),fabs(e*y))) {
             e = d; d = x/y;
         } else {
             d = xm; e = d;
         }
     } else {
         d = xm; e = d;
     }
     a = b; fa = fb;
     if (fabs(d) > brent_critical) {
         b += d;
     } else {
         // sign of xm, value of brent_critical
         b += ((brent_critical<xm)-(xm<brent_critical)) 
              * brent_critical; // SIGN(brent_critical,xm)
     }
     fb = calculate_qf(b,p);
 }
 error->all(FLERR,"Brent solver limit reached");
}

/* ---------------------------------------------------------------------- */
double FixQEqCtip::calculate_qf(double c, double *p)
{

 // qf = qi + sigma_k * p  (sigma_k = c here)
 int inum = list->inum, *ilist = list->ilist;
 int *numneigh = list->numneigh, **firstneigh = list->firstneigh;
 int ii, jj, i, j;
 double g[inum];

 for (ii = 0; ii < inum; ii++) {
     i = ilist[ii];
     if (atom->mask[i] & groupbit)
         qf[i] = atom->q[i]+c*p[ii];
 }
 pack_flag = 1;
 comm->forward_comm_fix(this);
 for (ii = 0; ii < inum; ii++) {
     g[ii] = 0.0;
     i = ilist[ii];
     if (atom->mask[i] & groupbit)
         g[ii] = -dUdq_self(qf[ii],atom->type[i]);
 }
 for (ii = 0; ii < inum; ii++) {
   i = ilist[ii];
   if (atom->mask[i] & groupbit) {
     for (jj = 0; jj < numneigh[i]; jj++) {
       j = firstneigh[i][jj];
       j &= NEIGHMASK;
       g[ii] -= dUdq_two_body(i,j,qf[j]);
     }
   }
 }
 
 double fc = 0.0, fc_all = 0.0;
 for (ii = 0; ii < inum; ii++) fc -= g[ii]*p[ii];
 MPI_Allreduce(&fc,&fc_all,1,MPI_DOUBLE,MPI_SUM,world);
 return fc_all;

}

/* ---------------------------------------------------------------------- */

double FixQEqCtip::dUdq_two_body_ewald(double zei, double zej, double zj,
        double r, double qj)
{

  double exp2zir = exp(-2.0*zei*r);
  double zei2 = zei*zei;
  double zei4 = zei2*zei2;
  double zei6 = zei2*zei4;

  double exp2zjr = exp(-2.0*zej*r);
  double zej2 = zej*zej;
  double zej4 = zej2*zej2;
  double zej6 = zej2*zej4;

  double rinv = 1.0/r;
  double ci_jfi = -zei*exp2zir - rinv*exp2zir;

  double ci_fifj;
  if (zei == zej) {
    double sm1 = 11.0/8.0;// * zei;
    double sm2 = 3.00/4.0;// * zei2;
    double sm3 = 1.00/6.0;// * zei*zei2;
    ci_fifj = -exp2zir*(rinv + zei*(sm1 + sm2*zei*r + sm3*zei2*r*r));
  } else {
    double e1 = zei*zej4/((zei+zej)*(zei+zej)*(zei-zej)*(zei-zej));
    double e2 = zej*zei4/((zei+zej)*(zei+zej)*(zej-zei)*(zej-zei));
    double e3 = (3.0*zei2*zej4-zej6) /
         ((zei+zej)*(zei+zej)*(zei+zej)*(zei-zej)*(zei-zej)*(zei-zej));
    double e4 = (3.0*zej2*zei4-zei6) /
         ((zei+zej)*(zei+zej)*(zei+zej)*(zej-zei)*(zej-zei)*(zej-zei));
    ci_fifj = -exp2zir*(e1+e3/r) - exp2zjr*(e2+e4/r);
  }

  // k_c * (  Z([j|f_i]-[f_i|f_j]) + [f_i|f_j]  )
  return force->qqrd2e * (  zj*(ci_jfi - ci_fifj) + 
                            qj*(ci_fifj + erfc(wolfewald*r)*rinv)  );
}

/* ---------------------------------------------------------------------- */

double FixQEqCtip::dUdq_two_body_wolf(double zei, double zej, double zj,
        double r, double qj)
{

  double exp2zir = exp(-2.0*zei*r);
  double zei2 = zei*zei;
  double zei4 = zei2*zei2;
  double zei6 = zei2*zei4;

  double exp2zjr = exp(-2.0*zej*r);
  double zej2 = zej*zej;
  double zej4 = zej2*zej2;
  double zej6 = zej2*zej4;

  double rinv = 1.0/r;
  double rc = cutoff;
  double rcinv = 1.0/rc;
  double rcinv2 = rcinv*rcinv;
  double exp2zirsh = exp(-2.0*zei*rc);
  double exp2zjrsh = exp(-2.0*zej*rc);

  double eshift, fshift;

  eshift = -zei*exp2zirsh - rcinv*exp2zirsh;
  fshift = 2.0*zei2*exp2zirsh + rcinv2*exp2zirsh + 2.0*zei*rcinv*exp2zirsh;
  double ci_jfi = -zei*exp2zir - rinv*exp2zir - eshift - (r-rc)*fshift;

  double ci_fifj;
  if (zei == zej) {
    double sm1 = 11.0/8.0;// * zei;
    double sm2 = 3.00/4.0;// * zei2;
    double sm3 = 1.00/6.0;// * zei*zei2;
    eshift = -exp2zirsh*(rcinv + zei*(sm1 + sm2*zei*rc + sm3*zei2*rc*rc));
    fshift =  exp2zirsh*(rcinv2 + 2.0*zei*rcinv 
              + zei2*(2.0 + 7.0/6.0*zei*rc + 1.0/3.0*zei2*rc*rc));
    //fshift =  exp2zirsh*(rcinv2 + 2.0*zei*rcinv 
    //          + zei2*zei*(11.0/4.0 + 3.0/2.0*zei2*rc + 1.0/3.0*zei2*zei2*rc*rc - 3.0/4.0*zei - 1.0/3.0*zei2*zei*rc));
    ci_fifj = -exp2zir*(rinv + zei*(sm1 + sm2*zei*r + sm3*zei2*r*r))
	          - eshift - (r-rc)*fshift;
  } else {
    double e1 = zei*zej4/((zei+zej)*(zei+zej)*(zei-zej)*(zei-zej));
    double e2 = zej*zei4/((zei+zej)*(zei+zej)*(zej-zei)*(zej-zei));
    double e3 = (3.0*zei2*zej4-zej6) /
                ((zei+zej)*(zei+zej)*(zei+zej)*(zei-zej)*(zei-zej)*(zei-zej));
    double e4 = (3.0*zej2*zei4-zei6) /
                ((zei+zej)*(zei+zej)*(zei+zej)*(zej-zei)*(zej-zei)*(zej-zei));
    eshift = -exp2zirsh*(e1+e3/rc) - exp2zjrsh*(e2+e4/rc);
    fshift = (exp2zirsh*(2.0*zei*(e1+e3/rc) + e3*rcinv2)
	         + exp2zjrsh*(2.0*zej*(e2+e4/rc) + e4*rcinv2));
    ci_fifj = -exp2zir*(e1+e3/r) - exp2zjr*(e2+e4/r)
	          - eshift - (r-rc)*fshift;
  }

  // k_c * (  Z([j|f_i]-[f_i|f_j]) + [f_i|f_j]  )
  return force->qqrd2e * (  zj*(ci_jfi - ci_fifj) + 
                            qj*(erfc(wolfewald*r)/r - erfc(wolfewald*rc)/rc + ci_fifj) );
}

/* ---------------------------------------------------------------------- */

double FixQEqCtip::memory_usage()
{
  double bytes = atom->natoms * sizeof(double);
  return bytes;
}

/* ---------------------------------------------------------------------- */

int FixQEqCtip::pack_forward_comm(int n, int *list, double *buf,
                                  int pbc_flag, int *pbc)
{
  int m;

  if (pack_flag == 1)
    for (m = 0; m < n; m++) buf[m] = qf[list[m]];
  else if (pack_flag == 2)
    for (m = 0; m < n; m++) buf[m] = atom->q[list[m]];

  return m;
}

/* ---------------------------------------------------------------------- */

void FixQEqCtip::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m;

  if (pack_flag == 1)
    for (m = 0, i = first; m < n; m++, i++) qf[i] = buf[m];
  else if (pack_flag == 2)
    for (m = 0, i = first; m < n; m++, i++) atom->q[i] = buf[m];
}

/* ---------------------------------------------------------------------- */

void FixQEqCtip::pre_force(int vflag)
{
  if (update->ntimestep % nevery) return;
  conjugate_gradient();
  if (force->kspace) force->kspace->qsum_qsq();
}

