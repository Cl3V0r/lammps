// clang-format off
 /* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_wall_reflect_cylinder.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "lattice.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cstring>
#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixWallReflectCylinder::FixWallReflectCylinder(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  nwall(0)
{
  if (narg < 4) error->all(FLERR,"Illegal fix wall/reflect command");

  // let child class process all args

  //if (strcmp(arg[2],"wall/reflect/stochastic") == 0) return;

  dynamic_group_allow = 1;

  // parse args

  nwall = 0;
  int scaleflag = 1;

  int iarg = 3;
  R = utils::numeric(FLERR,arg[iarg],false,lmp);
  yc = utils::numeric(FLERR,arg[iarg+1],false,lmp);
  zc = utils::numeric(FLERR,arg[iarg+2],false,lmp);
  //std::cout << "R: " << R << " yc: " << yc << " zc: " << zc << std::endl;

  //while (iarg < narg) {
    //if ((strcmp(arg[iarg],"xlo") == 0) || (strcmp(arg[iarg],"xhi") == 0) ||
    //    (strcmp(arg[iarg],"ylo") == 0) || (strcmp(arg[iarg],"yhi") == 0) ||
    //    (strcmp(arg[iarg],"zlo") == 0) || (strcmp(arg[iarg],"zhi") == 0)) {
    //  if (iarg+2 > narg) error->all(FLERR,"Illegal fix wall/reflect command");
//
    //  int newwall;
    //  if (strcmp(arg[iarg],"xlo") == 0) newwall = XLO;
    //  else if (strcmp(arg[iarg],"xhi") == 0) newwall = XHI;
    //  else if (strcmp(arg[iarg],"ylo") == 0) newwall = YLO;
    //  else if (strcmp(arg[iarg],"yhi") == 0) newwall = YHI;
    //  else if (strcmp(arg[iarg],"zlo") == 0) newwall = ZLO;
    //  else if (strcmp(arg[iarg],"zhi") == 0) newwall = ZHI;
//
    //  for (int m = 0; (m < nwall) && (m < 6); m++)
    //    if (newwall == wallwhich[m])
    //      error->all(FLERR,"Wall defined twice in fix wall/reflect command");
//
    //  wallwhich[nwall] = newwall;
    //  if (strcmp(arg[iarg+1],"EDGE") == 0) {
    //    wallstyle[nwall] = EDGE;
    //    int dim = wallwhich[nwall] / 2;
    //    int side = wallwhich[nwall] % 2;
    //    if (side == 0) coord0[nwall] = domain->boxlo[dim];
    //    else coord0[nwall] = domain->boxhi[dim];
    //  } else if (utils::strmatch(arg[iarg+1],"^v_")) {
    //    wallstyle[nwall] = VARIABLE;
    //    varstr[nwall] = utils::strdup(arg[iarg+1]+2);
    //  } else {
    //    wallstyle[nwall] = CONSTANT;
    //    coord0[nwall] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    //  }
//
    //  nwall++;
    //  iarg += 2;
//
    //} else if (strcmp(arg[iarg],"units") == 0) {
    //  if (iarg+2 > narg) error->all(FLERR,"Illegal wall/reflect command");
    //  if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
    //  else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
    //  else error->all(FLERR,"Illegal fix wall/reflect command");
    //  iarg += 2;
//
    //} else error->all(FLERR,"Illegal fix wall/reflect command");
  //}

  // error check

  //if (nwall == 0) error->all(FLERR,"Illegal fix wall command");

  //for (int m = 0; m < nwall; m++) {
  //  if ((wallwhich[m] == XLO || wallwhich[m] == XHI) && domain->xperiodic)
  //    error->all(FLERR,"Cannot use fix wall/reflect in periodic dimension");
  //  if ((wallwhich[m] == YLO || wallwhich[m] == YHI) && domain->yperiodic)
  //    error->all(FLERR,"Cannot use fix wall/reflect in periodic dimension");
  //  if ((wallwhich[m] == ZLO || wallwhich[m] == ZHI) && domain->zperiodic)
  //    error->all(FLERR,"Cannot use fix wall/reflect in periodic dimension");
  //}
  if ((domain->zperiodic) || (domain->yperiodic)){
    error->all(FLERR,"Cannot use fix wall/reflect in periodic dimension");
  } 

  if(domain->dimension == 2){
    error->all(FLERR,"Cannot use fix wall/reflect in 2D simulation");
  }

  //for (int m = 0; m < nwall; m++)
  //  if ((wallwhich[m] == ZLO || wallwhich[m] == ZHI) && domain->dimension == 2)
  //    error->all(FLERR,
  //               "Cannot use fix wall/reflect zlo/zhi for a 2d simulation");
//
  // scale factors for CONSTANT and VARIABLE walls

  //int flag = 0;
  //for (int m = 0; m < nwall; m++)
  //  if (wallstyle[m] != EDGE) flag = 1;
//
  //if (flag) {
  //  if (scaleflag) {
  //    xscale = domain->lattice->xlattice;
  //    yscale = domain->lattice->ylattice;
  //    zscale = domain->lattice->zlattice;
  //  }
  //  else xscale = yscale = zscale = 1.0;
//
  //  for (int m = 0; m < nwall; m++) {
  //    if (wallstyle[m] != CONSTANT) continue;
  //    if (wallwhich[m] < YLO) coord0[m] *= xscale;
  //    else if (wallwhich[m] < ZLO) coord0[m] *= yscale;
  //    else coord0[m] *= zscale;
  //  }
  //}
//
  //// set varflag if any wall positions are variable
//
  //varflag = 0;
  //for (int m = 0; m < nwall; m++)
  //  if (wallstyle[m] == VARIABLE) varflag = 1;
}

/* ---------------------------------------------------------------------- */

FixWallReflectCylinder::~FixWallReflectCylinder()
{
  if (copymode) return;

  //for (int m = 0; m < nwall; m++)
  //  if (wallstyle[m] == VARIABLE) delete [] varstr[m];
}

/* ---------------------------------------------------------------------- */

int FixWallReflectCylinder::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallReflectCylinder::init()
{
  //for (int m = 0; m < nwall; m++) {
  //  if (wallstyle[m] != VARIABLE) continue;
  //  varindex[m] = input->variable->find(varstr[m]);
  //  if (varindex[m] < 0)
  //    error->all(FLERR,"Variable name for fix wall/reflect does not exist");
  //  if (!input->variable->equalstyle(varindex[m]))
  //    error->all(FLERR,"Variable for fix wall/reflect is invalid style");
  //}

  int nrigid = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->rigid_flag) nrigid++;

  if (nrigid && comm->me == 0)
    error->warning(FLERR,"Should not allow rigid bodies to bounce off "
                   "relecting walls");
}

/* ---------------------------------------------------------------------- */

void FixWallReflectCylinder::post_integrate()
{
  //double coord;

  // coord = current position of wall
  // evaluate variable if necessary, wrap with clear/add

  if (varflag) modify->clearstep_compute();

  //for (int m = 0; m < nwall; m++) {
  //  if (wallstyle[m] == VARIABLE) {
  //    coord = input->variable->compute_equal(varindex[m]);
  //    if (wallwhich[m] < YLO) coord *= xscale;
  //    else if (wallwhich[m] < ZLO) coord *= yscale;
  //    else coord *= zscale;
  //  } else coord = coord0[m];
//
  //}
  wall_particle(R,yc,zc);

  if (varflag) modify->addstep_compute(update->ntimestep + 1);
}

/* ----------------------------------------------------------------------
   this method may be overwritten by a child class
------------------------------------------------------------------------- */

void FixWallReflectCylinder::wall_particle(double R, double yc, double zc)
{
  
  int i, dim; //,dim,side;

  double **x = atom->x;
  double **v = atom->v;
  double **vest = atom->vest; //vest is needed for DPD 
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  double R_2 = R*R;
  double dy, dz, w;

  for (i = 0; i < nlocal; i++) { //nlocal is the number of atoms of this processor
    //setmask() sets on bits POST_INTEGRATE and POST_INTEGRATE_RESPA
    //checks if mask and groupbit are the same -> if they are, then the atom is in the group specified by the fix
    if (mask[i] & groupbit) { 
        dy = x[i][1] - yc; 
        dz = x[i][2] - zc;
        if (dy*dy + dz*dz > R_2) {
          w = 2*(1-R/sqrt(dy*dy + dz*dz)); //calculate projection
          x[i][1] -= w*dy;
          x[i][2] -= w*dz;

          for (dim=0; dim<3; dim++) {
              v[i][dim] = -v[i][dim];
              vest[i][dim] = -vest[i][dim];
          }
        }
    }
  }
}

