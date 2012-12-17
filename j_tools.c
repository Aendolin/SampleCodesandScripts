#include <typedefs.h>
#include "smalloc.h"
#include "macros.h"
#include "math.h"
#include "xvgr.h"
#include "copyrite.h"
#include "statutil.h"
#include "string2.h"
#include "vec.h"
#include "rdgroup.h"
#include "pbc.h"
#include "gmx_fatal.h"
#include "futil.h"
#include "gstat.h"
#include "pbc.h"
#include "j_tools.h"

real find_midplane(int k, int *isize, rvec *x, atom_id **index, t_pbc *pbc, matrix box)
{
  /* input: surface atom group */
  real upperavg=0, loweravg=0, midavg=0;
  int count=0,i;
  
  for (i=0;i<isize[k];i++)
    {
      midavg+=x[index[k][i]][2];
      count++;
    }
  
  midavg /= count;
  
  return midavg;
  
}

void define_limits(Limits *limits, int k, int *isize, rvec *x, atom_id **index,t_pbc *pbc, matrix box)
{
  int i;
  
  limits->Xmin=LARGE;
  limits->Xmax=0;
  limits->Ymin=LARGE;
  limits->Ymax=0;
  
  for(i=0;i<isize[k];i++)
    {
      if (x[index[k][i]][0] > limits->Xmax)
	limits->Xmax = x[index[k][i]][0];
      if (x[index[k][i]][0] < limits->Xmin)
	limits->Xmin = x[index[k][i]][0];
      if (x[index[k][i]][1] > limits->Ymax)
	limits->Ymax = x[index[k][i]][1];
      if (x[index[k][i]][1] < limits->Ymin)
	limits->Ymin = x[index[k][i]][1];
    }

  limits->Xmax -= 0.1;
  limits->Xmin += 0.1;
  limits->Ymax -= 0.1;
  limits->Xmin += 0.1;

  return;
}

real within_limits(Limits limits, int i, rvec *x, atom_id **index,t_pbc *pbc, matrix box)
{
  int valid=0;
  //printf("valid is %d\n",valid);
  if (x[index[0][i]][0] < limits.Xmax && x[index[0][i]][0] > limits.Xmin && x[index[0][i]][1] < limits.Ymax && x[index[0][i]][1] > limits.Ymin)
    valid = 1;
  //printf("valid is now %d\n",valid);
  return valid;
  
}

void print_limits(Limits limits)
{
  printf("Xmin is %f\n",limits.Xmin);
  printf("Xmax is %f\n",limits.Xmax);
  printf("Ymin is %f\n",limits.Ymin);
  printf("Ymax is %f\n",limits.Ymax);
}

real calculate_op(int i,rvec *x,atom_id **index,t_pbc *pbc, matrix box, int vectortype, int orderparameter, int zdistancetype)
{
  real op, dipmag, cangle;
  rvec vec1, vec2, vec3;

  if (vectortype == OH)
      pbc_dx(pbc,x[index[0][i]+1],x[index[0][i]],vec3);
  if (vectortype == DIPOLE)      
    {
      pbc_dx(pbc,x[index[0][i]+1],x[index[0][i]],vec1);
      pbc_dx(pbc,x[index[0][i]+2],x[index[0][i]],vec2);
      rvec_add(vec1,vec2,vec3);
    }
  
  cangle = vec3[2]/norm(vec3);
  
  if (orderparameter == 1)
    op = cangle;
  if (orderparameter == 2)
    op = 0.5*(3*cangle*cangle - 1);
  
  return op;
  
}

real calculate_z_distance(int i, rvec *x, atom_id **index, int *isize, matrix box, t_pbc *pbc, int zdistancetype, int *closest, int *valid)
{
  int j;
  real min = 1000.0, d = 0, z_distance=-1000;
  rvec dx;
  
  *valid = 0;
  *closest = -1;
  
  if (zdistancetype == SMOOTH)
    {
      z_distance = x[index[0][i]][2];
      *valid = 1;
    }
  if (zdistancetype == ROUGH)
    {
      for(j=0;j<isize[1];j++)
	{
	  pbc_dx(pbc,x[index[0][i]],x[index[1][j]],dx);
	  d = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
	  
	  if (d < min && fabs(x[index[0][i]][2]-x[index[1][j]][2]) < 2.5) 
	    {
	      *closest = index[1][j];
	      min = d;
	    }
	}
      
      if (*closest != -1)
      //if (fabs(x[index[0][i]][2]-x[*closest][2]) < 2.0) // make sure water is at most 2 nm away from surface atom
	{
	  //printf("closest is %d\n",*closest);
	  //printf("closest z: %f\n",x[*closest][2]);
	  //printf("water z: %f\n", x[index[0][i]][2]);
	  
	  *valid = 1;
	  if (x[*closest][2] > box[2][2]/2)
	    z_distance = x[index[0][i]][2]-x[*closest][2];
	  else 
	    //op1 *= -1; make sure to add this to main code!!
	    z_distance = x[*closest][2]-x[index[0][i]][2];
	}
    } // end of ROUGH
  
  //printf("valid is %d\n",*valid);
  
  return z_distance;
}
