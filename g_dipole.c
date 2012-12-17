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

int main(int argc, char *argv[])
{
  static char *desc[] = {
    "g_dipole calculates the average dipole of a water molecule witht the z-axis.",
    "SPCE water is assumed, as well as OW,H,H ordering in the gro/trajectory files.",
    "Choose 'smooth' or 'rough' distance determination." 
    "This program also assumes the hydrophobic core is in the middle of the box;",
    "you must shift your system if it is not!"
  };
  
  t_topology *top=NULL;
  rvec *x=NULL, dx;
  matrix box;
  int ePBC,i,j,k,n,natoms,bin,valid;
  t_pbc *pbc;
  int teller=0;
  real *time=NULL, d;
  int status;
  real t, min, z_distance;
  atom_id **index;
  int *isize, closest, slice; 
  char **grpname;
  FILE *fp;
  real *OP, *OP_count, *OP_sq, op_stderr,normalizer=0;
  real slab_height=0, avg_slab_height=0, slab_volume=0, avg_slab_volume=0, op=0;
  double oxygen_charge = -0.8476, hydrogen_charge = 0.4238;
  double water_mass=18.01528;
  static char *string1="op"; 
  static int nslices=100;
  static int mods=1;
  static int zdistancetype=SMOOTH;
  static int orderparameter=2;
  static int vectortype=OH;
  
  static t_pargs pa[] = {
    { "-nslices", FALSE, etINT, {&nslices},
      "Divide the box z-wise into this many slices" },
    { "-distance", FALSE, etINT, {&zdistancetype},
      "choose between 'smooth' or 'rough' distance" },
    { "-out", FALSE, etSTR, {&string1},
      "Output file" },
    { "-mods", FALSE, etINT, {&mods},
      "How many frames until you calculate nearest neighbors again"},
    { "-op", FALSE, etINT, {&orderparameter},
      "The order parameter (1 or 2)"},
    { "-vector", FALSE, etINT, {&vectortype},
      "The vector type OH(0) or DIPOLE(1)"}
  };
  
#define NPA asize(pa)
  
  t_filenm fnm[] = {
    { efTRX, "-f", NULL, ffREAD },
    { efTPX, NULL, NULL, ffREAD },
    { efNDX, NULL, NULL, ffOPTRD },
  };
  
#define NFILE asize(fnm)
  
  CopyRight(stderr,argv[0]);
  
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_BE_NICE,
		    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL);
  
  top=read_top(ftp2fn(efTPX,NFILE,fnm),&ePBC);
  snew(pbc,1);     
  
  snew(OP,nslices);
  snew(OP_sq,nslices);
  snew(OP_count,nslices);
  
  natoms=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
  
  if (zdistancetype == SMOOTH) 
    {
      snew(index,1);
      snew(isize,1);
      snew(grpname,1); 
      printf("Please enter an index file with the oxygens of the water\n");
      get_index(&top->atoms,ftp2fn_null(efNDX,NFILE,fnm),1,isize,index,grpname);
    }
  else
    {
      snew(index,2);
      snew(isize,2);
      snew(grpname,2); 
      printf("Doing rough surface calculations\n");
      printf("Please enter an index file with the oxygens of the water, and then a surface atom index\n");
      get_index(&top->atoms,ftp2fn_null(efNDX,NFILE,fnm),2,isize,index,grpname);
    }
  
  do {
    set_pbc(pbc,ePBC,box);
	
	if (zdistancetype == SMOOTH)
	  slab_height = box[2][2] / nslices;	
	else
	  slab_height = 5.0 / nslices;
	
	avg_slab_height += slab_height;
	avg_slab_volume += box[0][0]*box[1][1]*slab_height;
	
	for (i=0;i<isize[0];i++) // cycle through all water oxygens
	  {
	    
	    op = calculate_op(i,x,index,pbc,box,vectortype,orderparameter,zdistancetype);
	    z_distance = calculate_z_distance(i,x,index,isize,box,pbc,zdistancetype,&closest,&valid);
	    
	    if (valid)
	      {
		if (zdistancetype == SMOOTH)
		  slice = (int)(z_distance/slab_height);
		else
		  slice = (int)(z_distance/slab_height+nslices/2.0);
		
		if (zdistancetype == ROUGH && x[closest][2] < box[2][2]/2 && orderparameter == 1)
		  op *= -1;
		
		OP[slice]+=op;
		OP_sq[slice]+=op*op;
		OP_count[slice]+=1;
	      }
	  }
	teller++;
	
  } while (read_next_x(status,&t,natoms,x,box));
  
  fp = fopen(string1,"w");
  
  avg_slab_height /= teller;
  //avg_slab_volume /= teller;
  
  //normalizer = water_mass / (teller*100.0*6.022*avg_slab_volume); // 100
  
  for (i=0;i<nslices;i++)
    {
      if (OP_count[i] != 0)
	{
	  OP[i] /= OP_count[i];
	  OP_sq[i] /= OP_count[i];
	}
    }  
  
  for (i=0;i<nslices;i++)
    {
      if(OP_count[i] != 0)
	op_stderr = sqrt(OP_sq[i]-pow(OP[i],2))/sqrt(OP_count[i]);
      else
	op_stderr = 0;
      
      if (zdistancetype == SMOOTH)
	fprintf(fp, "%f %f %f\n",avg_slab_height*i,OP[i],op_stderr);
      else
        fprintf(fp, "%f %f %f\n",avg_slab_height*i-2.5,OP[i],op_stderr);
      
    }
  
  fclose(fp);
  close_trj(status);
  thanx(stderr);
  return 0;
}
