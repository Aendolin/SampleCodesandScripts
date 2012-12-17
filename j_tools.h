/*
This source code is for the j_hydration program. 

These functions assume that SPCE water is used, and that HW1 and HW2 come directly after OW.
Only the OW then need be specified in the index file.

A vector type of '0' corresponds to the OH bond, '1' corresponds to the water dipole. 
*/

#define LARGE 100000

enum { OH, DIPOLE };
enum { SMOOTH, ROUGH };

typedef struct
{
  real Xmin;
  real Xmax;
  real Ymin;
  real Ymax;
} Limits;

void define_limits(Limits *limits, int k, int *isize, rvec *x, atom_id **index,t_pbc *pbc, matrix box);
real within_limits(Limits limits, int i, rvec *x, atom_id **index,t_pbc *pbc, matrix box);
void print_limits(Limits limits);
real calculate_op(int i,rvec *x,atom_id **index, t_pbc *pbc, matrix box, int vectortype, int orderparameter, int zdistancetype);
real calculate_z_distance(int t, rvec *x, atom_id **index, int *isize, matrix box, t_pbc *pbc, int zdistancetype, int *closest, int *valid);
real find_midplane(int k,int *isize, rvec *x, atom_id **index, t_pbc *pbc, matrix box);
