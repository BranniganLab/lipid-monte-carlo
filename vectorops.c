
#include "structures.h"
#include <math.h>
#include "vectorops.h"
//inline void vsub_periodic(vector vr, vector v1, vector v2, vector box_length);
inline void testinline()
{
	while (1 == 2)
 	{break;}
}


//inline void vsub_periodic(vector vr, vector v1, vector v2, vector box_length)
//{
//  int i;
//  extern int pboundary;
//      vr[X] = v1[X] - v2[X];
//      vr[Y] = v1[Y] - v2[Y];
//      vr[Z] = v1[Z] - v2[Z];
//     if (pboundary != NONE)
//	{if (vr[X] > box_length[X]/2)
//	  vr[X] = - (box_length[X] - vr[X]);
//	else if (vr[X] < -box_length[X]/2)
//	  vr[X] =  (box_length[X] + vr[X]);
//
//
//    	if (vr[Y] > box_length[Y]/2)
//	  vr[Y] = - (box_length[Y] - vr[Y]);
//	else if (vr[Y] < -box_length[Y]/2)
//	  vr[Y] =  (box_length[Y] + vr[Y]);
//
//
//    if (pboundary != PLANE)
//    	{if (vr[Z] > box_length[Z]/2)
//	  vr[Z] = - (box_length[Z] - vr[Z]);
//	else if (vr[Z] < -box_length[Z]/2)
//	  vr[Z] =  (box_length[Z] + vr[Z]);}
//   }
//
//}
/********************************************************/
double dsign(double a, double b)
{
  double answer;
  if (b > 0)
    answer = fabs(a);
  else
    answer = -1 * fabs(a);
  return answer; 		
}
/********************************************************/
void vprint(vector v)
{
  printf("x: %f   y: %f   z: %f  \n", v[X], v[Y], v[Z]);
}
/********************************************************/
void vmult(vector vr, vector v1, vector v2)
{
  vr[X] = v1[X] * v2[X];
  vr[Y] = v1[Y] * v2[Y];
  vr[Z] = v1[Z] * v2[Z];
}
/************************************************************/
void vclear(vector vect)
{
  vect[X] = 0;
  vect[Y] = 0;
  vect[Z] = 0;
}
/******************************************************************/
void rotate_vect( vector vect, vector theta)
/*rotates vector by angle in EITHER x, y, or z direction*/ 
{

  vector r_matrix[DIMENSION];
  vector v;
  vector transform;
  vclear(transform);
  vclear(r_matrix[X]);
  vclear(r_matrix[Y]);
  vclear(r_matrix[Z]);
		
  if (theta[Z] != 0)
    {(r_matrix[X])[X] = cos(theta[Z]);
    (r_matrix[Y])[Y] = cos(theta[Z]);
    (r_matrix[X])[Y] = sin(theta[Z]);
    (r_matrix[Y])[X] = -sin(theta[Z]);
    (r_matrix[Z])[Z] = 1;

    }
		 	
  else  if (theta[X] != 0)
    {(r_matrix[X])[X] = 1;
    (r_matrix[Y])[Y] = cos(theta[X]);
    (r_matrix[Y])[Z] = sin(theta[X]);
    (r_matrix[Z])[Y] = -sin(theta[X]);
    (r_matrix[Z])[Z] = cos(theta[X]);
    }
		
  else if (theta[Y] != 0)
    {(r_matrix[X])[X] = cos(theta[Y]);
    (r_matrix[Y])[Y] = 1;
    (r_matrix[X])[Z] = sin(theta[Y]);
    (r_matrix[Z])[X] = -sin(theta[Y]);
    (r_matrix[Z])[Z] = cos(theta[Y]);
    }
  else
    {(r_matrix[X])[X] = 1;
    (r_matrix[Y])[Y] = 1;
    (r_matrix[Z])[Z] = 1;
    }
					
  v[X] = vdot(r_matrix[X], vect);
  v[Y] = vdot(r_matrix[Y], vect);
  v[Z] = vdot(r_matrix[Z], vect);
  if (fabs((vsquare(vect) -vsquare(v))) > 0.00001)
    {printf("discrepancy:  %f %f \n", vsquare(v) ,    vsquare(vect)    );
    printf("theta: \n");
    vprint(theta);
    printf("Matrix: \n");
    vprint(r_matrix[X]);
    vprint(r_matrix[Y]);
    vprint(r_matrix[Z]);
    }
  vcopy(vect, v);
}
		 	

/*********************************************************************/
double vdistance(vector v1, vector v2)
{
  //returns the square of the difference vector
  double radius;
  vector temp;
  vsub(temp, v1, v2);
  radius = vsquare(temp);
  return radius;
}

/**********************************************************************/	
double vdistance_periodic(vector v1, vector v2, vector box_length)
{
  //returns the square of the difference vector, including boundary conditions
  double radius;
  vector vr;


      vr[X] = v1[X] - v2[X];
      vr[Y] = v1[Y] - v2[Y];
      vr[Z] = v1[Z] - v2[Z];

	if (vr[X] > box_length[X]/2)
	  vr[X] = - (box_length[X] - vr[X]);
	else if (vr[X] < -box_length[X]/2)
	  vr[X] =  (box_length[X] + vr[X]);

    	if (vr[Y] > box_length[Y]/2)
	  vr[Y] = - (box_length[Y] - vr[Y]);
	else if (vr[Y] < -box_length[Y]/2)
	  vr[Y] =  (box_length[Y] + vr[Y]);


    	if (vr[Z] > box_length[Z]/2)
	  vr[Z] = - (box_length[Z] - vr[Z]);
	else if (vr[Z] < -box_length[Z]/2)
	  vr[Z] =  (box_length[Z] + vr[Z]);


  radius = vsquare(vr);
  return radius;
}

/*********************************************************************/
void vcopy(vector hold, vector copy)
{
  //copies one vector to another
  hold[X] = copy[X];
  hold[Y] = copy[Y];
  hold[Z] = copy[Z];
}

void vadd(vector vr, vector v1, vector v2)
{
  vr[X] = v1[X] + v2[X];
  vr[Y] = v1[Y] + v2[Y];
  vr[Z] = v1[Z] + v2[Z];
}

void vsub(vector vr, vector v1, vector v2)
{
  vr[X] = v1[X] - v2[X];
  vr[Y] = v1[Y] - v2[Y];
  vr[Z] = v1[Z] - v2[Z];
}

void vmult_scal(vector vr,vector v1,double coeff)
{
  vr[X] = v1[X]*coeff;
  vr[Y] = v1[Y]*coeff;
  vr[Z] = v1[Z] *coeff;
}

//#define psub(a,b, l)



double vdot(vector v1, vector v2)
{double total = 0;
 total = v1[X]*v2[X] +v1[Y]*v2[Y]+v1[Z]*v2[Z];
 return total;
}

double vsquare(vector v1)
{double total = 0;
 total = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
 return total;
}

void vunit(vector vr, vector v1)
{
	double coeff;
	coeff = 1/ sqrt(vsquare(v1));
	vr[X] = v1[X]*coeff;
  	vr[Y] = v1[Y]*coeff;
  	vr[Z] = v1[Z] *coeff;
}


void vresize(vector vr, vector v1, double length)
{
  vmult_scal(vr, v1, length/sqrt(vsquare(v1)));
}	



double myacos(vector v1, vector v2)
{
	vector rotate, a, b;
 	double theta;

	rotate[X] = 1 - v1[X];
 	rotate[Y] = -v1[Y];
  	rotate[Z] =  -v1[Z];
 	vadd(a, v1, rotate);
  	vadd(b, v2, rotate);
   	vunit(a,a); vunit(b,b);
   	theta = acos(vdot(a,b));
    	if (b[Y] < 0)
     	theta = PI - theta;
     return theta;
}

/**********************************/
void vcross(vector vr, vector v1, vector v2)
{
    vr[X] = v1[Y]*v2[Z] - v1[Z]*v2[Y];
    vr[Y] = -v1[X]*v2[Z] + v1[Z] * v2[X];
    vr[Z] = v1[X]*v2[Y] - v1[Y]*v2[X];
}

/*********************************/
void vrotate(vector vr, vector r, vector n, double phi)
{
    vector temp;
    double start, end;
    start = vsquare(r);
    vmult_scal(vr, r, cos(phi));
    vmult_scal(temp, n, vdot(n, r) *(1 - cos(phi)));
    vadd(vr, vr, temp);
    vcross(temp, r, n);
    vmult_scal(temp, temp, sin(phi));
    vadd(vr, vr, temp);
    end = vsquare(vr); 
}

