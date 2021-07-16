#include <stdlib.h>
#include <stdio.h>
#include "structures.h"
#include "constraints.h"
#include "vectorops.h"
#include <math.h>

int enforce_constraints(particle *the_membrane, vector box_length)
{
	int i; 
	chainindex ci; 
	extern int pnumber_user_constraints; 
	extern uconstraint *puser_constraints; 
	vector translation, new_point,r; 
	double radius; 
	
	for (i = 0; i < pnumber_user_constraints; i++)
		{
		ci.m = puser_constraints[i].ci.m; 
		ci.b = puser_constraints[i].ci.b; 
		if (puser_constraints[i].type == FIXED_POINT)
			move_molecule_to_point(the_membrane, ci, puser_constraints[i].point, box_length); 
		else if (puser_constraints[i].type == FIXED_ANGLE)
			{
				vsub_periodic(r, the_membrane[ci.m].chain[ci.b].position, puser_constraints[i].point,box_length);
				radius = sqrt(r[X]*r[X]+ r[Y]*r[Y]); 
				vmult_scal(new_point,puser_constraints[i].arrow,radius);
				vadd(new_point, new_point, puser_constraints[i].point); 
				new_point[Z] = the_membrane[ci.m].chain[ci.b].position[Z]; 				
				move_molecule_to_point(the_membrane,ci,new_point, box_length); 
			}
		}			
}

int move_molecule_to_point(particle *the_membrane,chainindex ci, vector new_point, vector box_length)
{
	vector translation; 
	int i; 
		vsub(translation, new_point, the_membrane[ci.m].chain[ci.b].position);
		translation[Z] = 0; 	
		for (i = 0; i < the_membrane[ci.m].chain_length; i++)
			vadd(the_membrane[ci.m].chain[i].position, the_membrane[ci.m].chain[i].position,translation);  
		
		}

/*************************************************/

int print_constraints(particle *the_membrane, mparticle *the_types, vector box_length)
{
	extern int pnumber_user_constraints;
	extern uconstraint *puser_constraints; 
	int i; 
	vector position; 
	FILE *fp; 
	int ibindex; 
	fp = fopen("mc.constraints","a"); 
	for (i = 0; i < pnumber_user_constraints; i++)
	{
		ibindex = the_types[the_membrane[puser_constraints[i].ci.m].model_index].ibindex; 
		if (puser_constraints[i].type == FIXED_POINT)
			vcopy(position,the_membrane[puser_constraints[i].ci.m].chain[ibindex].position); 
		else
			vcopy(position,the_membrane[puser_constraints[i].ci.m].chain[puser_constraints[i].ci.b].position); 
		fprintf(fp,"%lf %lf %lf\n", position[X], position[Y], position[Z]); 
		}
	fclose(fp); 
	
}
