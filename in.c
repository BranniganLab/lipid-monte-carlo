#include <stdlib.h>
#include <stdio.h>
#include "structures.h"
#include "in.h"
#include "vectorops.h"
#include "ghosts.h"
#include "math.h"
/******************************************************************/
int get_one_ensemble(FILE *fp, particle *the_membrane, mparticle *the_types)
{/*reads in one configuration of the membrane, from either a *.trj or *.in file*/
extern long pnum_particles;

extern double ptorsion;
extern double plength, pdiameter;
vector temp;
double a2,b2,c2;


double x, y, z, a, b, c,r,l;
int i,j;

//printf("reading xyz file in 5.0 format: beads are read in order"); 
if (fp == NULL)
{printf("bad file pointer!\n");
    return 1;}
if (read_model(the_membrane) == NO)
{
	for (i = 0; i < pnum_particles; i++)
	{
		the_membrane[i].model_index = 0; 
		the_membrane[i].domain_index = 0;
		the_membrane[i].chain_length = plength+2; 
	}
	the_types[0].cbend = ptorsion;
	
	the_types[0].chain_length = plength+2;
 	
	the_types[0].chain[0].type = HEAD; 
	the_types[0].chain[1].type = INTERFACE; 
	for (j = 2; j < the_types[0].chain_length; j++)
		the_types[0].chain[j].type = TAIL; 
	the_types[0].chain1start = 2; 
	the_types[0].chain2start = the_membrane[0].chain_length + 1; 
	the_types[0].ibindex = 1; 
}


for (i =0; i< pnum_particles; i++)
{
    if (feof(fp))
	  {
        if (i > 0)
            printf("WARNING: %ld particles specified but only %d read, may be bad file format. \n", pnum_particles, i);
        return 1;}
    r = 0.0;
    
    
    the_membrane[i].chain_length = the_types[the_membrane[i].model_index].chain_length; 
    fscanf(fp,"%lf %lf %lf %lf %lf %lf", &x, &y, &z, &a, &b, &c);
    the_membrane[i].chain[1].position[X]= x;
    the_membrane[i].chain[1].position[Y]= y;
    the_membrane[i].chain[1].position[Z] = z;
    //the_membrane[i].chain_length = plength+2;
    the_membrane[i].chain[1].ci.b = 1;
	the_membrane[i].chain[1].ci.m = i;
 	the_membrane[i].chain[1].user = NO;
	the_membrane[i].chain[1].typeptr = &(the_types[the_membrane[i].model_index].chain[1]); 

	//the_membrane[i].chain[1].type = INTERFACE; 
    /*if the particle is out of bounds, move it in bounds*/
    //check_bounds(the_membrane[i].crossings, the_membrane[i], pbox_length);
    //vsub(the_membrane[i].position,the_membrane[i].position,the_membrane[i].crossings);
    the_membrane[i].chain[0].position[X]= a;
    the_membrane[i].chain[0].position[Y]= b;
    the_membrane[i].chain[0].position[Z]= c;
    the_membrane[i].chain[0].ci.b = 0;
	the_membrane[i].chain[0].ci.m = i;
  	the_membrane[i].chain[0].user = NO;
	the_membrane[i].chain[0].typeptr = &(the_types[the_membrane[i].model_index].chain[0]); 
	//printf("first bead type:%d\n", the_membrane[i].chain[0].typeptr->type); 
	//the_membrane[i].chain[0].type = HEAD; 
    for (j = 2; j < the_membrane[i].chain_length;j++)
    		{fscanf(fp,"%lf %lf %lf ", &a2, &b2, &c2);
		the_membrane[i].chain[j].position[X]= a2;
		the_membrane[i].chain[j].position[Y]= b2;
		the_membrane[i].chain[j].position[Z]= c2;
		the_membrane[i].chain[j].ci.b = j;
		the_membrane[i].chain[j].ci.m = i;
   	the_membrane[i].chain[j].user = NO;
		the_membrane[i].chain[j].typeptr = &(the_types[the_membrane[i].model_index].chain[j]); 
		//the_membrane[i].chain[j].type = TAIL; 
			}
	
	fscanf(fp,"%lf %lf \n", &r, &l);
	
	for (j =1; j < the_membrane[i].chain_length;j++)
    		{
		int bondindex; 
		bondindex = find_bond(the_types[the_membrane[i].model_index], j-1, j); 
		if (bondindex < 0) 
		{printf("Error! Bad bondindex!\n"); 
		return;}
		if (the_types[the_membrane[i].model_index].bonds[bondindex].constrained == YES)
		{
		double radius = the_types[the_membrane[i].model_index].chain[j-1].radius + 
		the_types[the_membrane[i].model_index].chain[j].radius;
		//printf("radius:%lf\n", radius); 
		vsub(temp, the_membrane[i].chain[j].position, the_membrane[i].chain[j-1].position);
		vresize(temp, temp, radius);
		vadd(the_membrane[i].chain[j].position, the_membrane[i].chain[j-1].position, temp);}
			}
    the_membrane[i].crossings[X] = 0;
    the_membrane[i].crossings[Y] = 0;
    the_membrane[i].crossings[Z] = 0;

}

return 0;

}

int find_bond(mparticle type, int index1, int index2)
{
	int j;
	for (j = 0; j < type.num_bonds; j++)
		if ((type.bonds[j].b1 == index1) && (type.bonds[j].b2 == index2))
			return j; 
		else if ((type.bonds[j].b2 == index1) && (type.bonds[j].b1 == index2))
			return j; 
	return -1; 
}	


/*************************************************************************/
int read_model(particle *the_membrane)
{
	FILE *fp;
 	int i,j, type, lasttype;
  	extern long pnum_particles;
  	double  b;
   	int a,c;
	int b1,b2,b3; 
	int num_angles; 
	double cbend, equil; 
	fp = fopen("mc.model", "r");
 	if (fp == NULL)
	  {printf("no mc.model found\n");
		return NO; }
	for (i = 0; i < pnum_particles; i++)
	  {	fscanf(fp, "%d %d \n", &a, &b);
	  	the_membrane[i].model_index = a;
		the_membrane[i].domain_index = b; 
	  }
	fclose(fp);   		
 	return YES;
}

/*************************************************************************/
int read_top(mparticle *the_types)
{
	FILE *fp;
 	int i,j, type, lasttype;
  	extern long pnum_types;
  	double  b;
   	int a,c;
	int b1,b2,b3; 
	int num_angles, num_bonds; 
	double cbend, equil, cbond; 
	fp = fopen("mc.top", "r");
 	if (fp == NULL)
	  {printf("no mc.top found\n");
		return NO; }
	for (i = 0; i < pnum_types; i++)
	  {	fscanf(fp, "%d %lf", &a, &b);
		the_types[i].chain_length = a;
		the_types[i].cbend = b;
		the_types[i].ibindex = -1; 
		the_types[i].chain1start = -1; 
		the_types[i].chain2start = -1;
		lasttype = TAIL2 + 1;   
		the_types[i].num_chains = 2; 
		for (j = 0; j < a; j++)
			{
				fscanf(fp, "%d ", &type); 
				if (type == HEAD)
				   lasttype = HEAD; 
				if (type == INTERFACE)
				{
					the_types[i].ibindex = j;
					 lasttype = INTERFACE;
				}
				if (type == TAIL)
				{	if (((lasttype == HEAD) || (lasttype == INTERFACE))|| (j == 0))
						the_types[i].chain1start = j; 
					lasttype = TAIL; 
				}
				if (type == TAIL2)
				{
					if (lasttype == TAIL)
						the_types[i].chain2start = j; 
					lasttype = TAIL2; 
					type = TAIL; 
				}
				the_types[i].chain[j].type = type;
				if ((j == 0) && (the_types[i].chain[j].type != HEAD))
					printf("WARNING: first bead is not a head bead!\n"); 
				the_types[i].chain[j].radius = 0.5;
	//			printf("%d ", type);  
			}
		if (the_types[i].ibindex == -1)
			the_types[i].ibindex = 1; 
		if (the_types[i].chain1start == -1)
			the_types[i].chain1start = a + 1; 
		if (the_types[i].chain2start == -1)
			{the_types[i].chain2start = a + 1; 
			the_types[i].num_chains = 1; }
			
//		printf("\n"); 
		fscanf(fp, "\n"); 
		fscanf(fp, "%d\n", &num_angles); 
		the_types[i].num_angles = num_angles; 
		for (j = 0; j < num_angles; j++)
			{
				fscanf(fp, "%d %d %d %lf %lf\n", &b1, &b2, &b3, &cbend, &equil); 
				the_types[i].bondangles[j].b1 = b1; 
				the_types[i].bondangles[j].b2 = b2;
				the_types[i].bondangles[j].b3 = b3;
				the_types[i].bondangles[j].cbend = cbend;
				the_types[i].bondangles[j].equil = equil*PI;
			}
		fscanf(fp, "\n"); 
		fscanf(fp, "%d\n", &num_bonds); 
		the_types[i].num_bonds = num_bonds; 
		for (j = 0; j < num_bonds; j++)
			{
				fscanf(fp, "%d %d %d %lf %lf\n", &c, &b1, &b2, &cbond, &equil); 
				the_types[i].bonds[j].b1 = b1; 
				the_types[i].bonds[j].b2 = b2; 
				the_types[i].bonds[j].constrained = c; 
				the_types[i].bonds[j].cbond = cbond; 
				the_types[i].bonds[j].equil = equil; 
			}
	  }
	fclose(fp);   	

 	return YES;
}

int print_types(mparticle *the_types)
{
	int i,j; 
	extern long pnum_types; 
	for (i = 0; i < pnum_types; i++)
	{	
		for (j = 0; j < the_types[i].chain_length; j++)
			printf("%d ", the_types[i].chain[j].type); 
		printf("\n"); 
	}
}