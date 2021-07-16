#include <stdlib.h>
#include <math.h>
#include "structures.h"
#include "model.h"
#include "vectorops.h"
/*********************************************************************************/
double get_bent(particle *part, mparticle *type)
{
  	double cosangle, hold;
	int b1, b2, b3; 
	double total = 0;
	int i, chain1stop, count; 
	extern double total_bend_energy; 
	extern double total_bend_samples; 
	chain1stop = type->chain2start; 
	if (chain1stop > part->chain_length)
		chain1stop = part->chain_length; 
	count = 0;
	for (i = 0; i < type->num_angles; i++)
		{
			b1 = type->bondangles[i].b1; 
			b2 = type->bondangles[i].b2; 
			b3 = type->bondangles[i].b3; 
			cosangle = get_angle(*part, b1, b2, b3); 
			//part->cosangle[i] = cosangle; 
			//hold = type->bondangles[i].cbend*mypow(acos(cosangle) - type->bondangles[i].equil,2)/2.0; 
			//hold = -type->bondangles[i].cbend*(cos(acos(cosangle)-type->bondangles[i].equil)); 
			//hold = type->bondangles[i].cbend*(1-cos(acos(cosangle) - type->bondangles[i].equil));
			hold = type->bondangles[i].cbend*(1+cosangle); 
			//printf("cosangle[%d]:%lf\n", i, cosangle ); 
			total = total + hold; 
			count++; 
			total_bend_energy = total_bend_energy + hold; 
			total_bend_samples = total_bend_samples + 1.0; 
		} 
	/*	
	for (i = part->chain1start - 1; i < chain1stop - 2;i++)
	  { 
		cosangle = get_angle(*part, i, i+ 1, i+ 2);
		part->cosangle[count] = cosangle;
		hold = part->cbend*(1 + cosangle);
		total = total + hold;
		count++; 
	  }
	if (part->chain_length > i)
		{cosangle = get_angle(*part, part->chain1start - 1, part->chain2start, part->chain2start + 1);
		part->cosangle[count] = cosangle; 
		hold = part->cbend*(1+cosangle); 
		total = total + hold;
		count++; }
	for (i = part->chain2start; i < part->chain_length - 2; i++)
	   {
		cosangle = get_angle(*part, i, i+ 1, i+ 2);
		part->cosangle[count] = cosangle;
		hold = part->cbend*(1 + cosangle);
		total = total + hold;
		count++; 
	   }*/
	return total;	
}
/*********************************************************************************/
double get_bond_energy(particle *part, mparticle *type, vector box_length)
{
	int i; 
	double distance, energy, total; 
	extern double total_bond_energy; 
	extern double total_bond_samples; 
	total = 0.0; 
	for (i = 0; i < type->num_bonds; i++)
	{	if (type->bonds[i].constrained == NO)
	{
		distance = vdistance_periodic(part->chain[type->bonds[i].b1].position, part->chain[type->bonds[i].b2].position, box_length); 
		energy = type->bonds[i].cbond/2.0*pow(sqrt(distance) - type->bonds[i].equil, 2);
		total_bond_energy = total_bond_energy + energy;  
		total_bond_samples = total_bond_samples + 1.0; 
		total = total + energy; 
	}
		//if (energy > 0.000001) 
		//{printf("i:%d b1:%d b2:%d distance:%lf, energy:%lf\n", i, type->bonds[i].b1, type->bonds[i].b2, distance, energy); 
		//vprint(part->chain[type->bonds[i].b1].position); 
		//vprint(part->chain[type->bonds[i].b2].position); }
	}
		
	return total; 
}
/*********************************************************************************/
int adjacent(bead bead1, bead bead2)
{
	int chain1, chain2, cindex1, cindex2;
	int headstop, separation;
	separation = 4; 
	
	
	
	chain1 = 1; 
	chain2 = 1;  
	cindex1 = bead1.typeptr->chainpos; 
	cindex2 = bead2.typeptr->chainpos; 
	
	if (bead1.typeptr->chainpos < 0)
		{chain1 = 0; cindex1 = -bead1.typeptr->chainpos;} 
	else if (bead1.typeptr->chainpos != bead1.ci.b)
		chain1 = 2; 
	if (bead2.typeptr->chainpos < 0)
		{chain2 = 0; cindex2 = -bead2.typeptr->chainpos;} 
	else if (bead2.typeptr->chainpos != bead2.ci.b)
		chain2 = 2; 
		
	if ((chain1 == chain2)||((chain1 == 0) || (chain2 == 0)))
		separation = fabs(cindex1 - cindex2);
	else
		return NO; 
	//printf("separation: %d cindex1: %d cindex2:%d\n", separation, cindex1, cindex2); 
	if (separation < 4)
		return YES; 
	return NO; 
}
/*********************************************************************************/
double poly_energy_noshift(double radius, potential model, double core)
{
	double cutoff;
	int i;
	double temp=0;
	
	cutoff = (model.range) * (model.range);
/*	if (radius > cutoff)
	  {return 0.0;}*/
	for (i = 0; i< model.num_terms; i++)
      {
		temp = temp - model.coefficients[i][0] *(mypow(core*core/radius,(model.coefficients[i][1]) /2));
	  }	
	return temp;
}
/*********************************************************************************/
double poly_energy(double radius, potential model, double core)
{
	double cutoff;
	int i;
	double temp=0;
	
	cutoff = (model.range) * (model.range);
//		if (model.coefficients[0][0] != 0)
//	printf("core:%lf cutoff:%lf, coefficient:%lf, exponent: %lf, shift: %lf\n", core,
//	model.range, model.coefficients[0][0],model.coefficients[0][1],model.coefficients[0][2]); 	if (radius > cutoff)
	if (radius > cutoff)
	  {return 0.0;}
	for (i = 0; i< model.num_terms; i++)
      {
		temp = temp - model.coefficients[i][0] *(mypow(core*core/radius,(model.coefficients[i][1]) /2));
		temp = temp - model.coefficients[i][2];
      }
	
	return temp;
}
/***********************************************************************************/
double poly_attract_and_repulse(double radius, mtablesite model, double core)
{
	double cutoff;
	int i;
	double temp=0;
	double attract, repulse; 
	
	cutoff = (model.attract.range) * (model.attract.range);
	if (radius > cutoff)
	  {return 0.0;}
		attract = -model.attract.coefficient *(mypow(core*core/radius,(model.attract.exponent) /2));
		attract = attract - model.attract.shift; 
		repulse = -model.repulse.coefficient *(mypow(core*core/radius,(model.repulse.exponent) /2));
		repulse = repulse - model.repulse.shift; 
	return attract + repulse;
}


/***********************************************************************************/
int make_model_table()
{
	extern potential pmodel[]; 
	extern mtablesite pmodeltable[TAIL][TAIL];
	int i, j; 
	for (i = 0; i <= TAIL; i++)
		for (j = 0; j <= TAIL; j++)
		{
			pmodeltable[i][j].repulse.coefficient = pmodel[HSOFTCORE].coefficients[0][0]; 
			pmodeltable[i][j].repulse.exponent = pmodel[HSOFTCORE].coefficients[0][1]; 
			pmodeltable[i][j].repulse.shift = pmodel[HSOFTCORE].coefficients[0][2]; 
			pmodeltable[i][j].repulse.range = pmodel[HSOFTCORE].range; 
		}
		
	for (i = 0; i <= TAIL; i++)
		{
			pmodeltable[i][HEAD].attract.coefficient = 0; 
			pmodeltable[i][HEAD].attract.exponent = 0; 
			pmodeltable[i][HEAD].attract.shift = 0; 
			pmodeltable[i][HEAD].attract.range = pmodel[HSOFTCORE].range; 
			pmodeltable[HEAD][i].attract.coefficient = 0; 
			pmodeltable[HEAD][i].attract.exponent = 0; 
			pmodeltable[HEAD][i].attract.shift = 0; 
			pmodeltable[HEAD][i].attract.range = pmodel[HSOFTCORE].range; 
		}
	
			pmodeltable[INTERFACE][INTERFACE].attract.coefficient = pmodel[HINTERFACE].coefficients[0][0]; 
			pmodeltable[INTERFACE][INTERFACE].attract.exponent = pmodel[HINTERFACE].coefficients[0][1]; 
			pmodeltable[INTERFACE][INTERFACE].attract.shift = pmodel[HINTERFACE].coefficients[0][2]; 
			pmodeltable[INTERFACE][INTERFACE].attract.range = pmodel[HINTERFACE].range; 
			
			pmodeltable[TAIL][INTERFACE].attract.coefficient = pmodel[HTAIL].coefficients[0][0]; 
			pmodeltable[TAIL][INTERFACE].attract.exponent = pmodel[HTAIL].coefficients[0][1]; 
			pmodeltable[TAIL][INTERFACE].attract.shift = pmodel[HTAIL].coefficients[0][2]; 
			pmodeltable[TAIL][INTERFACE].attract.range = pmodel[HTAIL].range; 
			
			pmodeltable[INTERFACE][TAIL].attract.coefficient = pmodel[HTAIL].coefficients[0][0]; 
			pmodeltable[INTERFACE][TAIL].attract.exponent = pmodel[HTAIL].coefficients[0][1]; 
			pmodeltable[INTERFACE][TAIL].attract.shift = pmodel[HTAIL].coefficients[0][2]; 
			pmodeltable[INTERFACE][TAIL].attract.range = pmodel[HTAIL].range; 
			
			pmodeltable[TAIL][TAIL].attract.coefficient = pmodel[HTAIL].coefficients[0][0]; 
			pmodeltable[TAIL][TAIL].attract.exponent = pmodel[HTAIL].coefficients[0][1]; 
			pmodeltable[TAIL][TAIL].attract.shift = pmodel[HTAIL].coefficients[0][2]; 
			pmodeltable[TAIL][TAIL].attract.range = pmodel[HTAIL].range; 
	
	
}

/*************************************************************************************/
double mypow(double x, double y)
{
	int i = 0;
 	double temp = x;
  	for (i = 2; i < y + 1; i++)
   		temp =  x*temp;
	return temp;
}

double get_angle(particle part, int i, int j, int k)
{
    vector head, tail;
    double costheta;
	vsub(tail,part.chain[j].position,part.chain[i].position);
	vsub(head, part.chain[j].position, part.chain[k].position);
    costheta = (vdot(head, tail)/sqrt(vsquare(head) * vsquare(tail)));
    return costheta;
}
