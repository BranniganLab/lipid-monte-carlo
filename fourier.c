#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "structures.h"
#include "ghost_energy.h"
#include "ghosts.h"
#include "fourier.h"
#include "vectorops.h"
#include "energy.h"
#include "debug.h"
#include "accept.h"

int num_modesx; 
int num_modesy; 
int max_mode_radius_squared; 

/************************************/
double fourier_step(  particle *the_membrane,   site *the_lattice, mparticle *the_types, int *accept,
                       movesize *fdelta, double previous_energy, vector box_length)
{
    /*carries out one perturbation of the volume, following surface tension algorithm */
	
    extern double ptension;
    extern int abort_move;
    extern vector last_update; 
	int recheck, p, q; 
    double new_energy, weight, energy_diff, shift, phase, slow_energy;
	int r2; 
	
	abort_move = NO_OVERLAP;
    /*generate the perturbation*/
    shift = generate_fourier_move(&p, &q, &phase, fdelta); 
	fourier_shift(the_membrane, p, q, phase, shift, box_length); 
	//assign_all_ghosts(the_membrane, box_length); 
	new_energy = get_total_energy(the_membrane, the_lattice, the_types, box_length);
	//slow_energy = slow_total_energy(the_membrane, box_length); 
	//if (fabs(slow_energy - new_energy)>0.0001)
	//	printf("error! slow_energy:%lf new_energy:%lf\n", slow_energy, new_energy); 
    energy_diff = new_energy - previous_energy;
    weight = (energy_diff);
    /*test acceptance*/
	*accept = test_acceptance(weight);
	r2 = p*p + q*q; 
	fdelta[r2].attempts = fdelta[r2].attempts + 1; 
	 
    if (*accept == 1)
	  {
		fdelta[r2].accepts = fdelta[r2].accepts + 1; 
		printf("fourier move accepted: p: %d q:%d delta: %lf, energy_diff: %lf\n", p, q, shift, energy_diff); 
		return energy_diff; 
		}
    else /*if it's not accepted*/
	  {
	  printf("fourier move rejected: p: %d q:%d delta: %lf, energy diff: %lf\n", p, q, shift, energy_diff); 
        /*clear the energy difference*/
        energy_diff = 0;
		fourier_shift(the_membrane, p, q, phase, -shift, box_length); 
		//assign_all_ghosts(the_membrane, box_length); 
	  }
	return energy_diff;
}

/*******************************************************/
double generate_fourier_move(int *p, int *q, double *phase, movesize *fdelta)
{
	extern long seed; 
	long *base = &seed; 
	extern int num_modesx, num_modesy, max_mode_radius_squared; 
	double delta, shift; 
	int r2; 
	
	*p = (int)(ran3(base)* num_modesx); 
	*q = (int)(ran3(base)* num_modesy); 
	if (*p == num_modesx) *p = num_modesx - 1; 
	if (*q == num_modesy) *q = num_modesy - 1; 

	while (((*p)*(*p) + (*q)*(*q) > max_mode_radius_squared) || ((*p)*(*p) + (*q)*(*q)==0))
	{
		*p = (int)(ran3(base)* num_modesx); 
		*q = (int)(ran3(base)* num_modesy); 
		if (*p == num_modesx) *p = num_modesx - 1; 
		if (*q == num_modesy) *q = num_modesy - 1; 
	}
	r2 = (*p)*(*p) + (*q)*(*q); 
	delta = fdelta[r2].size; 
	shift = delta*(ran3(base)-0.5); 
	
  	*phase = (ran3(base))/(2 * PI);
	
	return shift; 
}

/***************************************************/
void fourier_shift(particle *the_membrane, int p, int q, double phase, double shift, vector box_length)
{
	int i; 
	vector r, translation;
	extern long pnum_particles; 
	extern int num_modesx, num_modesy; 	
	vclear(translation); 
	for (i = 0; i < pnum_particles; i++)
	{
		vcopy(r, the_membrane[i].chain[0].position); 
		translation[Z] = shift*cos(2 * PI * (r[X] * p/num_modesx 
+ r[Y]*q/num_modesy) + phase); 
		translate_periodic(&(the_membrane[i]), translation, box_length); 
	}
}

/***********************************************/
int read_in_fdelta(movesize *fdelta)
{
	int i,r;
	double d;  
	FILE *fp; 
	extern int max_mode_radius_squared; 
	

	for (i = 0; i <= max_mode_radius_squared; i++)
		{fdelta[i].size = 00000.0;
		fdelta[i].accepts = 0;  
		fdelta[i].attempts = 0; } 

	fp = fopen("fourier.pref","r"); 
	if (fp == NULL)
		{
			printf("I can't find fourier.pref - generating now\n"); 
			generate_fdelta(fdelta);
			write_fdelta(fdelta);
			return NO; 
		}
	while (!feof(fp))
	{
		fscanf(fp,"%d %lf\n", &r, &d); 
		
		if (r <= max_mode_radius_squared)
			fdelta[r].size = d; 
	
	}
	fclose(fp); 
	return YES; 	
}
/*************************************/

void generate_fdelta(movesize *fdelta)
{
	int i,j, r2; 
	extern int num_modesx, num_modesy, max_mode_radius_squared; 
	for (i = 0; i < num_modesx; i++)
		for (j = 0; j < num_modesy; j++)
			{
				r2 = i*i + j*j; 
				if (r2 <= max_mode_radius_squared)
					fdelta[r2].size = 0.02; 
			}
}

/************************************/
void write_fdelta(movesize *fdelta)
{
	int i; 
	extern int max_mode_radius_squared; 
	FILE *fp; 
	fp = fopen("fourier.pref", "w"); 
	for (i = 0; i < max_mode_radius_squared; i++)
		if (fdelta[i].size > 0)
			{fprintf(fp, "%d %lf\n", i, fdelta[i].size); }
	fclose(fp); 
}

/***********************************/
void adjust_fdelta(movesize *fdelta)
{
	int i; 
	extern int max_mode_radius_squared; 
	double ratio; 
	
	for (i = 0; i < max_mode_radius_squared; i++)
		if (fdelta[i].attempts > 0)
			{
				ratio = ((double)(fdelta[i].accepts))/((double)fdelta[i].attempts); 
				if (ratio < 0.2)
					fdelta[i].size = fdelta[i].size * 0.95; 
				else if (ratio > 0.5)
					fdelta[i].size = fdelta[i].size /0.95; 
				fdelta[i].accepts = 0; 
				fdelta[i].attempts = 0; 
			}
}

/********************************/
void set_fourier_constants(vector box_length)
{
	num_modesx = (int)(box_length[X]/2.0); 
	num_modesy = (int)(box_length[Y]/2.0); 
	max_mode_radius_squared = 20; 
	printf("num_modesx: %d num_modesy:%d, max:%d\n", num_modesx, num_modesy, max_mode_radius_squared); 
}
