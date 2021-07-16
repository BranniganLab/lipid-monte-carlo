#include <stdlib.h>
#include <math.h>
#include "structures.h"
#include "volume.h"
#include "vectorops.h"
#include "lattice.h"
#include "energy.h"
#include "accept.h"
#include "domain.h"
#include "ghosts.h"
#include "ghost_energy.h"
#include "preset.h"

int preshaping_method = POSSIBLY_CONSTRAINED;
double total_bond_energy; 
double total_bond_samples; 
 double total_bend_energy; 
double total_bend_samples; 
/************************************/
double volume_step(  particle *the_membrane,   site *the_lattice, mparticle *the_types, int *accept,
                       double delta, double previous_energy, vector box_length)
{
    /*carries out one perturbation of the volume, following surface tension algorithm */
    extern int iteration, pdo_run_preset, pgenerate_preset; 
	extern step *preset_run;
    extern double ptension;
    extern int abort_move;
    extern vector last_update; 
	int recheck; 
	vector stretch, new_box_length;
	double old_bond_energy, new_bond_energy, bond_energy_diff; 
		double old_bend_energy, new_bend_energy, bend_energy_diff; 
    double new_energy, weight, energy_diff, area_diff, old_energy;
    
	
    /***************REMOVE*******************/
    total_bond_energy = 0; 
    total_bond_samples = 0; 
       total_bend_energy = 0; 
    total_bend_samples = 0; 
    old_energy = get_total_energy(the_membrane, the_lattice, the_types, box_length); 
    if (fabs(old_energy - previous_energy) > 0.00001)
    	printf("old_energy: %lf previous_energy:%lf\n", old_energy, previous_energy); 
    //printf("total_bend_energy:%lf total_bond_energy:%lf\n", total_bend_energy, total_bond_energy); 
    old_bond_energy = total_bond_energy; 
    old_bend_energy = total_bend_energy; 
    total_bond_energy = 0; 
    total_bond_samples = 0;
        total_bend_energy = 0; 
    total_bend_samples = 0;
    vcopy(new_box_length, box_length);
	abort_move = NO_OVERLAP;
    /*generate the perturbation*/
    generate_box_shape_change(new_box_length, box_length, delta);
    reshape_membrane(the_membrane, the_lattice, the_types, box_length, new_box_length);
    vcopy(stretch, new_box_length);
	recheck =check_checklists(the_lattice, last_update, new_box_length); 
	if (recheck == YES)
	  {
		get_checklists(the_lattice, new_box_length); 
		vcopy(last_update, new_box_length); 
	  }
	//assign_all_ghosts(the_membrane, new_box_length); 
	new_energy = get_total_energy(the_membrane, the_lattice, the_types, new_box_length);
	new_bond_energy = total_bond_energy; 
	bond_energy_diff = new_bond_energy - old_bond_energy; 
		new_bend_energy = total_bend_energy; 
	bend_energy_diff = new_bend_energy - old_bend_energy; 
    energy_diff = new_energy - old_energy;
    /*add the rest (non internal energy terms) to the hamiltonian*/
    area_diff = (new_box_length[X] * new_box_length[Y]) - (box_length[X]*box_length[Y]);
    weight = (energy_diff - ptension * area_diff);
    /*test acceptance*/
	*accept = test_acceptance(weight);
    if (*accept == 1)
	  {
        /*copy the new size into the global variable*/
	// printf("dl: %lf accept:%d bond energy diff:%lf bend energy diff: %lf energy diff:%lf weight:%lf\n", new_box_length[X] - box_length[X], *accept, bond_energy_diff, bend_energy_diff, energy_diff, weight); 
        vcopy(box_length,new_box_length);
			
	  }
    else /*if it's not accepted*/
	  {
        /*clear the energy difference*/
		// printf("dl: %lf accept:%d bond energy diff:%lf bend energy diff: %lf energy diff:%lf weight:%lf\n", new_box_length[X] - box_length[X], *accept, bond_energy_diff, bend_energy_diff, energy_diff, weight); 
        energy_diff = 0;
        /*stretch membrane back to original size*/
        reshape_membrane( the_membrane, the_lattice, the_types, new_box_length, box_length);
		if (recheck == YES)
		  {
			get_checklists(the_lattice, box_length); 
			vcopy(last_update, box_length); 
		  }
	//	assign_all_ghosts(the_membrane, box_length); 
	  }
if ((pdo_run_preset == YES) || (pgenerate_preset == YES))
	     output_preset_results(preset_run[iteration]);
	iteration = iteration + 1; 
        return energy_diff;
}

/************************************/
double generate_box_shape_change(vector new_box, vector old_box, double delta)
{
	extern int pdo_run_preset, pgenerate_preset, iteration; 
   extern step *preset_run; 
    extern long seed;
    long *base = &seed;
    double size, new_length, new_height, volume,dA;
    size= delta*(ran3(base)-0.5);
    if (pdo_run_preset == YES)
        size = preset_run[iteration].scalar; 
	if (pgenerate_preset == YES)
		preset_run[iteration].scalar = size; 
    new_length = sqrt(pow(old_box[X], 2.0) + size);
   	volume = old_box[X] * old_box[Y] * old_box[Z];
	new_height = volume/(new_length*new_length);
	new_box[X] = new_length;
	new_box[Y] = new_length;
	new_box[Z] = new_height;
	dA = (new_box[X] * new_box[Y] - old_box[X] * old_box[Y]); 
	return dA; 
	
	
}
/************************************/
void reshape_membrane(particle *the_membrane, site *the_lattice, mparticle *the_types,  vector vold, vector vnew)
{
    //stretches the membrane in all directions by vnew/vold
    int i, ibindex;

    extern long pnum_particles;

    double temp;
    vector old, translate;

    particle oldparticle; 
    int mi;
    domain *domain_list[MAX_DOMAINS];
	domain *the_domain;
 	the_domain = NULL;
  	for (i = 0; i < MAX_DOMAINS; i++)
   		domain_list[i] = NULL;

    for (i = 0; i<pnum_particles; i++)
	  {
	  	mi = the_membrane[i].model_index; 
	  	ibindex = the_types[the_membrane[i].model_index].ibindex; 
		if (the_membrane[i].domain_index == 0)
			{
			CopyParticleContents(&the_membrane[i], &oldparticle); 
			vcopy(old, the_membrane[i].chain[ibindex].position);
			the_membrane[i].chain[ibindex].position[X] = the_membrane[i].chain[1].position[X] * vnew[X] / vold[X];
			the_membrane[i].chain[ibindex].position[Y] = the_membrane[i].chain[1].position[Y] * vnew[Y] / vold[Y];
			the_membrane[i].chain[ibindex].position[Z] = the_membrane[i].chain[1].position[Z] * vnew[Z]/vold[Z];
			vsub(translate, the_membrane[i].chain[ibindex].position, old);
			stretch_molecule(&the_membrane[i], &oldparticle, &(the_types[the_membrane[i].model_index]), vnew, vold); 
			
			}
		else
  			{
         			domain_list[the_membrane[i].domain_index] = the_membrane[i].domainptr;
         		}
		
	  }

   	for (i = 0; i < MAX_DOMAINS; i++)
    		{	
          	if (domain_list[i] != NULL)
			     {
            			the_domain = domain_list[i];
					vcopy(old, the_domain->center);
					the_domain->center[X] = the_domain->center[X] * vnew[X] / vold[X];
				the_domain->center[Y] = the_domain->center[Y] * vnew[Y] / vold[Y];
				the_domain->center[Z] = the_domain->center[Z] * vnew[Z] / vold[Z];
				vsub(translate, the_domain->center, old);
					translate_domain(the_domain, the_membrane, the_types, translate, vnew);
     			}
     	}

	
}

/************************************/
int stretch_molecule(particle *part, particle *oldpart, mparticle *type, vector vnew, vector vold)
{
	extern double percent_change; 
    	extern double percent_change_samples; 
	extern int preshaping_method;
	vector oldbond, newbond;  
	int  b1, b2,ibindex,j; 
	ibindex = type->ibindex; 
	for (j = 0; j < type->num_bonds; j++)
	{
		b1 = type->bonds[j].b1; 
		b2 = type->bonds[j].b2; 
		if (b1 < ibindex)
		{
			vsub_periodic(oldbond, oldpart->chain[b1].position, oldpart->chain[b2].position, vold);
			if ((type->bonds[j].constrained == YES) || (preshaping_method == DEFINITELY_CONSTRAINED))
			vadd(part->chain[b1].position, part->chain[b2].position, oldbond); 
			else
			{
				newbond[X] = oldbond[X] * vnew[X]/vold[X]; 
				newbond[Y] = oldbond[Y] * vnew[Y]/vold[Y]; 
				newbond[Z] = oldbond[Z] * vnew[Z]/vold[Z]; 
				vadd(part->chain[b1].position, part->chain[b2].position, newbond);	percent_change = percent_change + fabs(sqrt(vsquare(newbond)) - sqrt(vsquare(oldbond)))/sqrt(vsquare(oldbond)); 	percent_change_samples = percent_change_samples + 1.0; 
			}
		}
		else
		{
			vsub_periodic(oldbond, oldpart->chain[b2].position, oldpart->chain[b1].position, vold);
			if ((type->bonds[j].constrained == YES) || (preshaping_method == DEFINITELY_CONSTRAINED)) 
			{
				vadd(part->chain[b2].position, part->chain[b1].position, oldbond); 
			}
			else
			{
				newbond[X] = oldbond[X] * vnew[X]/vold[X]; 
				newbond[Y] = oldbond[Y] * vnew[Y]/vold[Y]; 
				newbond[Z] = oldbond[Z] * vnew[Z]/vold[Z]; 
				vadd(part->chain[b2].position, part->chain[b1].position, newbond); 
				percent_change = percent_change + fabs(sqrt(vsquare(newbond)) - sqrt(vsquare(oldbond)))/sqrt(vsquare(oldbond)); 
				percent_change_samples = percent_change_samples + 1.0;
			}
		}
	}
}