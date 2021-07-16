#include <stdlib.h>
#include <stdio.h>
#include "structures.h"
#include "energy.h"
#include "model.h"
#include "vectorops.h"
#include "domain.h"
#include "lattice.h"
#include <math.h>
#include "perturb.h"
#include "accept.h"
#include "bounds.h"
#include "ghosts.h"
#include "ghost_energy.h"
/*******************************************************************/
double step_domain(domain *the_domain, particle *the_membrane, site *the_lattice, mparticle *the_types, double delta, vector box_length)
{
	extern double domain_attempts, domain_accepts, total_system_energy; 
	vector old_position, old_arrow;
	double old_energy, new_energy, energy_diff; 
	double old_intra, new_intra, final_intra; /*should all be the same*/ 
	int accept, i, beadi; 
	vector translation, rotation; 
	particle *domain_molecules; 
	
//	printf("Starting: calculated energy: %lf\n", get_total_energy_with_ghosts(the_membrane, the_lattice, box_length)); 
//	printf("total_system_energy: %lf\n", total_system_energy); 
	
	vcopy(old_position, the_domain->center); 
	vcopy(old_arrow, the_domain->arrow); 
	domain_attempts = domain_attempts + 1; 
	old_energy = get_domain_energy(the_domain, the_membrane, the_lattice, the_types, box_length);
	old_energy = old_energy + inter_domain_energy(the_membrane, box_length);  
	domain_molecules = store_domain_positions(the_domain, the_membrane); 

	perturb_domain(the_domain, the_membrane, the_types,  translation, rotation, delta, box_length);
	for (i = 0; i < the_domain->size; i++)
		for (beadi = 0; beadi < the_membrane[the_domain->partlist[i]].chain_length; beadi++)
		{update_beads_lattice_site(&(the_membrane[the_domain->partlist[i]].chain[beadi]), the_lattice, box_length); 
		the_membrane[the_domain->partlist[i]].chain[beadi].num_ghosts = assign_ghosts(the_membrane[the_domain->partlist[i]].chain[beadi], box_length); 	
		}
	new_energy = get_domain_energy(the_domain, the_membrane, the_lattice, the_types, box_length); 
	new_energy = new_energy + inter_domain_energy(the_membrane, box_length);  
	energy_diff = 2*(new_energy - old_energy); 
	accept = test_acceptance(energy_diff); 
	domain_accepts = domain_accepts + accept;  
	if (accept == YES)
	{	
		if (domain_molecules != NULL)
			free(domain_molecules);
		return energy_diff;} 
		
	if (accept == NO)
	{
		vcopy(the_domain->center, old_position);
		vcopy(the_domain->arrow, old_arrow); 
		retrieve_domain_positions(the_domain, the_membrane, domain_molecules); 
		for (i = 0; i < the_domain->size; i++)
		for (beadi = 0; beadi < the_membrane[the_domain->partlist[i]].chain_length; beadi++)
		{update_beads_lattice_site(&(the_membrane[the_domain->partlist[i]].chain[beadi]), the_lattice, box_length); 
		the_membrane[the_domain->partlist[i]].chain[beadi].num_ghosts = assign_ghosts(the_membrane[the_domain->partlist[i]].chain[beadi], box_length); 	
		}	
		if (domain_molecules != NULL)
			free(domain_molecules);		
		return 0.0; 
	}
}

/**************************************/
//~ int fill_domain(domain *the_domain, particle *the_membrane, site *the_lattice, vector box_length)
//~ {
	//~ vector tmpx, tmpy, center, position, xaxis, yaxis, n, newy, newx, chain_shift; 
	//~ int i, bead; 
	//~ double a, b; 

	//~ vclear(xaxis); xaxis[X] = 1; 
	//~ vclear(yaxis); yaxis[Y] = 1; 
	//~ vcopy(n, the_domain->arrow); 
	//~ vcross(newy, n, xaxis); 
	//~ vunit(newy, newy); 
	//~ vcross(newx, newy, n); 
	//~ vunit(newx, newx); 
	
	//~ vunit(chain_shift, the_domain->arrow); 
	//~ i = 0; 
	//~ while (i < the_domain->size)
		//~ {
			//~ a = the_domain->radius * cos(the_domain->spacing * i/2.0); 
			//~ b = the_domain->radius * sin(the_domain->spacing * i/2.0); 
			//~ vmult_scal(tmpx, newx, a); 
			//~ vmult_scal(tmpy, newy, b); 
			//~ vadd(center, tmpx, tmpy);
			//~ vadd(center, the_domain->center, center);  
			//~ for (bead = 0; bead < the_domain->chain_length; bead++)
			//~ {
				//~ vmult_scal(position, chain_shift, -0.5 + (the_domain->chain_length - bead)); 
				//~ vadd(position, center, position); 
				//~ vcopy(the_membrane[the_domain->partlist[i]].chain[bead].position, position); 
				//~ vmult_scal(position, chain_shift, +0.5 - (the_domain->chain_length - bead)); 
				//~ vadd(position, center, position); 
				//~ vcopy(the_membrane[the_domain->partlist[i+1]].chain[bead].position, position); 
			//~ }
			//~ i = i + 2; 
		//~ }		

	//~ for (i = 0; i < the_domain->size; i++)
	//~ {
		//~ for (bead = 0; bead < the_domain->chain_length; bead++)
			//~ update_beads_lattice_site(&(the_membrane[the_domain->partlist[i]].chain[bead]), the_lattice, box_length); 
			
	//~ }
//~ }

/*****************************************/
int translate_domain(domain *the_domain, particle *the_membrane, mparticle *the_types, vector  translation, vector box_length)
{
	int i; 
	for (i = 0; i < the_domain->size; i++)
		translate_periodic(&(the_membrane[the_domain->partlist[i]]), &(the_types[the_membrane[the_domain->partlist[i]].model_index]), translation, box_length); 
}
/******************************************/
int rotate_domain(domain *the_domain, particle *the_membrane, vector  rotation, vector box_length)
{
	vector new_arrow, old_arrow, unitn, newr, chain_n,r; 
	int i, j, last_bead; 

	vcopy(old_arrow, the_domain->arrow); 
	rotate_vect(the_domain->arrow, rotation); 
	
	vunit(unitn, new_arrow); 
	for (i = 0; i < the_domain->size; i++)
		{
			last_bead = the_membrane[the_domain->partlist[i]].chain_length - 1; 
			vsub_periodic(r,  the_membrane[the_domain->partlist[i]].chain[last_bead].position,the_domain->center, box_length); 
			rotate_vect(r, rotation);
			vcopy(newr, r); 
			vsub_periodic(chain_n, the_membrane[the_domain->partlist[i]].chain[last_bead-1].position,
					the_membrane[the_domain->partlist[i]].chain[last_bead].position, box_length); 
			vadd(the_membrane[the_domain->partlist[i]].chain[last_bead].position,  the_domain->center, newr); 
			rotate_vect(chain_n, rotation); 
			for (j = last_bead; j > 0; j--)
				vadd(the_membrane[the_domain->partlist[i]].chain[j-1].position,  the_membrane[the_domain->partlist[i]].chain[j].position, chain_n); 
		}
}

/**************************************************************/
int initialize_domains(particle *the_membrane, site *the_lattice, mparticle *the_types, vector box_length)
{
	extern vector pdfixed; 
	extern int pdomain_fixed; 
	extern long pnum_particles; 
	domain *the_domain; 
	int i,r, already_allocated;
	double a,b,c, x, y, z; 
	domain *domain_list[MAX_DOMAINS] ;
 	FILE *fp;

	for (i = 0; i < MAX_DOMAINS; i++)
    		domain_list[i] = NULL;
	for (i = 0; i < pnum_particles; i++)
		if (the_membrane[i].domain_index == 0)
			the_membrane[i].domainptr = NULL; 
		else if (domain_list[the_membrane[i].domain_index] == NULL)
  			{
         				the_domain = (domain *)(malloc(sizeof(domain)));
     				the_membrane[i].domainptr = the_domain;
         				domain_list[the_membrane[i].domain_index] = the_domain;
             			the_domain->partlist[0] = i;
             			the_domain->size = 1;
         					
           }
         else
         		{
           			the_domain = domain_list[the_membrane[i].domain_index];
              		the_membrane[i].domainptr = the_domain;	
                		the_domain->partlist[the_domain->size] = i;
                  		the_domain->size = the_domain->size + 1;
                }

     for (i = 0; i < MAX_DOMAINS; i++)
     	if (domain_list[i] != NULL)
			get_domain_center_and_tilt(domain_list[i], the_membrane, the_types, box_length);
	if ((pdomain_fixed > 0) && (pdomain_fixed < MAX_DOMAINS))
		vcopy(pdfixed, domain_list[pdomain_fixed]->center); 
}

/*************************************************************/
double get_domain_energy(domain *the_domain, particle *the_membrane, site *the_lattice, mparticle *the_types, vector box_length)
{
	int i; 
	double total = 0; 
	for (i = 0; i < the_domain->size; i++)
		total = total + get_molecule_energy(the_domain->partlist[i], the_membrane, the_lattice, the_types, box_length); 
	return total; 
} 

/*************************************************************/
double inter_domain_energy(particle *the_membrane, vector box_length)
{
	int i, j; 
	extern long pnum_particles; 
	double radius; 
	extern double pumbrella_coeff, pumbrella_center; 
	double total = 0; 
	domain *domain_list[MAX_DOMAINS];
	
    	for (i = 0; i < MAX_DOMAINS; i++)
     		domain_list[i] = NULL;

    for (i = 0; i<pnum_particles; i++)
		if (the_membrane[i].domain_index != 0)
  			{if (the_membrane[i].domainptr != NULL)
				domain_list[the_membrane[i].domain_index] = the_membrane[i].domainptr;
               else
               	printf("Non zero domain index has NULL pointer!\n");
               }
   	for (i = 0; i < MAX_DOMAINS; i++)
    	for (j= i+1; j < MAX_DOMAINS; j++)
			if ((domain_list[i]!=NULL) && (domain_list[j]!=NULL))
				{
					radius = sqrt(vdistance_periodic(domain_list[i]->center, domain_list[j]->center, box_length));
					total = total + pumbrella_coeff*pow(radius - pumbrella_center, 2.0);
				}
	return total; 
}

/*************************************************************/
double get_intra_domain_energy(domain *the_domain, particle *the_membrane, site *the_lattice, vector box_length)
{
/**********fix to include user energies**********/
	int i, j; 
	double total = 0; 
	for (i = 0; i < the_domain->size; i++)
		for (j = i + 1; j < the_domain->size; j++)
			total = total + get_molecule_pair_energy(the_membrane[the_domain->partlist[i]], 
					the_membrane[the_domain->partlist[j]], the_lattice, box_length); 
	return total; 
} 
/**************************************************************/
int get_domain_center_and_tilt(domain *the_domain, particle *the_membrane,  mparticle *the_types, vector box_length)
{
/****only works for domain not crossing boundary conditions*****/
	int i;
	double largest_r; 
	int ibindex; 
	vector center, radius, chain_n, arrow; 
	vclear(center); 
	vclear(arrow); 
	for (i = 0; i < the_domain->size; i++)
		{
			ibindex = the_types[the_membrane[the_domain->partlist[i]].model_index].ibindex; 
			vadd(center, center, the_membrane[the_domain->partlist[i]].chain[0].position); 
			vsub(chain_n, the_membrane[the_domain->partlist[i]].chain[0].position, the_membrane[the_domain->partlist[i]].chain[ibindex].position); 
			chain_n[Z] = fabs(chain_n[Z]); 
			vadd(arrow, arrow, chain_n); 
		}
	vmult_scal(center, center, 1.0/((double)the_domain->size));
	vmult_scal(arrow, arrow, 1.0/((double)the_domain->size));
	vunit(arrow, arrow); 
	vcopy(the_domain->arrow, arrow); 
	vcopy(the_domain->center, center); 
	largest_r = 0; 
	for (i = 0; i < the_domain->size; i++)
	{
		vsub_periodic(radius, center, the_membrane[the_domain->partlist[i]].chain[0].position, box_length); 
		if (radius[X]*radius[X] + radius[Y]*radius[Y] > largest_r)
			largest_r = radius[X]*radius[X] + radius[Y]*radius[Y]; 
	}
	the_domain->radius = sqrt(largest_r); 	
}

/********************************************************/
particle *store_domain_positions(domain *the_domain, particle *the_membrane)
{
	int i,j; 
	particle *domain_molecules; 
	domain_molecules = (particle *)(malloc(the_domain->size * sizeof(particle))); 
	for (i = 0; i < the_domain->size; i++)
		for (j = 0; j < the_membrane[the_domain->partlist[i]].chain_length; j++)
		{
			vcopy(domain_molecules[i].chain[j].position, the_membrane[the_domain->partlist[i]].chain[j].position);
		}
	return domain_molecules; 
}

/********************************************************/
int retrieve_domain_positions(domain *the_domain, particle *the_membrane, particle *domain_molecules)
{
	int i,j; 
	for (i = 0; i < the_domain->size; i++)
		for (j = 0; j < the_membrane[the_domain->partlist[i]].chain_length; j++)
		{
			vcopy(the_membrane[the_domain->partlist[i]].chain[j].position,domain_molecules[i].chain[j].position);
		}
}
