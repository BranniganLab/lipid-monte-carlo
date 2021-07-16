#include <stdlib.h>
#include <stdio.h>
#include "structures.h"
#include "energy.h"
#include "ghost_energy.h"
#include "model.h"
#include "vectorops.h"
#include "domain.h"
#include <math.h>

/***********************************************/
double get_bead_energy_with_ghosts(bead the_bead, particle *the_membrane, site *the_lattice, vector box_length)
{
	extern int num_to_check;
 extern int pnumber_user_interactions;
 extern upotential *puser_potentials;
	int lattice_site, i, to_check,j, found_ghost; 
	cintnode *ref_bead; 
	chainindex cci; 
	double hold; 
	double total = 0.0; 

	lattice_site = the_bead.site_index; 
	for (i = 0; i < num_to_check; i++)
	  {
		to_check = the_lattice[lattice_site].checklist[i]; 
		ref_bead = the_lattice[to_check].bead_list; 
		while (ref_bead != NULL)
		  {
			cci.m = ref_bead->contents.m; 
			cci.b = ref_bead->contents.b; 
			if (cci.m != the_bead.ci.m)
			{
			hold = 0.5*get_bead_pair_energy_with_ghosts(the_bead,the_membrane[ref_bead->contents.m].chain[ref_bead->contents.b], box_length); 
			total = total + hold; 
			found_ghost = NO; 
			if (hold == 0.0)
				for (j = 0; j < the_membrane[cci.m].chain[cci.b].num_ghosts; j++)
				{
					hold = .25 * get_bead_pair_energy_with_ghosts(the_bead, the_membrane[cci.m].chain[cci.b].ghosts[j], box_length); 
					total = total + hold; 
					//if (found_ghost == NO) total = total + hold; 
					//if (hold != 0) found_ghost = YES; 
				}
				for (j = 0; j < the_bead.num_ghosts; j++)
				{
					total = total + .25 * get_bead_pair_energy_with_ghosts(the_bead.ghosts[j], the_membrane[cci.m].chain[cci.b], box_length); 
				}
			}
			ref_bead = ref_bead->next; 
		  }
	  }
     if (the_bead.user == YES)
     	for (i = 0; i < pnumber_user_interactions; i++)
      		{if ((puser_potentials[i].ci1.m == the_bead.ci.m) && (puser_potentials[i].ci1.b == the_bead.ci.b))
      			total = total + get_bead_user_energy(i,the_membrane,  box_length);
         		else if ((puser_potentials[i].ci2.m == the_bead.ci.m) && (puser_potentials[i].ci2.b == the_bead.ci.b))
                 	total = total + get_bead_user_energy(i, the_membrane, box_length);}
	return total; 
}

/***********************************************/
double get_bead_pair_energy_with_ghosts(bead bead1, bead bead2, vector box_length)
{
	extern potential pmodel[]; 
	extern int puser; 
	double core, r, ex_volume, interf,tail, user_energy; 
	double tmp; 
	extern mtablesite pmodeltable[3][3]; 
	if (bead1.ci.m == bead2.ci.m) 
		if (adjacent(bead1, bead2)==YES)
			return 0.0; 
	core = bead1.typeptr->radius + bead2.typeptr->radius; 
	r = vdistance(bead1.position, bead2.position); 
	
	//tmp = poly_attract_and_repulse(r, pmodeltable[bead1.type][bead2.type], core); 
	
	ex_volume = poly_energy(r, pmodel[HSOFTCORE], core);
	if ((bead1.typeptr->type == HEAD) || (bead2.typeptr->type == HEAD)) return ex_volume;
	if ((bead1.typeptr->type == INTERFACE) && (bead2.typeptr->type == INTERFACE))
	  {
		interf = poly_energy(r, pmodel[HINTERFACE], core); 
		return interf + ex_volume;
	  }
	tail = poly_energy(r, pmodel[HTAIL], core); 
	return tail + ex_volume;
	return tmp; 
}



/***********************************************/
double get_molecule_pair_energy_with_ghosts(particle part1, particle part2, site *the_lattice, vector box_length)
{
/*does not include user energies. bad. */
	int i, j;
	double total = 0;
	for (i = 0; i < part1.chain_length; i++)
		for (j = 0; j < part2.chain_length; j++)
			{total = total + get_bead_pair_energy_with_ghosts(part1.chain[i], part2.chain[j], box_length);}
	return total;
}
/***********************************************/

double get_molecule_energy_with_ghosts(int index, particle *the_membrane, site *the_lattice, mparticle *the_types, vector box_length)
{
	double bending; 
	double total = 0; 
	int i; 
	bending = get_bent(&the_membrane[index], &the_types[the_membrane[index].model_index]); 
	for (i = 0; i < the_membrane[index].chain_length; i++)
		total = total + get_bead_energy_with_ghosts(the_membrane[index].chain[i], the_membrane, the_lattice, box_length); 
	return bending + total; 
}

/***********************************************/
double get_total_energy_with_ghosts(particle *the_membrane, site *the_lattice, mparticle *the_types, vector box_length)
{
	int i; 
	double total = 0; 
	extern long pnum_particles; 
	for (i = 0; i < pnum_particles; i++)
		total = total + get_molecule_energy_with_ghosts(i, the_membrane, the_lattice, the_types, box_length); 
	total = total + inter_domain_energy(the_membrane, box_length); 
	
	return total; 
}

/***********************************************/
double get_total_energy_with_ghosts_slow(particle *the_membrane, site *the_lattice, mparticle *the_types, vector box_length)
{
	int i,j,k,l,m; 
	double total = 0; 
	double scale; 
	double hold; 
	extern long pnum_particles; 
	
	for (i = 0; i < pnum_particles; i++)
		total = total + get_bent(&the_membrane[i], &the_types[the_membrane[i].model_index]);
	
	for (i = 0; i < pnum_particles; i++)
		for (j = i; j < pnum_particles; j++)
			for (k = 0; k < the_membrane[i].chain_length; k++)
				for (l = 0; l < the_membrane[j].chain_length; l++)
				{
					if (i == j) scale = 0.5; else scale = 1.0; 
					hold = scale*get_bead_pair_energy_with_ghosts(the_membrane[j].chain[l],the_membrane[i].chain[k], box_length); 
			 		total = total + hold; 
					for (m = 0; m < the_membrane[i].chain[k].num_ghosts; m++)
					{
						hold = scale*0.5*get_bead_pair_energy_with_ghosts(the_membrane[j].chain[l], the_membrane[i].chain[k].ghosts[m], box_length); 
						total = total + hold; 
					}
					for (m = 0; m < the_membrane[j].chain[l].num_ghosts; m++)
					{
						hold = scale*0.5*get_bead_pair_energy_with_ghosts(the_membrane[i].chain[k], the_membrane[j].chain[l].ghosts[m], box_length); 
						total = total + hold; 
					}
				}
	total = total + inter_domain_energy(the_membrane, box_length); 
	
	return total; 
}
