#include <stdlib.h>
#include <stdio.h>
#include "structures.h"
#include "energy.h"
#include "model.h"
#include "vectorops.h"
#include "domain.h"
#include <math.h>

int tail_pairs = 0; 

/***********************************************/
double get_bead_energy(bead the_bead, particle *the_membrane, site *the_lattice, vector box_length)
{
	extern int num_to_check;
 extern int pnumber_user_interactions;
 extern upotential *puser_potentials;
	int lattice_site, i, to_check; 
	cintnode *ref_bead; 
	double total = 0.0; 
	
	lattice_site = the_bead.site_index; 
	for (i = 0; i < num_to_check; i++)
	  {
		to_check = the_lattice[lattice_site].checklist[i]; 
		ref_bead = the_lattice[to_check].bead_list; 
		while (ref_bead != NULL)
		  {
			total = total + 0.5*get_bead_pair_energy(the_bead,the_membrane[ref_bead->contents.m].chain[ref_bead->contents.b], box_length); 
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
double get_molecule_self_energy(particle part, vector box_length)
{
	int k,l; 
	double total = 0.0; 
	for (k = 0; k < part.chain_length; k++)
		for (l = k+1; l < part.chain_length; l++)
				total = total + get_bead_pair_energy(part.chain[k], part.chain[l], box_length);
	return total; 
}
/***********************************************/
double get_bead_pair_energy(bead bead1, bead bead2, vector box_length)
{
	extern potential pmodel[]; 
	extern int puser; 
	extern int tail_pairs;
	extern double global_r; 
	double core, r, ex_volume, intf,tail, user_energy; 
	double scale = 1.0; 
	FILE *fp; 

	if (bead1.ci.m == bead2.ci.m) 
	{
		if (adjacent(bead1, bead2)==YES)
			{ 
			return 0.0; }
		
		//scale = 0.5;
		}

	core = bead1.typeptr->radius + bead2.typeptr->radius; 
	r = vdistance_periodic(bead1.position, bead2.position,box_length); 
	global_r = r; 
	ex_volume = scale*poly_energy(r, pmodel[HSOFTCORE], core);

	if ((bead1.typeptr->type == HEAD) || (bead2.typeptr->type == HEAD)) return ex_volume;
	if ((bead1.typeptr->type == INTERFACE) && (bead2.typeptr->type == INTERFACE))
	  {
		intf = scale*poly_energy(r, pmodel[HINTERFACE], core); 
		return intf + ex_volume;
	  }
	tail = scale*poly_energy(r, pmodel[HTAIL], core);
	
	if (bead1.ci.m == bead2.ci.m) 
	{
		printf("tail beads on the same molecule interacting: %d %d %lf \n", bead1.ci.b, bead2.ci.b, sqrt(r)); 
		adjacent(bead1, bead2);
		}
	tail_pairs++;  
	
	return tail + ex_volume;
}

/***********************************************/

double get_bead_user_energy(int index,particle *the_membrane, vector box_length)
{
	double user_energy = 0; 
	vector my_radius; 
	extern int pnumber_user_interactions; 
	extern upotential *puser_potentials;
 chainindex ci1, ci2;
double rho; 
    ci1.m = puser_potentials[index].ci1.m;
    ci2.m =  puser_potentials[index].ci2.m;
     ci1.b =  puser_potentials[index].ci1.b;
    ci2.b =  puser_potentials[index].ci2.b; 
     vsub_periodic(my_radius, the_membrane[ci1.m].chain[ci1.b].position, 
the_membrane[ci2.m].chain[ci2.b].position, box_length);
	rho = my_radius[X] * my_radius[X] + my_radius[Y] * my_radius[Y]; 
	user_energy = calculate_user_energy(rho, 
puser_potentials[index]);
	return user_energy; 
}

/***********************************************/
double calculate_user_energy(double radius, upotential the_potential)
{
	if (the_potential.type == HARMONIC)
		return 0.5*the_potential.coefficient*pow(sqrt(radius) - the_potential.r0,the_potential.exponent); 
	else if ((the_potential.type == LJ) && (sqrt(radius) < the_potential.range))
		return the_potential.coefficient*pow(the_potential.r0/sqrt(radius), the_potential.exponent); 
}

/***********************************************/
double get_molecule_pair_energy(particle part1, particle part2, site *the_lattice, vector box_length)
{
/*does not include user energies. bad. */
	int i, j;
	double total = 0;
	for (i = 0; i < part1.chain_length; i++)
		for (j = 0; j < part2.chain_length; j++)
			{total = total + get_bead_pair_energy(part1.chain[i], part2.chain[j], box_length);}
	return total;
}
/***********************************************/

double get_molecule_energy(int index, particle *the_membrane, site *the_lattice, mparticle *the_types, vector box_length)
{
	double bending,hold, bonding; 
	double total = 0; 
	int i; 
	bending = get_bent(&the_membrane[index], &(the_types[the_membrane[index].model_index])); 
	bonding = get_bond_energy(&the_membrane[index], &(the_types[the_membrane[index].model_index]), box_length); 
	//printf("bending: %lf, bonding %lf\n", bending, bonding); 
	for (i = 0; i < the_membrane[index].chain_length; i++)
		{
			hold = get_bead_energy(the_membrane[index].chain[i], the_membrane, the_lattice, box_length); 
		//	if (fabs(hold) > 0)
		//		printf("i:%d, hold:%lf\n", i, hold); 
			total = total + hold; 
		}
	//total = total - 0.5*get_molecule_self_energy(the_membrane[index], box_length); 
	//printf("total:%lf\n", bending+total); 
	return bending + bonding+ total; 
}

/***********************************************/
double get_total_energy(particle *the_membrane, site *the_lattice, mparticle *the_types, vector box_length)
{
	int i; 
	double total = 0; 
	extern int tail_pairs; 
	double hold; 
	FILE *fp; 
	extern long pnum_particles;
	tail_pairs = 0;  

	fp = fopen("energylist","w"); 
	for (i = 0; i < pnum_particles; i++)
	{
	
		hold = get_molecule_energy(i, the_membrane, the_lattice, the_types, box_length); 
		total = total + hold; 
		fprintf(fp, "%d %lf\n", i, hold); 
		}
	total = total + inter_domain_energy(the_membrane, box_length); 
	//printf("total energy:%lf\n", total); 
	//printf("tail pairs: %d\n", tail_pairs); 
	return total; 
}
/*********************************************/
