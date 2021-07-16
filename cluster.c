#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "structures.h"
#include "accept.h"
#include "ghost_energy.h"
#include "cluster.h"
#include "energy.h"
#include "domain.h"
#include "perturb.h"
#include "ghosts.h"
#include "lattice.h"
#include "vectorops.h"
#include "particle.h"
#include "bounds.h"
#include "out.h"
#include "volume.h"
double cluster_move(domain *the_domain, particle *the_membrane, site *the_lattice, mparticle *the_types, int domain_index, vector box_length)
{
	extern double total_system_energy; 
	int lipid_index, cluster_size, i, j, beadi, accept, ibindex; 
	double cradius, new_energy, old_energy, energy_diff; 
	double center_separation; 
	int the_cluster[MAX_CLUSTER_SIZE]; 
	int the_shell[MAX_SHELL_SIZE]; 
	int shell_size; 
	lipid_index = GetPertIndex(); 
	while (the_membrane[lipid_index].domainptr != NULL)
		lipid_index = GetPertIndex(); 
	printf("in cluster move.\n"); 
	cradius = the_domain->radius +2; 
	ibindex = the_types[the_membrane[lipid_index].model_index].ibindex; 
	center_separation = vdistance_periodic(the_membrane[lipid_index].chain[ibindex].position, the_domain->center, box_length); 
	if (center_separation < (2*cradius)*(2*cradius))
		return 0.0; 
	cluster_size = find_cluster_slow(the_membrane, lipid_index, cradius, the_cluster, box_length); 
	shell_size = find_coordination_shell_slow(the_membrane, the_domain, cradius, the_shell, box_length); 
	printf("cluster_size: %d\n", cluster_size); 
	printf("shell size: %d\n", shell_size); 
	printf("radius: %lf\n", cradius); 
	if (cluster_size == 0)
		return 0.0; 
	//old_energy = get_domain_energy(the_domain, the_membrane, the_lattice, box_length);
	//old_energy = old_energy + inter_domain_energy(the_membrane, box_length);  
	//for (i = 0; i < cluster_size; i++)
	//	old_energy = old_energy + get_molecule_energy_with_ghosts(the_cluster[i], the_membrane, the_lattice, box_length); 
	old_energy = total_system_energy; 
	cluster_swap(the_membrane, the_types, the_domain, the_cluster, cluster_size, the_shell, shell_size, box_length); 
	
	write_crd("mc.crd", the_membrane, "a"); 
	printf("wrote to crd. \n"); 
	for (i = 0; i < cluster_size; i++)
		for (j = 0; j < the_membrane[the_cluster[i]].chain_length; j++)
		{	
			update_beads_lattice_site(&(the_membrane[the_cluster[i]].chain[j]), the_lattice, box_length);	
			the_membrane[the_cluster[i]].chain[j].num_ghosts = assign_ghosts(the_membrane[the_cluster[i]].chain[j], box_length); 
		}
		
	for (i = 0; i < shell_size; i++)
		for (j = 0; j < the_membrane[the_shell[i]].chain_length; j++)
		{	
			update_beads_lattice_site(&(the_membrane[the_shell[i]].chain[j]), the_lattice, box_length);	
			the_membrane[the_shell[i]].chain[j].num_ghosts = assign_ghosts(the_membrane[the_shell[i]].chain[j], box_length); 
		}
		
	for (i = 0; i < the_domain->size; i++)
		for (beadi = 0; beadi < the_membrane[the_domain->partlist[i]].chain_length; beadi++)
			{
				update_beads_lattice_site(&(the_membrane[the_domain->partlist[i]].chain[beadi]), the_lattice, box_length); 	
				the_membrane[the_domain->partlist[i]].chain[beadi].num_ghosts = assign_ghosts(the_membrane[the_domain->partlist[i]].chain[beadi], box_length); 	
			}
			
	printf("finished updating domain: \n"); 
	//new_energy = get_domain_energy(the_domain, the_membrane, the_lattice, box_length);
	//new_energy = new_energy + inter_domain_energy(the_membrane, box_length);  
	//for (i = 0; i < cluster_size; i++)
	//	new_energy = new_energy + get_molecule_energy_with_ghosts(the_cluster[i], the_membrane, the_lattice, box_length); 
	new_energy = get_total_energy_with_ghosts(the_membrane, the_lattice, the_types, box_length); 
	energy_diff = new_energy - old_energy; 
	accept = test_acceptance(energy_diff); 
	printf("accept: %d\n", accept); 
	printf("energy diff: %lf\n", energy_diff); 

	if (accept == YES)
		return energy_diff; 
		
	cluster_swap(the_membrane, the_types, the_domain, the_cluster, cluster_size, the_shell, shell_size, box_length); 
	for (i = 0; i < cluster_size; i++)
		for (j = 0; j < the_membrane[the_cluster[i]].chain_length; j++)
		{	
			update_beads_lattice_site(&(the_membrane[the_cluster[i]].chain[j]), the_lattice, box_length);	
			the_membrane[the_cluster[i]].chain[j].num_ghosts = assign_ghosts(the_membrane[the_cluster[i]].chain[j], box_length); 
		}
		
	for (i = 0; i < shell_size; i++)
		for (j = 0; j < the_membrane[the_shell[i]].chain_length; j++)
		{	
			update_beads_lattice_site(&(the_membrane[the_shell[i]].chain[j]), the_lattice, box_length);	
			the_membrane[the_shell[i]].chain[j].num_ghosts = assign_ghosts(the_membrane[the_shell[i]].chain[j], box_length); 
		}
		
	for (i = 0; i < the_domain->size; i++)
		for (beadi = 0; beadi < the_membrane[the_domain->partlist[i]].chain_length; beadi++)
			{
				update_beads_lattice_site(&(the_membrane[the_domain->partlist[i]].chain[beadi]), the_lattice, box_length); 
				the_membrane[the_domain->partlist[i]].chain[beadi].num_ghosts = assign_ghosts(the_membrane[the_domain->partlist[i]].chain[beadi], box_length); 
			}
	//printf("calculated energy: %lf\n", get_total_energy_with_ghosts(the_membrane, the_lattice, box_length)); 
	//printf("total_system_energy: %lf\n", total_system_energy); 
	return 0.0; 
}

double hybrid_cluster_volume_move( domain *the_domain, particle *the_membrane,   site *the_lattice, mparticle *the_types, int *accept,
                       double delta,  vector box_length)
{
	extern double total_system_energy;
    extern double ptension;
    extern int abort_move;
    extern vector last_update; 
	int recheck; 
	vector stretch, new_box_length;
    double new_energy, weight, energy_diff, area_diff;
	double previous_energy; 
	int i, j,beadi, lipid_index, ibindex; 
	double cradius, center_separation; 
	int shell_size, cluster_size; 
	int the_cluster[MAX_CLUSTER_SIZE], the_shell[MAX_SHELL_SIZE]; 
	
	previous_energy = total_system_energy; 
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
	
	lipid_index = GetPertIndex(); 
	while (the_membrane[lipid_index].domainptr != NULL)
		lipid_index = GetPertIndex(); 
	printf("in cluster move.\n"); 
	cradius = the_domain->radius +2; 
	ibindex = the_types[the_membrane[lipid_index].model_index].ibindex; 
	center_separation = vdistance_periodic(the_membrane[lipid_index].chain[ibindex].position, the_domain->center, box_length); 
	if (center_separation < (2*cradius)*(2*cradius))
		return 0.0; 
	cluster_size = find_cluster_slow(the_membrane, lipid_index, cradius, the_cluster, box_length); 
	shell_size = find_coordination_shell_slow(the_membrane, the_domain, cradius, the_shell, box_length); 
 
	cluster_swap(the_membrane, the_types, the_domain, the_cluster, cluster_size, the_shell, shell_size, box_length); 
	
	write_crd("mc.crd", the_membrane, "a"); 
	printf("wrote to crd. \n"); 
	for (i = 0; i < cluster_size; i++)
		for (j = 0; j < the_membrane[the_cluster[i]].chain_length; j++)
			update_beads_lattice_site(&(the_membrane[the_cluster[i]].chain[j]), the_lattice, box_length);	
	for (i = 0; i < shell_size; i++)
		for (j = 0; j < the_membrane[the_shell[i]].chain_length; j++)	
			update_beads_lattice_site(&(the_membrane[the_shell[i]].chain[j]), the_lattice, box_length);	
	for (i = 0; i < the_domain->size; i++)
		for (beadi = 0; beadi < the_membrane[the_domain->partlist[i]].chain_length; beadi++)
				update_beads_lattice_site(&(the_membrane[the_domain->partlist[i]].chain[beadi]), the_lattice, box_length); 
	assign_all_ghosts(the_membrane, new_box_length); 
	
	new_energy = get_total_energy_with_ghosts(the_membrane, the_lattice, the_types, new_box_length);
    energy_diff = new_energy - previous_energy;
    /*add the rest (non internal energy terms) to the hamiltonian*/
    area_diff = (new_box_length[X] * new_box_length[Y]) - (box_length[X]*box_length[Y]);
    weight = (energy_diff - ptension * area_diff);
    /*test acceptance*/
	*accept = test_acceptance(weight);
    if (*accept == 1)
	  {
        /*copy the new size into the global variable*/
        vcopy(box_length,new_box_length);
	  }
    else /*if it's not accepted*/
	  {
        /*clear the energy difference*/
        energy_diff = 0;
		cluster_swap(the_membrane, the_types, the_domain, the_cluster, cluster_size, the_shell, shell_size, box_length); 

        /*stretch membrane back to original size*/
        reshape_membrane( the_membrane, the_lattice, the_types, new_box_length, box_length);
		if (recheck == YES)
		  {
			get_checklists(the_lattice, box_length); 
			vcopy(last_update, box_length); 
		  }
		for (i = 0; i < cluster_size; i++)
			for (j = 0; j < the_membrane[the_cluster[i]].chain_length; j++)
				update_beads_lattice_site(&(the_membrane[the_cluster[i]].chain[j]), the_lattice, box_length);	
		for (i = 0; i < shell_size; i++)
			for (j = 0; j < the_membrane[the_shell[i]].chain_length; j++)	
				update_beads_lattice_site(&(the_membrane[the_shell[i]].chain[j]), the_lattice, box_length);	
		for (i = 0; i < the_domain->size; i++)
			for (beadi = 0; beadi < the_membrane[the_domain->partlist[i]].chain_length; beadi++)
				update_beads_lattice_site(&(the_membrane[the_domain->partlist[i]].chain[beadi]), the_lattice, box_length); 

		assign_all_ghosts(the_membrane, box_length); 
	  }
        return energy_diff;
}

int cluster_swap(particle *the_membrane, mparticle *the_types, domain *the_domain, int the_cluster[], int cluster_size, int the_shell[], int shell_size, vector box_length)
{
	vector cluster_center, domain_center, translation; 
	int i; 
	int ibindex = the_types[the_membrane[the_cluster[0]].model_index].ibindex; 
	vcopy(cluster_center, the_membrane[the_cluster[0]].chain[ibindex].position); 
	vcopy(domain_center, the_domain->center); 
	vsub_periodic(translation, cluster_center, domain_center, box_length); 
	vmult_scal(translation, translation, -1); 
	translation[Z] = 0; 
	for (i = 0; i < cluster_size; i++)
		translate_periodic(&(the_membrane[the_cluster[i]]), &(the_types[the_membrane[the_cluster[i]].model_index]),translation, box_length); 
	vmult_scal(translation, translation, -1);
	for (i = 0; i < shell_size; i++)
		translate_periodic(&(the_membrane[the_shell[i]]),&(the_types[the_membrane[the_shell[i]].model_index]), translation, box_length); 
	translate_domain(the_domain, the_membrane, the_types, translation, box_length); 	
	vcopy(the_domain->center, cluster_center);
	
}

	

int find_cluster_slow(particle *the_membrane, int center_index, double cradius, int the_cluster[MAX_CLUSTER_SIZE], vector box_length)
{
	extern long pnum_particles; 
	int i; 
	int count = 0; 
	the_cluster[count] = center_index; 
	count = count + 1; 
	for (i = 0; i < pnum_particles; i++)
		if (i != center_index)
			if (test_same_cluster(the_membrane[i], the_membrane[center_index], cradius, box_length) == YES)
			{ 
				if (the_membrane[i].domainptr!= NULL)
					return 0; 
				if (count < MAX_CLUSTER_SIZE)
					{
						the_cluster[count] = i; 
						count = count + 1; 
					}
			}			
	return count; 			
}

int find_coordination_shell_slow(particle *the_membrane, domain *the_domain, double cradius, int the_shell[MAX_SHELL_SIZE], vector box_length)
{
	extern long pnum_particles; 
	vector radius; 
	double r; 
	int num; 
	int i,j; 
	int count = 0;  
	for (i = 0; i < pnum_particles; i++)
		{
			if (the_membrane[i].domainptr == NULL)
			{
			num = 0; 
			for (j = 0; j < the_membrane[i].chain_length; j++)
			{
				vsub_periodic(radius, the_domain->center, the_membrane[i].chain[j].position, box_length); 
				r = radius[X] * radius[X] + radius[Y] * radius[Y]; 
				if (r < cradius*cradius) num = num + 1; }
			if (num > 2)
				{
					the_shell[count] = i; 
					count = count + 1; 
				}
		}}
	return count; 		
}

int test_same_cluster(particle part1, particle part2, double cradius, vector box_length)
{
	double r; 
	int i; 
	int num = 0; 
	vector radius;  
	int ibindex = 1; 
	for (i = 0; i < part2.chain_length; i++)
//i = 1; 
		{vsub_periodic(radius, part1.chain[ibindex].position, part2.chain[i].position, box_length); 
		r = radius[X] * radius[X] + radius[Y] * radius[Y]; 
		if (r < cradius*cradius)
			num = num + 1;  }
	if (num > 2) return YES; 
	return NO; 
}
