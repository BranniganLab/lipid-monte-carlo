#include <stdlib.h>
#include "structures.h"
#include "energy.h"
#include "model.h"
#include "debug.h"
#include "vectorops.h"
#include "ghost_energy.h"
#include "volume.h"

double global_r = 0; 

/********************************************************************/
int check_bond_energy_difference(particle *the_membrane, mparticle *the_types, vector box_length)
{
	extern long pnum_particles; 
	int i;
	double hypo, before, after, diff, dA, bond_hypo, before_bond, after_bond, bond_diff; 
	particle oldpart; 
	vector new_box_length;  
	
	dA = generate_box_shape_change(new_box_length, box_length, 0.01);
	
	for (i = 0; i < pnum_particles; i++)
		{
			hypo = get_hypothetical_da_with_bond(&(the_membrane[i]), &(the_types[the_membrane[i].model_index]), box_length, dA);
			bond_hypo = get_hypothetical_bond_stretch(&(the_membrane[i]), &(the_types[the_membrane[i].model_index]), box_length, dA);
			before = get_bent(&(the_membrane[i]), &(the_types[the_membrane[i].model_index]));
			before_bond = get_bond_energy(&(the_membrane[i]), &(the_types[the_membrane[i].model_index]), box_length); 
			CopyParticleContents(&the_membrane[i], &oldpart);
			stretch_molecule(&(the_membrane[i]), &oldpart, &(the_types[the_membrane[i].model_index]), new_box_length, box_length); 
			after = get_bent(&(the_membrane[i]), &(the_types[the_membrane[i].model_index])); 
			after_bond = get_bond_energy(&(the_membrane[i]), &(the_types[the_membrane[i].model_index]), new_box_length);
			diff = after - before; 
			bond_diff = after_bond - before_bond; 
			if (fabs(diff - hypo)>0.0001)
				printf("diff:%lf, hypo:%lf, percent: %lf, ratio:%lf\n", diff, hypo, (diff - hypo)/hypo, diff/hypo);
			//if (fabs(bond_diff - bond_hypo)>0.0001)
				printf("BOND: diff:%lf, hypo:%lf, percent: %lf, ratio:%lf\n", bond_diff, bond_hypo, (bond_diff - bond_hypo)/bond_hypo, bond_diff/bond_hypo);	
				
			CopyParticleContents(&the_membrane[i], &oldpart);
			stretch_molecule(&(the_membrane[i]), &oldpart, &(the_types[the_membrane[i].model_index]), box_length, new_box_length);  
		}
}
double get_hypothetical_bond_stretch(particle *part, mparticle *type, vector box_length, double dA)
{
	int i; 
	double area, r, duda, du; 
	vector radius; 
	du = 0;  
	area = box_length[X] * box_length[Y];  
	for (i = 0; i < type->num_bonds; i++)
	{
		vsub_periodic(radius, part->chain[type->bonds[i].b1].position, part->chain[type->bonds[i].b2].position, box_length);  
		r = sqrt(vsquare(radius)); 
		duda = r*r/area - 3*radius[Z]*radius[Z]/area; 
		duda = duda * (r - type->bonds[i].equil)/r; 
		duda = duda * type->bonds[i].cbond/2.0; 
		du = du + duda*dA; 		
	}
	 
	return du; 
}
/********************************************************************/
double get_hypothetical_da_with_bond(particle *part, mparticle *type, vector box_length, double dA)
{
	int i, b1, b2, b3; 
	double area, cosangle, r1, r2, dcostheta; 
	vector radius1, radius2; 
	double du = 0; 
	double scale1 = 3.0; 
	double scale2 = 3.0/2.0; 
	area = box_length[X]*box_length[Y]; 
	for (i = 0; i < type->num_angles; i++)
		{
			b1 = type->bondangles[i].b1; 
			b2 = type->bondangles[i].b2; 
			b3 = type->bondangles[i].b3; 
			cosangle = get_angle(*part, b1, b2, b3); 
			vsub_periodic(radius1, part->chain[b1].position, part->chain[b2].position, box_length); 		
			vsub_periodic(radius2, part->chain[b3].position, part->chain[b2].position, box_length); 
			r1 = sqrt(vsquare(radius1)); 
			r2 = sqrt(vsquare(radius2)); 
			dcostheta = -scale1 * radius1[Z]*radius2[Z]/(r1*r2); 
			dcostheta = dcostheta + scale2*radius1[Z]*radius1[Z]*cosangle/(r1*r1); 
			dcostheta = dcostheta + scale2*radius2[Z]*radius2[Z]*cosangle/(r2*r2);
			dcostheta = dcostheta/area; 
			du = du + type->bondangles[i].cbend*dcostheta; 
		} 
	du = du * dA; 
	return du; 
}

int check_bond_lengths(particle *part, mparticle *type)
{
	extern vector pbox_length; 
	int i; 
	double distance; 
	for (i = 0; i < type->num_bonds; i++)
	{
		distance = vdistance_periodic(part->chain[type->bonds[i].b1].position, part->chain[type->bonds[i].b2].position, pbox_length); 
		if (fabs(distance - type->bonds[i].equil)>0.0001)
			printf("b1:%d b2:%d distance:%lf\n", type->bonds[i].b1, type->bonds[i].b2, distance); 
	}
}
/************************************************************/

int lattice_tests(particle *the_membrane, site *the_lattice, vector box_length)
{
	check_num_nodes(the_membrane, the_lattice); 
}

int check_bead_index(particle *the_membrane, site *the_lattice)
{
	extern long pnum_particles; 
	int i, j, found; 
	cintnode *the_node; 
	
	for (i = 0; i < pnum_particles; i++)
		for (j = 0; j < the_membrane[i].chain_length; j++)
		  {
			found = NO; 
			the_node = the_lattice[the_membrane[i].chain[j].site_index].bead_list; 
			while ((the_node != NULL) && (found == NO))
			  {
				if ((the_node->contents.m == i) && (the_node->contents.b == j))
					found = YES; 
				the_node = the_node->next; 
			  }
			if (found == NO)
				printf("error!\n"); 
		  }
}

int check_num_nodes(particle *the_membrane, site *the_lattice)
{
	extern long pnum_particles;
	extern int nsites; 
	int i, total_nodes,total_beads; 
	cintnode *the_node;
	
	total_nodes = 0; 
	total_beads = 0; 
	
	for (i = 0; i < pnum_particles; i++)
		total_beads = total_beads + the_membrane[i].chain_length; 
	
	for (i = 0; i < nsites*nsites; i++)
	  {
		the_node = the_lattice[i].bead_list; 
		while (the_node!=NULL)
		  {the_node = the_node->next; total_nodes = total_nodes + 1;}
	  }
	if (total_beads != total_nodes)
		printf("total_beads: %d total_nodes: %d\n", total_beads, total_nodes); 
}

double energy_tests(particle *the_membrane, site *the_lattice, mparticle *the_types, vector box_length)
{
	double fast_energy, slow_energy, ghost_energy; 
	fast_energy = get_total_energy(the_membrane, the_lattice, the_types, box_length); 
	slow_energy = slow_total_energy(the_membrane, the_types, box_length); 
	//ghost_energy = get_total_energy_with_ghosts_slow(the_membrane, the_lattice, box_length); 
	check_range(the_membrane, the_lattice, box_length); 
	if (fabs(slow_energy-fast_energy)>0.0001)
	  {
		printf("error! fast: %lf slow: %lf\n", fast_energy, slow_energy); 
		check_range(the_membrane, the_lattice, box_length); 
	  }
	// if (fabs(slow_energy-ghost_energy)>0.0001)
	  {
	//	printf("error! ghost: %lf slow: %lf\n", ghost_energy, slow_energy); 
	  }
	return slow_energy; 
}
/*******************************************************/
int check_range(particle *the_membrane, site *the_lattice, vector box_length)
{
	extern long pnum_particles; 
	extern int num_to_check;
	int i,j,k,l,m, found, ref_site, to_find; 
	double r2,x1,x2,y1,y2,z1,z2; 
	for (i = 0; i < pnum_particles; i++)
		for (j = i + 1; j < pnum_particles; j++)
			for (k = 0; k < the_membrane[i].chain_length; k++)
				for (l = 0; l < the_membrane[i].chain_length; l++)
				  {
					r2 = vdistance_periodic(the_membrane[i].chain[k].position, the_membrane[j].chain[l].position, box_length); 
					if (r2 < 3.0*3.0)
					  {
						ref_site = the_membrane[i].chain[k].site_index; 
						to_find = the_membrane[j].chain[l].site_index; 
						found = NO; 
						for (m = 0; m < num_to_check; m++)
							if (the_lattice[ref_site].checklist[m] == to_find)
								found = YES; 
						if (found == NO)
						  {
							x1 = the_membrane[i].chain[k].position[X]; 
							x2 = the_membrane[j].chain[l].position[X]; 
							y1 = the_membrane[i].chain[k].position[Y]; 
							y2 = the_membrane[j].chain[l].position[Y]; 
							z1 = the_membrane[i].chain[k].position[Z]; 
							z2 = the_membrane[j].chain[l].position[Z]; 
							r2 = vdistance_periodic(the_membrane[i].chain[k].position, the_membrane[j].chain[l].position, box_length); 
						  }
						  }
				  }
}
/*******************************************************/
double slow_total_energy(particle *the_membrane, mparticle *the_types, vector box_length)
{
	extern long pnum_particles; 
	double bend =0; 
	double total = 0; 
	double bond = 0; 
	int i,j,k,l; 
	double hold; 
	FILE *fp; 
	fp = fopen("energymatrix","w"); 
 
	for (i = 0; i < pnum_particles; i++)
		bend = bend + get_bent(&the_membrane[i], &the_types[the_membrane[i].model_index]); 
	for (i = 0; i < pnum_particles; i++)
		bond = bond + get_bond_energy(&the_membrane[i], &the_types[the_membrane[i].model_index], box_length); 
	for (i = 0; i < pnum_particles; i++)
		{
		for (k = 0; k < the_membrane[i].chain_length; k++)
			for (l = k+1; l < the_membrane[i].chain_length; l++)
				total = total + get_bead_pair_energy(the_membrane[i].chain[k], the_membrane[i].chain[l], box_length);
		for (j = i + 1; j < pnum_particles; j++)
			for (k = 0; k < the_membrane[i].chain_length; k++)
				for (l = 0; l < the_membrane[j].chain_length; l++)
				{	
					hold = get_bead_pair_energy(the_membrane[i].chain[k], the_membrane[j].chain[l], box_length); 
					total = total + hold; 
					fprintf(fp, "%d %d %d %d %lf %lf\n", i, k, j, l, global_r, hold); 
					}
		}
fclose(fp);
 return total + bend + bond;
}

/*******************************************************/
double slow_total_energy_with_ghosts(particle *the_membrane, mparticle *the_types, vector box_length)
{
	extern long pnum_particles; 
	double bend =0; 
	double total = 0; 
	int i,j,k,l; 
	for (i = 0; i < pnum_particles; i++)
		bend = bend + get_bent(&the_membrane[i], &the_types[the_membrane[i].model_index]); 
	for (i = 0; i < pnum_particles; i++)
		for (j = i + 1; j < pnum_particles; j++)
			for (k = 0; k < the_membrane[i].chain_length; k++)
				for (l = 0; l < the_membrane[j].chain_length; l++)
					total = total + get_bead_pair_energy_with_ghosts(the_membrane[i].chain[k], the_membrane[j].chain[l], box_length); 	
 return total + bend;
}

/******************************************/
int check_all_ghosts(particle *the_membrane, site **the_lattice, vector box_length)
{
	extern long pnum_particles; 
	int i, j, k; 
	int correct = YES; 
	int old_num_ghosts, new_num_ghosts; 
	vector r; 
	double shouldbe; 
	int found; 
	int gx, gy, gz; 
	double range = 4.0; 
	
	for (i = 0; i < pnum_particles; i++)
		for (j = 0; j < the_membrane[i].chain_length; j++)
			{
			vcopy(r, the_membrane[i].chain[j].position); 
			
			if (r[X] < range) gx = -1; 
			else if (r[X]+range > box_length[X]) gx = 1; 
			else gx = 0; 
			
			if (r[Y] < range) gy = -1; 
			else if (r[Y]+range > box_length[Y]) gy = 1; 
			else gy = 0; 
			
			if (r[Z] < range) gz = -1; 
			else if (r[Z]+range > box_length[Z]) gz = 1; 
			else gz = 0; 
			
			if (gx != 0)
				{shouldbe = r[X] - gx*box_length[X]; 
				found = NO; 
				for (k = 0; k < the_membrane[i].chain[j].num_ghosts; k++ )
					if (the_membrane[i].chain[j].ghosts[k].position[X] == shouldbe)
						found = YES; }
			if (found == NO)
				{printf("particle %d bead %d xposition %lf has no ghosts\n", i,j,r[X]); 
				correct = NO; }
			}
	return correct; 
}

/******************************************/
int check_particle_ghosts(particle *the_membrane, int i, vector box_length)
{
	extern long pnum_particles; 
	int j, k; 
	int correct = YES; 
	int old_num_ghosts, new_num_ghosts; 
	vector r; 
	double shouldbe; 
	int found; 
	int gx, gy, gz; 
	double range = 4.0; 
	
		for (j = 0; j < the_membrane[i].chain_length; j++)
			{
			vcopy(r, the_membrane[i].chain[j].position); 
			
			if (r[X] < range) gx = -1; 
			else if (r[X]+range > box_length[X]) gx = 1; 
			else gx = 0; 
			
			if (r[Y] < range) gy = -1; 
			else if (r[Y]+range > box_length[Y]) gy = 1; 
			else gy = 0; 
			
			if (r[Z] < range) gz = -1; 
			else if (r[Z]+range > box_length[Z]) gz = 1; 
			else gz = 0; 
			
			if (gx != 0)
				{shouldbe = r[X] - gx*box_length[X]; 
				found = NO; 
				for (k = 0; k < the_membrane[i].chain[j].num_ghosts; k++ )
					if (the_membrane[i].chain[j].ghosts[k].position[X] == shouldbe)
						found = YES; 
			if (found == NO)
				{printf("particle %d bead %d xposition %lf has no ghosts\n", i,j,r[X]); 
				correct = NO; }
				}
			}
	return correct; 
}

/******************************************/
int check_all_ghost_pairs(particle *the_membrane, site **the_lattice, vector box_length)
{
	extern long pnum_particles; 
	int i, j, k, l, m, found; 
	int correct = YES; 
	vector r1, r2, r3; 
	double d1, d2, d3; 
	for (i = 0; i < pnum_particles; i++)
		for (j = i+1; j < pnum_particles; j++)
			for (k = 0; k < the_membrane[i].chain_length; k++)
				for (l = 0; l < the_membrane[j].chain_length; l++)
				{
					vcopy(r1, the_membrane[i].chain[k].position); 
					vcopy(r2, the_membrane[j].chain[l].position); 
					d1 = vdistance(r1, r2); 
					d2 = vdistance_periodic(r1, r2, box_length); 
					if ((d1 > 9.0 ) && (d2 < 9.0))
					{
						found = NO; 					
						for (m = 0; m < the_membrane[i].chain[k].num_ghosts; m++)
						{
							vcopy(r3, the_membrane[i].chain[k].ghosts[m].position); 
							d3 = vdistance(r2, r3); 
							if (d3 < 9.0)
							{
								if (fabs(d3 - d2) > 0.0001)
									printf("%d %d %d %d Ghost within range but not correct distance\n", i,k,j,l); 
								found = YES; 			
							}
						}
						for (m = 0; m < the_membrane[j].chain[l].num_ghosts; m++)
						{
							vcopy(r3, the_membrane[j].chain[l].ghosts[m].position); 
							d3 = vdistance(r1, r3); 
							if (d3 < 9.0)
							{
								if (fabs(d3 - d2) > 0.0001)
									printf("%d %d %d %d Ghost within range but not correct distance\n", i,k,j,l); 
								found = YES; 			
							}
						}
						if (found == NO)
						{	printf("%d %d %d %d No ghost within range\n", i, k, j,l); 
							printf("bead %d %d position:\n", i, k); 
							vprint(the_membrane[i].chain[k].position); 
							printf("bead %d %d position:\n", j,l); 
							vprint(the_membrane[j].chain[l].position); 
							correct = NO; 
						}
					}
				}
return correct; 
				}
/******************************************************/
double check_molecule_energy(particle *the_membrane, site *the_lattice, mparticle *the_types,  int index, vector box_length)
{
	extern long pnum_particles; 
	int i, j, k; 
	double hold, fast_bead_energy, slow_bead_energy; 
	double total = 0; 
	double total2 = 0; 
	for (i = 0; i < the_membrane[index].chain_length; i++)
		{
			slow_bead_energy = 0; 
			for (j = 0; j < pnum_particles; j++)
				for (k = 0; k < the_membrane[j].chain_length; k++)
					{hold = 0.5*get_bead_pair_energy( the_membrane[index].chain[i], the_membrane[j].chain[k], box_length); 
					slow_bead_energy = slow_bead_energy + hold; 
					if (fabs(hold)>0.0000001)
						printf("i:%d, k:%d energy:%lf\n", i, k, hold); 
					total = total + hold; 
					}
			fast_bead_energy = get_bead_energy(the_membrane[index].chain[i], the_membrane, the_lattice, box_length); 
			total2 = total2 + fast_bead_energy;
			//printf("i:%d, fast_bead_energy:%lf\n", i, fast_bead_energy);  
			if (fabs(fast_bead_energy - slow_bead_energy) > 0.00001)
				printf("Molecule %d bead %d fast_energy: %lf slow_energy:%lf\n", index, i, fast_bead_energy, slow_bead_energy); 
		}
	total = total + get_bent(&(the_membrane[index]), &the_types[the_membrane[index].model_index]); 
	return total2; 
	
}