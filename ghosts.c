#include <stdlib.h>
#include <math.h>
#include "structures.h"
#include "vectorops.h"
#include "ghosts.h"

int assign_all_ghosts(particle *the_membrane, vector box_length)
{
	int i, j; 
	extern long pnum_particles; 
	for (i = 0; i < pnum_particles; i++)
		for (j = 0; j < the_membrane[i].chain_length; j++)
			assign_ghosts(the_membrane[i].chain[j], box_length); 
			
}

int assign_ghosts(bead the_bead, vector box_length)
{
	double range = 4.0; 
	vector r; 
	int gx, gy, gz; 
	int count = 0; 
	
	vcopy(r, the_bead.position); 
	
	if (r[X] < range) gx = -1; 
	else if (r[X]+range > box_length[X]) gx = 1; 
	else gx = 0; 
		
	if (r[Y] < range) gy = -1; 
	else if (r[Y]+range > box_length[Y]) gy = 1; 
	else gy = 0; 
	
	if (r[Z] < range) gz = -1; 
	else if (r[Z]+range > box_length[Z]) gz = 1; 
	else gz = 0; 
	
	if ((gx == 0) && ((gy == 0) && (gz == 0))) return 0;
	
	if (gx != 0)
	{
		bead_copy(&(the_bead.ghosts[count]), &the_bead);
		the_bead.ghosts[count].position[X] = r[X] - gx * box_length[X]; 
		count++; 
		if (gy != 0)
		{
			bead_copy(&(the_bead.ghosts[count]), &the_bead);
			the_bead.ghosts[count].position[X] = r[X] - gx * box_length[X]; 
			the_bead.ghosts[count].position[Y] = r[Y] - gy * box_length[Y]; 
			count++; 
			if (gz != 0) 
			{
				bead_copy(&(the_bead.ghosts[count]), &the_bead);
				the_bead.ghosts[count].position[X] = r[X] - gx * box_length[X]; 
				the_bead.ghosts[count].position[Y] = r[Y] - gy * box_length[Y]; 
				the_bead.ghosts[count].position[Z] = r[Z] - gz * box_length[Z]; 
				count++; 
			}
		}
		if (gz !=0)
		{
				bead_copy(&(the_bead.ghosts[count]), &the_bead);
				the_bead.ghosts[count].position[X] = r[X] - gx * box_length[X]; 
				the_bead.ghosts[count].position[Z] = r[Z] - gz * box_length[Z]; 
				count++; 
		}
	}
	if (gy != 0)
		{
			bead_copy(&(the_bead.ghosts[count]), &the_bead);
			the_bead.ghosts[count].position[Y] = r[Y] - gy * box_length[Y]; 
			count++; 
			if (gz != 0) 
			{
				bead_copy(&(the_bead.ghosts[count]), &the_bead);
				the_bead.ghosts[count].position[Y] = r[Y] - gy * box_length[Y]; 
				the_bead.ghosts[count].position[Z] = r[Z] - gz * box_length[Z]; 
				count++; 
			}
		}
		if (gz !=0)
		{
				bead_copy(&(the_bead.ghosts[count]), &the_bead); 
				the_bead.ghosts[count].position[Z] = r[Z] - gz * box_length[Z]; 
				count++; 
		}
		
	return count; 
}

int bead_copy(bead *hold, bead *copy)
{
	vcopy(hold->position, copy->position); 
	hold->ci.m = copy->ci.m; 
	hold->ci.b = copy->ci.b; 
	hold->site_index = copy->site_index; 
	hold->num_ghosts = 0;  
	hold->user = copy->user; 
}
