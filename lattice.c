#include <stdlib.h>
#include "structures.h"
#include "lattice.h"
#include "vectorops.h"
/****************/
int check_checklists(site *the_lattice, vector old_box_length, vector new_box_length)
{
	extern int nsites; 
	vector old_spacing, new_spacing, old_range, new_range; 
	old_spacing[X] = old_box_length[X]/((double)nsites); 
	old_spacing[Y] = old_box_length[Y]/((double)nsites); 
	new_spacing[X] = new_box_length[X]/((double)nsites); 
	new_spacing[Y] = new_box_length[Y]/((double)nsites); 
	
	old_range[X] = (int)((3.0 + old_spacing[X] + 0.5)/old_spacing[X]); 
	old_range[Y] = (int)((3.0 + old_spacing[X] + 0.5)/old_spacing[Y]);
	
	new_range[X] = old_range[X] * new_spacing[X]; 
	new_range[Y] = old_range[Y] * new_spacing[Y]; 
	
	if (new_range[X] < (3.0 + old_spacing[X] + 0.0))
		return YES; 
	if (new_range[Y] < (3.0 + old_spacing[Y] + 0.0))
		return YES; 
	if (new_range[X] > (4.0 + old_spacing[X] + 0.0))
		return YES; 
	if (new_range[Y] > (4.0 + old_spacing[Y] + 0.0))
		return YES; 
	return NO; 
}
/****************/
int get_checklists(site *the_lattice, vector box_length)
{
	extern int nsites; 
	extern int num_to_check; 
	vector spacing, v1, v2; 
	double range, radius; 
	int count, i, j, ref_site_x, ref_site_y; 
	
	spacing[X] = box_length[X]/((double)nsites); 
	spacing[Y] = box_length[Y]/((double)nsites); 
	range = 3.0 + spacing[X] + 0.5; 
	for (ref_site_x = 0; ref_site_x < nsites; ref_site_x++)
		for (ref_site_y = 0; ref_site_y < nsites; ref_site_y++)
		  {
			//printf("x:%d, y:%d\n", ref_site_x,ref_site_y); 
			count = 0; 
			for (i = 0; i < nsites; i++)
				for (j = 0; j < nsites; j++)
				  {
					v1[X] = (ref_site_x+0.5) * spacing[X]; v2[X] = (i+0.5)*spacing[X]; 
					v1[Y] = (ref_site_y+0.5) * spacing[Y]; v2[Y] = (j+0.5)*spacing[Y]; 
					v1[Z] = 0; v2[Z] = 0; 
					radius = vdistance_periodic(v1,v2, box_length); 
					if (radius < range*range)
					  {
						the_lattice[ref_site_x+nsites*ref_site_y].checklist[count] = i + nsites * j; 
						count++; 
					  }
				  }
		  }
			num_to_check = count; 
			return count; 
				
}

/****************/
int initialize_lattice(particle *the_membrane, site *the_lattice, vector box_length)
{
	extern int nsites; 
	extern long pnum_particles; 
	int i, j, xindex, yindex, site_index; 
	double xspacing, yspacing; 
	vector r; 
	chainindex ci; 
	int temp; 
	temp = nsites; 
	
	xspacing = box_length[X]/nsites; 
	yspacing = box_length[Y]/nsites; 
	
	for (i = 0; i < pnum_particles; i++)
		for (j = 0; j < the_membrane[i].chain_length; j++)
		  {
			vcopy(r, the_membrane[i].chain[j].position); 
			xindex = (int)(r[X]/xspacing);
			yindex = (int)(r[Y]/yspacing); 
			if (xindex < 0) xindex = nsites + xindex; 
			else if (xindex>=nsites) xindex = xindex - nsites; 
			if (yindex < 0) yindex = nsites + yindex; 
			else if (yindex>=nsites) yindex = yindex - nsites;
			site_index = xindex + yindex*nsites;
			if (fabs(site_index) > (nsites*nsites))
				{printf("Error: site_index: %d, xindex:%d yindex:%d, nsites*nsites: %d\n", site_index, xindex, yindex,  nsites*nsites); 
				printf("xspacing: %lf yspacing: %lf\n", xspacing, yspacing); 
				vprint(r); }
			ci.m = i; ci.b = j; 
			add_bead_to_lattice(ci, &(the_lattice[site_index])); 
			the_membrane[i].chain[j].site_index = site_index; 
		  }
			get_checklists(the_lattice, box_length); 
}

/****************/
int add_bead_to_lattice(chainindex ci,site *new_site)
{
	cintnode *tmp; 
	tmp = new_site->bead_list; 
	new_site->bead_list = (cintnode *)(malloc(sizeof(cintnode)));  
	new_site->bead_list->contents.m = ci.m; 
	new_site->bead_list->contents.b = ci.b; 
	new_site->bead_list->next = tmp;
	if (tmp != NULL) new_site->bead_list->next->previous = new_site->bead_list; 
	new_site->bead_list->previous = NULL; 
}

/****************/
int update_beads_lattice_site(bead *the_bead, site *the_lattice, vector box_length)
{
	extern int nsites;
	vector r; 
	double xspacing, yspacing; 
	int xindex, yindex, site_index; 
	vcopy(r, the_bead->position); 
	xspacing = box_length[X]/nsites; 
	yspacing = box_length[Y]/nsites; 
	xindex = (int)(r[X]/xspacing);
	yindex = (int)(r[Y]/yspacing); 
	if (xindex < 0) xindex = nsites + xindex; 
	else if (xindex>=nsites) xindex = xindex - nsites; 
	if (yindex < 0) yindex = nsites + yindex; 
	else if (yindex>=nsites) yindex = yindex - nsites; 
	site_index = xindex + yindex*nsites; 
	if (site_index == the_bead->site_index)
		return site_index; 
	move_bead_in_lattice(the_bead->ci, &(the_lattice[the_bead->site_index]),&(the_lattice[site_index])); 
	the_bead->site_index = site_index; 
	
}

/****************/
int move_bead_in_lattice(chainindex ci,site *old_site, site *new_site)
{
	cintnode *tmp, *tmp2; 
	extern int nsites; 
	int found = NO; 
	tmp = old_site->bead_list; 
	while ((found == NO) && (tmp != NULL))
	  {if ((tmp->contents.m == ci.m) && (tmp->contents.b == ci.b))
		  found = YES; 
		else tmp = tmp->next; }
	if (tmp == NULL) 
		return -1; 
	if (tmp->previous != NULL) 
	  {
		tmp->previous->next = tmp->next; }
	else 
		old_site->bead_list = tmp->next; 
	if (tmp->next !=NULL)
		tmp->next->previous = tmp->previous; 
	if (tmp->next!= NULL) 
		tmp->next->previous = tmp->previous; 
	tmp2 = new_site->bead_list; 
	new_site->bead_list = tmp; 
	if (tmp2 == NULL)
	  {
		tmp->next = NULL;}
	else
	  { tmp->next = tmp2; 
		tmp2->previous = tmp;} 
	tmp->previous = NULL;
	
}

/*********************************/
void free_node_list(cintnode *head)
{
    cintnode *nextnode, *thisnode;
    thisnode = head;
    while (thisnode != NULL)
	  {
        nextnode = thisnode->next;
        free(thisnode);
        thisnode = nextnode;
	  }
}
