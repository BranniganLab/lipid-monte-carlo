/*8/28/02- comments finished, tentative release version.  Some overflow checks left in.*/

#include <stdlib.h>

#include <stdio.h>
#include <math.h>
#include "structures.h"
#include "in.h"
#include "init.h"
#include "lattice.h"
#include "energy.h"
#include "out.h"
#include "vectorops.h"
#include "domain.h"
#include "ghosts.h"
#include "ghost_energy.h"
#include "fourier.h"

void init(particle *the_membrane, site *the_lattice, mparticle *the_types, double *total_energy, int total_accept[], double move_size[], long step)
{
    extern vector pbox_length;
    extern char pinfile[30];
	extern vector last_update;
 	extern int pnumber_user_interactions;
  extern upotential *puser_potentials;
	extern long pnum_particles; 
	FILE *fp;

 	int i,j;
    /*read the data with the xyz coordinates*/
    printf("reading file...\n");
    if ((fp = fopen(pinfile,"r"))==NULL)
    {printf("cannot open infile: %s\n",pinfile);
        exit(0);}
    else
    {
    read_top(the_types);
//    reorder_bonds(the_types); 

    get_one_ensemble(fp, the_membrane, the_types);
    set_chainpos(the_types); 
        fclose(fp); }
	//for (i = 0; i < pnum_particles; i++)
	//	for (j = 0; j < the_membrane[i].chain_length; j++)
	//	{
	//		the_membrane[i].chain[j].ghosts = (bead *)(malloc(MAX_GHOSTS * sizeof(bead))); 
	//		the_membrane[i].chain[j].num_ghosts = assign_ghosts(the_membrane[i].chain[j], pbox_length); 
	//	}
	enforce_constraints(the_membrane, pbox_length); 
	printf("enforced constraints...\n"); 
	get_checklists(the_lattice, pbox_length); 
	printf("got checklists....\n"); 
	vcopy(last_update, pbox_length); 
	initialize_lattice(the_membrane, the_lattice, pbox_length); 
	printf("initialized lattice\n"); 
	initialize_domains( the_membrane, the_lattice, the_types, pbox_length);
	printf("initialized domains\n"); 
    for (i =0 ; i < pnumber_user_interactions; i++)
    		{
      		the_membrane[puser_potentials[i].ci1.m].chain[puser_potentials[i].ci1.b].user = YES;
       		 the_membrane[puser_potentials[i].ci2.m].chain[puser_potentials[i].ci2.b].user = YES;
          }
	*total_energy =  get_total_energy(the_membrane, the_lattice, the_types, pbox_length);
    printf("Total energy: %12.6f\n", *total_energy) ;
	
    /*write out once*/
    write_interval( the_membrane, the_types, *total_energy,total_accept,
                    move_size,step, pbox_length) ;

}
/**************************************************************************/
void init_potential(potential U)
{
    int i;

    U.num_terms = 0;
    U.range = 0.0;
   	for (i = 0; i<5; i++)
        {U.coefficients[i][0] = 0.0;
            U.coefficients[i][1]=0.0;   }

}


void free_membrane(particle *the_membrane)
{
    extern long pnum_particles;
    int i,j;
    for (i = 0; i< pnum_particles; i++)
        for (j = 0; j < the_membrane[i].chain_length; j++)
            if (the_membrane[i].chain[j].ntree != NULL)
                free(the_membrane[i].chain[j].ntree);
    free(the_membrane);
}

int set_chainpos(mparticle *the_types)
{
	extern long pnum_types;
	int i, j;
	for (j = 0; j < pnum_types; j++)
	{ 
	if (the_types[j].num_chains == 1)
		for (i = 0; i < the_types[j].chain_length; i++)
			the_types[j].chain[i].chainpos = i; 
	else if (the_types[j].num_chains ==2)
		{
			for (i = 0; i < the_types[j].chain1start; i++)
				the_types[j].chain[i].chainpos = -i; 	
			for (i = the_types[j].chain1start; i < the_types[j].chain2start; i++)
				the_types[j].chain[i].chainpos = i; 	
			for (i = the_types[j].chain2start; i < the_types[j].chain_length; i++)
				the_types[j].chain[i].chainpos = i - (the_types[j].chain2start - the_types[j].chain1start) + 1;
		}
	}
}

/*int reorder_bonds(mparticle *the_types)
{
	extern long pnum_types; 
	int i; 
	for (i = 0; i < pnum_types; i++)
	{
		for (j = 0; j < the_types[i].num_bonds; j++)
		{
			if (the_types[i].bonds[j].b1 > the_types[i].bonds[j].b2)
				{
					temp = the_types[i].bonds[j].b2; 
					the_types[i].bonds[j].b2 = the_types[i].bonds[j].b1; 
					the_types[i].bonds[j].b1 = temp; 
				}
			neworder[j] = j; 
		}
		for (j = 0; j < the_types[i].num_bonds; j++)
		{
			if (the_types[i].bonds[j].b1 < the_types[i].ibindex)
				{
					for (k = 0; k < the_types[i].num_bonds; k++)
						if (the_types[i].bonds[j].b1 > the_types[i].bonds[neworder[k]].b1;
							{	
				}
				
		}
	}
}*/