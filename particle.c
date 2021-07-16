#include <stdlib.h>
#include "structures.h"
#include "particle.h"
#include "perturb.h"
#include "step.h"
#include "accept.h"
#include "energy.h"
#include "model.h"
#include "lattice.h"
#include "domain.h"
#include "ghosts.h"
#include "ghost_energy.h"
#include "debug.h"
#include "cluster.h"
#include "preset.h"
double particle_move(particle *the_membrane,  site *the_lattice, mparticle *the_types, int attempts[],int accepts[])
{
   extern int pdo_run_preset, pgenerate_preset, iteration; 
   extern step *preset_run; 

    extern int abort_move, pdomain_fixed;
    extern vector pbox_length;
    extern double pmove_size[NUM_DOF], total_system_energy;
    int pert_index, j, move, accept;
    int move_class[NUM_DOF] ={X,Y,ROTATION,ROTATION,ROTATION,VOLUME,VOLUME + Y, CLUSTER};
    double  energy_diff, temp_particle_energy,temp_double;
	if (pgenerate_preset == YES)
		move_class[1] = ROTATION; 
    abort_move = NO_OVERLAP;
    energy_diff = 0;
    /*find out which particle to perturb*/
    pert_index = GetPertIndex();
	if (the_membrane[pert_index].domain_index != 0)
		{
  			if (the_membrane[pert_index-1].domain_index == the_membrane[pert_index].domain_index )
  			//if (1 == 0)
				return 0.0;
			else if (the_membrane[pert_index].domain_index == pdomain_fixed)
				return 0.0;
		    else if (the_membrane[pert_index].domainptr != NULL)
			{
				int domain_move_type; 
				extern long seed;
    			long *base= &seed; 
				if (ran3(base) < 0)
					domain_move_type = 0; 
				else
					domain_move_type = 1; 
				if (pdo_run_preset == YES)
					domain_move_type = preset_run[iteration].domain_move_type; 
				if (pgenerate_preset == YES)
					preset_run[iteration].domain_move_type = domain_move_type; 

				if (domain_move_type == 0)
					return cluster_move(the_membrane[pert_index].domainptr, the_membrane, the_lattice, the_types, the_membrane[pert_index].domain_index, pbox_length);  
				else if (domain_move_type == 1)
					return step_domain(the_membrane[pert_index].domainptr, the_membrane, the_lattice, the_types, pmove_size[X], pbox_length);
						   if ((pdo_run_preset == YES) || (pgenerate_preset == YES))
	     output_preset_results(preset_run[iteration]);
	  		iteration = iteration + 1; 
			}
			else
     			{printf("Error! Non zero domain index has null domain pointer! \n"); return 0.0; }
       }
	the_membrane[pert_index].energy =get_molecule_energy(pert_index, the_membrane, the_lattice, the_types, pbox_length);
	
    /*get order of perturbations*/
    ShuffleMoveClass(move_class);
	
    /*loop through the DOF's to perturb*/
    for (j = 0; j <VOLUME; j++)
	  {
        abort_move = NO_OVERLAP;
        move = move_class[j];
	if (pdo_run_preset == YES)
    		move = preset_run[iteration].regular_move_type; 
	if (pgenerate_preset == YES)
			{preset_run[iteration].regular_move_type = move; 
			preset_run[iteration].m = pert_index;}
       // printf("move: %d, preset_run[iteration].regular_move_type: %d\n", move, preset_run[iteration].regular_move_type); 
        if ((pmove_size[move] > 0.0)){
            /*generate and perturb the particle, calculate the new energy,
			and decide acceptance*/


				if (move == X)
					energy_diff = energy_diff + step_dof(pert_index, the_membrane, the_lattice, the_types, &accept, pmove_size, pbox_length);
				
				else if (move == Y)
					{energy_diff = energy_diff + step_bond(pert_index, the_membrane, the_lattice, the_types,  &accept, pmove_size,  pbox_length);
					}
				else
					energy_diff = energy_diff + step_bead(pert_index, the_membrane, the_lattice, the_types, &accept, pmove_size, pbox_length);
				accepts[move] =accepts[move] + accept;
				attempts[move] = attempts[move] + 1;
				total_system_energy = total_system_energy + energy_diff; 
	   if ((pdo_run_preset == YES) || (pgenerate_preset == YES))
	     output_preset_results(preset_run[iteration]);
		iteration = iteration + 1; 
		
		}
	  }
	
	
	return energy_diff;
}


/**********************************************************************************/
double step_dof(int pert_index, particle *the_membrane, site *the_lattice, mparticle *the_types, int *accept, double *delta, vector box_length)
{
	extern vector pbox_length;
    particle *temp_part;
    double old_particle_energy, new_particle_energy, energy_diff,temp_double, right_molecule_energy;
	int bond,j, constraint_index;
	extern int pnumber_user_constraints; 
	extern uconstraint *puser_constraints; 
    chainindex ci; 
    vector translation, rotation;
//	    right_molecule_energy = check_molecule_energy(the_membrane, the_lattice,  pert_index, box_length); 
    //if (fabs(right_molecule_energy - the_membrane[pert_index].energy)>0.00001)
    //	printf("ERROR in stepdof beginning! right molecule %d energy:%lf , listed molecule %d energy: %lf\n", pert_index, right_molecule_energy, pert_index, the_membrane[pert_index].energy); 
	*accept = 0; 
	temp_part = (particle *)malloc(sizeof(particle));
    CopyParticleContents(&(the_membrane[pert_index]), temp_part);
	old_particle_energy = (the_membrane[pert_index]).energy;
	new_particle_energy = get_bent(&(the_membrane[pert_index]),&(the_types[the_membrane[pert_index].model_index])) + get_bond_energy(&(the_membrane[pert_index]),&(the_types[the_membrane[pert_index].model_index]), box_length); 
	ci.m = pert_index; ci.b = 0; 
	if (perturb_particle(&(the_membrane[pert_index]), &(the_types[the_membrane[pert_index].model_index]), ci, translation, delta[X], box_length, PERTURB_WHOLE_CHAIN) == NO)
	{free(temp_part); return 0; }
		
	for (j = 0; j < the_membrane[pert_index].chain_length; j++)
	{
		update_beads_lattice_site(&(the_membrane[pert_index].chain[j]), the_lattice, box_length);
	//	the_membrane[pert_index].chain[j].num_ghosts = assign_ghosts(the_membrane[pert_index].chain[j], box_length); 
	}
	for (j = 0; j < the_membrane[pert_index].chain_length; j++)
		new_particle_energy = new_particle_energy + get_bead_energy(the_membrane[pert_index].chain[j], the_membrane, the_lattice, box_length);

    energy_diff = 2 *(new_particle_energy - old_particle_energy);
	*accept = test_acceptance(energy_diff);
	
    if (*accept ==0)
	  {/*if there's a rejection */
		
        /*there will be no energy change in the move because no move will occur.    */
        energy_diff = 0;
        /*return perturbed particle to its original state.*/
        CopyParticleContents(temp_part, &(the_membrane[pert_index]));
        the_membrane[pert_index].energy = old_particle_energy;
		for (j = 0; j < the_membrane[pert_index].chain_length; j++)
		{
			update_beads_lattice_site(&(the_membrane[pert_index].chain[j]), the_lattice, box_length);
	 	//	the_membrane[pert_index].chain[j].num_ghosts = assign_ghosts(the_membrane[pert_index].chain[j], box_length); 
		}/*end if there's a rejection. 		*/
	}
    else if (*accept == 1)
        (the_membrane[pert_index]).energy = old_particle_energy + energy_diff/2.0;
   // right_molecule_energy = check_molecule_energy(the_membrane, the_lattice, pert_index, box_length); 
   // if (fabs(right_molecule_energy - the_membrane[pert_index].energy)>0.00001)
   // 	printf("ERROR in stepdof! right molecule %d energy:%lf , listed molecule %d energy: %lf, accept:%d\n", pert_index, right_molecule_energy, pert_index, the_membrane[pert_index].energy, *accept);  
    //printf("******checking ghosts in step_dof\n*******"); 
	//check_all_ghosts(the_membrane, the_lattice, box_length); 
	free(temp_part);
    return energy_diff;
}

/**********************************************************************************/
double step_bead(int pert_index, particle *the_membrane, site *the_lattice, mparticle *the_types,  int *accept, double *delta,  vector box_length)
{
	extern int pdo_run_preset, pgenerate_preset, iteration; 
   extern step *preset_run; 
    extern long inum_lipids;
    extern long seed; 
    extern int abort_move, pnum_cells;
    extern vector pbox_length;
    particle *temp_part;
    double old_bead_energy, new_bead_energy, energy_diff,temp_double;
    double old_bend_energy, new_bend_energy, old_particle_energy; 
    double old_self_energy, new_self_energy, right_molecule_energy; 
    double old_bond_energy, new_bond_energy; 
    int recheck;
    int bond, bead_index, nbeads;
    vector temp_pos, translation, rotation;
    long *base = &seed;
    chainindex ci;
	extern int pnumber_user_constraints; 
	extern uconstraint *puser_constraints; 
	int constraint_index,i; 
	*accept = 0; 
	//    right_molecule_energy = check_molecule_energy(the_membrane, the_lattice, pert_index, box_length); 
   // if (fabs(right_molecule_energy - the_membrane[pert_index].energy)>0.00001)
  //  	printf("stepbead beginning! right molecule %d energy:%lf , listed molecule %d energy: %lf\n", pert_index, right_molecule_energy, pert_index, the_membrane[pert_index].energy); 
    nbeads = (the_membrane[pert_index]).chain_length;
    bead_index = (int) (nbeads* (ran3(base)));
    if (pdo_run_preset == YES)
    	bead_index = preset_run[iteration].b; 
    if (pgenerate_preset == YES)
		preset_run[iteration].b = bead_index; 
    if (bead_index < 0)  bead_index = 0;
    if (bead_index >= nbeads) bead_index = bead_index - 1;

	ci.m = pert_index; 
	ci.b = bead_index; 
	temp_part = (particle *)malloc(sizeof(particle));
    CopyParticleContents(&(the_membrane[pert_index]), temp_part);

    old_particle_energy = (the_membrane[pert_index]).energy;
    ci.m = pert_index; ci.b = bead_index;
	old_bead_energy = get_bead_energy(the_membrane[pert_index].chain[bead_index], the_membrane, the_lattice, box_length);
	old_bend_energy = get_bent(&(the_membrane[pert_index]), &(the_types[the_membrane[pert_index].model_index]));
	old_self_energy = get_molecule_self_energy(the_membrane[pert_index], box_length); 
	old_bond_energy = get_bond_energy(&(the_membrane[pert_index]), &(the_types[the_membrane[pert_index].model_index]), box_length); 
	if (perturb_particle(&(the_membrane[pert_index]), &(the_types[the_membrane[pert_index].model_index]),ci, translation, delta[ROTATION], box_length, PERTURB_BEAD) == NO)
	{free(temp_part); return 0; }	

	update_beads_lattice_site(&(the_membrane[pert_index].chain[bead_index]), the_lattice, box_length);
	//for (i = 0; i < the_membrane[pert_index].chain_length; i++)
	//	the_membrane[pert_index].chain[i].num_ghosts = assign_ghosts(the_membrane[pert_index].chain[i], box_length); 
	new_bead_energy = get_bead_energy(the_membrane[pert_index].chain[bead_index], the_membrane, the_lattice, box_length);
    new_bend_energy = get_bent(&(the_membrane[pert_index]),&(the_types[the_membrane[pert_index].model_index]));
    new_self_energy = get_molecule_self_energy(the_membrane[pert_index], box_length); 
    new_bond_energy = get_bond_energy(&(the_membrane[pert_index]), &(the_types[the_membrane[pert_index].model_index]), box_length);
    	if (fabs(new_bond_energy - old_bond_energy)>0.0001)
		printf("new:%lf old%lf\n", new_bond_energy, old_bond_energy); 
    energy_diff = 2 *(new_bead_energy - old_bead_energy)  + (new_bend_energy - old_bend_energy);
	
	*accept = test_acceptance(energy_diff);
	//printf("energy diff: %lf, accept:%d\n", energy_diff, *accept); 
	if (*accept ==0)
	  {/*if there's a rejection */
        energy_diff = 0;
        /*return perturbed particle to its original state.*/
        CopyParticleContents(temp_part, &(the_membrane[pert_index]));
        the_membrane[pert_index].energy = old_particle_energy;
		update_beads_lattice_site(&(the_membrane[pert_index].chain[bead_index]), the_lattice, box_length);
		//for (i = 0; i < the_membrane[pert_index].chain_length; i++)
		//	the_membrane[pert_index].chain[i].num_ghosts = assign_ghosts(the_membrane[pert_index].chain[i], box_length); 
	  }/*end if there's a rejection. 		*/
	
    else if (*accept == 1)
	  {
		ran3(base);
        (the_membrane[pert_index]).energy = old_particle_energy + new_bend_energy + new_bead_energy - old_bead_energy - old_bend_energy + (new_self_energy - old_self_energy)/2.0;
	  }
	 
	//printf("******checking ghosts in step_bead\n*******"); 
	//if (check_all_ghosts(the_membrane, the_lattice, box_length)== NO)
//	if (check_all_ghosts(the_membrane, the_lattice, box_length) == NO)
//		{
//		printf("accept: %d, pert_index: %d, bead_index: %d\n", *accept, pert_index, bead_index); 
		//the_membrane[pert_index].chain[bead_index].num_ghosts = assign_ghosts(the_membrane[pert_index].chain[bead_index], box_length); 
//		}
	//check_bead_index(the_membrane,the_lattice);
 //   right_molecule_energy = check_molecule_energy(the_membrane, the_lattice, pert_index, box_length); 
  //  if (fabs(right_molecule_energy - the_membrane[pert_index].energy)>0.00001)
   // 	printf("ERROR in stepbead %d! right molecule %d energy:%lf , listed molecule %d energy: %lf, accept:%d, old_bead_energy:%lf, new_bead_energy:%lf\n", bead_index,pert_index, right_molecule_energy, pert_index, the_membrane[pert_index].energy, *accept, old_bead_energy, new_bead_energy); 
    free(temp_part);
    return energy_diff;
}
/************************************************/
/**********************************************************************************/
double step_bond(int pert_index, particle *the_membrane, site *the_lattice, mparticle *the_types,  int *accept, double *delta,  vector box_length)
{
	extern int pdo_run_preset, pgenerate_preset, iteration; 
   extern step *preset_run; 
    extern long inum_lipids;
    extern long seed; 
    extern int abort_move, pnum_cells;
    extern vector pbox_length;
    particle *temp_part;
    double old_bead_energy, new_bead_energy, energy_diff,temp_double;
    double old_bend_energy, new_bend_energy, old_particle_energy; 
    double old_self_energy, new_self_energy, right_molecule_energy; 
    double old_bond_energy, new_bond_energy; 
    int recheck;
    int bond, bond_index, nbonds;
    int b1, b2; 
    int mi; 
    vector temp_pos, translation, rotation;
    long *base = &seed;
    chainindex ci;
	extern int pnumber_user_constraints; 
	extern uconstraint *puser_constraints; 
	int constraint_index,i; 
	*accept = 0; 
    mi = the_membrane[pert_index].model_index; 
    nbonds = the_types[mi].num_bonds; 
    bond_index = (int) (nbonds* (ran3(base)));
    if (pdo_run_preset == YES)
    	bond_index = preset_run[iteration].bond_index; 
	if (pgenerate_preset == YES)
		preset_run[iteration].bond_index = bond_index; 
    b1 = the_types[mi].bonds[bond_index].b1; 
    b2 = the_types[mi].bonds[bond_index].b2; 
    if (bond_index < 0)  bond_index = 0;
    if (bond_index >= nbonds) bond_index = bond_index - 1;

	ci.m = pert_index; 
	ci.b = bond_index; 
	temp_part = (particle *)malloc(sizeof(particle));
    CopyParticleContents(&(the_membrane[pert_index]), temp_part);

    old_particle_energy = (the_membrane[pert_index]).energy;
	old_bead_energy = get_bead_energy(the_membrane[pert_index].chain[b1], the_membrane, the_lattice, box_length);
	old_bead_energy = old_bead_energy + get_bead_energy(the_membrane[pert_index].chain[b2], the_membrane, the_lattice, box_length);
	old_bend_energy = get_bent(&(the_membrane[pert_index]), &(the_types[the_membrane[pert_index].model_index]));
	old_self_energy = get_molecule_self_energy(the_membrane[pert_index], box_length); 
	old_bond_energy = get_bond_energy(&(the_membrane[pert_index]), &(the_types[the_membrane[pert_index].model_index]), box_length);
	if (perturb_particle(&(the_membrane[pert_index]), &(the_types[the_membrane[pert_index].model_index]),ci, translation, delta[Y], box_length, PERTURB_BOND) == NO)
	{free(temp_part); return 0; }	

	update_beads_lattice_site(&(the_membrane[pert_index].chain[b1]), the_lattice, box_length);
	update_beads_lattice_site(&(the_membrane[pert_index].chain[b2]), the_lattice, box_length);
	new_bead_energy = get_bead_energy(the_membrane[pert_index].chain[b1], the_membrane, the_lattice, box_length);
	new_bead_energy = new_bead_energy + get_bead_energy(the_membrane[pert_index].chain[b2], the_membrane, the_lattice, box_length);
    new_bend_energy = get_bent(&(the_membrane[pert_index]),&(the_types[the_membrane[pert_index].model_index]));
    new_self_energy = get_molecule_self_energy(the_membrane[pert_index], box_length); 
    new_bond_energy = get_bond_energy(&(the_membrane[pert_index]), &(the_types[the_membrane[pert_index].model_index]), box_length);
    energy_diff = 2 *(new_bead_energy - old_bead_energy)  + (new_bend_energy - old_bend_energy) + (new_bond_energy - old_bond_energy); 
	
	*accept = test_acceptance(energy_diff);
	
	if (*accept ==0)
	  {/*if there's a rejection */
        energy_diff = 0;
        /*return perturbed particle to its original state.*/
        CopyParticleContents(temp_part, &(the_membrane[pert_index]));
        the_membrane[pert_index].energy = old_particle_energy;
		update_beads_lattice_site(&(the_membrane[pert_index].chain[b1]), the_lattice, box_length);
		update_beads_lattice_site(&(the_membrane[pert_index].chain[b2]), the_lattice, box_length);
	  }/*end if there's a rejection. 		*/
	
    else if (*accept == 1)
	  {
	
        (the_membrane[pert_index]).energy = old_particle_energy + new_bend_energy + new_bond_energy - old_bond_energy + new_bead_energy - old_bead_energy - old_bend_energy + (new_self_energy - old_self_energy)/2.0;
	  }
    free(temp_part);
    return energy_diff;
}
