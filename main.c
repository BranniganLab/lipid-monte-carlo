#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "structures.h"
#include "pref.h"
#include "init.h"
#include "out.h"
#include "particle.h"
#include "volume.h"
#include "energy.h"
#include "debug.h"
#include "fourier.h"
#include "preset.h"
long seed;

/*long psample_interval = 1;

int new_clusters = 0;

int aborted = 0;
int stop_run=0;
double dT = 0.0;
average cluster_accepts, cluster_rejects;
int **checklist;
double dUdA[2];
int pnum_cells =NUM_CELLS;
double *lengths;
double *initlengths;
vector total_pressure;
vector old_pressure;
double plong_range = 0;*/
double pfourier_moves; 
double total_system_energy; 
double percent_change; 
double percent_change_samples; 
int stop_run = 0; 
int abort_move = 0; 
int samples;
vector last_update; 
int nsites =5;
int num_to_check; 
int iteration; 
step *preset_run;
int main (int argc, const char * argv[]) {
	get_sim_prefs();
    get_model_prefs();
    get_domain_prefs();
	get_user_prefs();
	RunSimulation();
    return 0;
}

/************************************************/
int RunSimulation()

{
    /**********************definitions*************************/
    	extern int pdo_run_preset, pgenerate_preset, iteration; 
   extern step *preset_run; 
    extern long pnum_steps, ptotal_steps, poutput_interval, pnum_particles, pnum_types;
    extern double domain_attempts, domain_accepts;
    extern long ptest_interval;
    extern int abort_move;
	extern int nsites, max_mode_radius_squared;
    extern double pvol_moves, pmove_size[NUM_DOF], pfourier_moves;
    extern vector  pbox_length;    
    extern double percent_change, percent_change_samples; 
    double total_energy, volume,random,temp_double, slow_energy;
	movesize *fdelta; 
    int accept, num_presets;
    int temp_accept[NUM_DOF] = {0,0,0,0,0,0,0,0};
    int total_accept[NUM_DOF]={0,0,0,0,0,0,0,0};
    int attempts[NUM_DOF]={0,0,0,0,0,0,0,0};
    int var_index = 0;
    long i,j;
    long next_test = ptest_interval*pnum_particles;
    long *base= &seed;
		int move_samples, total_particle, total_volume; 
    site *the_lattice;
    FILE *fp;
    particle *the_membrane;
    mparticle *the_types; 
    /***********************initialization*********************/
    iteration = 0; 
		total_particle = 0; 
	total_volume = 0; 
	move_samples = 0; 
    if (pdo_run_preset == YES)
    	{printf("going to read in presets.\n"); 
	num_presets =  read_in_preset_steps(&preset_run); 
	fp = fopen("preset_out","w"); 
	fclose(fp); 
	}
	if (pgenerate_preset == YES)
	{
		num_presets = pnum_steps * 10; 
		preset_run = (step *)(malloc(sizeof(step) * num_presets)); 
		for (i = 0; i < num_presets; i++)
			clear_preset_step(&(preset_run[i]));
			fp = fopen("preset_out","w"); 
	fclose(fp); 
	}
	nsites = (int)(pbox_length[X]*1);
	printf("nsites: %d, pbox_length[X]: %lf\n", nsites, pbox_length[X]); 
	
    abort_move = NO_OVERLAP;
    domain_accepts = 0; domain_attempts = 0;
	pfourier_moves = 0; 
    the_membrane = (particle *)malloc(pnum_particles * sizeof(particle));
    if ((the_membrane == NULL))
	  {
        printf("Error, not enough memory: \n");
        fp = fopen("errormsg","w");fclose(fp);return ERROR;
	  }
	
	the_lattice = (site *)malloc((nsites)*(nsites) * sizeof(site));
    for (i = 0; i < (nsites * nsites); i++)
	  {
        the_lattice[i].checklist =  (int *)malloc(nsites*nsites* sizeof(int));
        the_lattice[i].bead_list = NULL;
	  }
	  
    the_types = (mparticle *)malloc(pnum_types * sizeof(mparticle));     
	printf("Memory allocated.\n"); 
    init(the_membrane, the_lattice, the_types, &total_energy, total_accept,pmove_size, ptotal_steps);
	 total_system_energy=total_energy; 
	 slow_total_energy(the_membrane, the_types,  pbox_length);

	//energy_tests(the_membrane, the_lattice, pbox_length);
   // if (fabs(total_energy) >= INFNTY)
	  {
        /*if there's a bug or the initial configuration is bad, abort the simulation*/
     //   stop_run = 1; printf("Stopping due to huge energy of initial conditions!\n");
	  }
    /*store some important values*/
	set_fourier_constants(pbox_length); 
	fdelta = (movesize *)(malloc(sizeof(movesize) * (max_mode_radius_squared + 1))); 
	read_in_fdelta(fdelta); 
	volume = pbox_length[X] *pbox_length[Y]*pbox_length[Z];
	write_version();
	percent_change = 0; 
	percent_change_samples = 0;  
    /********************main loop*****************************/
    for (i = 0;i<pnum_steps;i++)
	  {int move_type; 
		
        if (stop_run == 0){			 
            ptotal_steps++;
            /*determine whether it will be a particle, cluster, or volume move*/
            random = ran3(base)*(pnum_particles-1 + pvol_moves + pfourier_moves);
	    if (random < pnum_particles)
	        { move_type = PARTICLE_MOVE; 
		   total_particle = total_particle + 1; }
	    else if (random < pnum_particles + pvol_moves)
	        {move_type = VOLUME_MOVE;
			total_volume = total_volume + 1; } 
			move_samples = move_samples + 1; 
	    if (pdo_run_preset == YES)
	        move_type = preset_run[iteration].particle_or_volume; 
		
	    if (pgenerate_preset == YES)
			preset_run[iteration].particle_or_volume=move_type; 
			
	    if (iteration > 9307)
			printf("iteration: %d\n", iteration); 
            /*************************particlemove******************************/
            if (move_type == PARTICLE_MOVE)
				{double hold; 
				
				total_energy =  total_energy + particle_move(the_membrane,  the_lattice, the_types, attempts,temp_accept);
				//printf("total_energy:%lf\n", total_energy); 
				//hold = get_total_energy(the_membrane, the_lattice, the_types, pbox_length); 
				//if (hold != total_energy)
					//printf("hold:%lf total_energy:%lf\n", hold, total_energy); 
					}
            /******************volume move******************/
            else if (move_type == VOLUME_MOVE)
			  {
			 // check_bond_energy_difference(the_membrane,the_types,pbox_length);
                abort_move = NO_OVERLAP;
                total_energy = total_energy + volume_step(the_membrane,the_lattice, the_types, &accept,pmove_size[VOLUME], total_energy,pbox_length);
                /*update*/			
		//printf("total_energy:%lf\n", total_energy); 
                temp_accept[VOLUME] = temp_accept[VOLUME] + accept;
                volume =pbox_length[X] * pbox_length[Y]*pbox_length[Z];
                attempts[VOLUME]++;
			  }
			else
			{
				total_energy = total_energy + fourier_step(the_membrane, the_lattice, the_types, &accept, fdelta, total_energy, pbox_length); 
			}
			if (iteration >= num_presets)
				stop_run = YES; 
			 total_system_energy = total_energy; 

            /***********************tests!******************************************/
            if (i == next_test)
			  {
                total_energy = get_total_energy(the_membrane, the_lattice, the_types, pbox_length);
                next_test =  i + ptest_interval*pnum_particles;
                /*adjust move size based on acceptance rates during the last interval*/
               Thermalize(temp_accept, total_accept, pmove_size, attempts); 			
			   adjust_fdelta(fdelta); 	
				write_fdelta(fdelta);
			  }
            /******************************output***************************************/
            if (((i % (poutput_interval * pnum_particles))==0) || (i ==( pnum_steps-1))||(stop_run == 1))
			  {
                write_interval( the_membrane, the_types, total_energy,
                                total_accept, pmove_size, i, pbox_length) ;
 
				slow_energy = energy_tests(the_membrane, the_lattice, the_types,pbox_length);
				printf("slow_energy: %lf\n", slow_energy);
			printf("iterative energy: %lf\n", total_energy);
		//	printf("percent change in bond length: %lf\n", percent_change/percent_change_samples); 
		//	printf("percent change samples: %lf\n", percent_change_samples); 
			printf("percentage of total moves particle moves: %lf\n", total_particle/((double)move_samples)); 
			printf("percentage of total moves volume moves: %lf\n", total_volume/((double)move_samples));
			percent_change = 0; 
			percent_change_samples = 0;               
			  }
        }
	  }/*end main loop*/
	if (pgenerate_preset == YES)
		write_out_preset_steps(preset_run, iteration);
    free_membrane(the_membrane);
    for (i = 0; i < nsites*nsites; i++)
	  {
        free_node_list(the_lattice[i].bead_list);
        free(the_lattice[i].checklist);
	  }
	free(the_lattice); 
	free(fdelta); 
 if ((pdo_run_preset == YES) || (pgenerate_preset == YES))
    	{free(preset_run); 
	}
    printf("Program Complete\n");
}
