#include "structures.h"
#include <stdlib.h>
#include "math.h"
#include "vectorops.h"
#include "accept.h"
#include "preset.h"

/***************************************************************************/
int test_acceptance(double energy_diff)
/*calculates boltzmann weight, compares to random number, and returns 1 (accept) if
the boltzmann factor is greater than the random number. otherwise returns 0 (don't
                                                                             accept)*/
{
	extern int pdo_run_preset, pgenerate_preset, iteration; 
   extern step *preset_run; 
    extern double pkt;
    double boltzmann_factor;
    int accept=0;
    extern long seed;
    long *base = &seed;
    double  random;
	
    random = ran3(base);
    if (pdo_run_preset == YES)
    	random = preset_run[iteration].acceptance_ran; 
	if (pgenerate_preset == YES)
    	preset_run[iteration].acceptance_ran = random; 
    if (energy_diff > 10000)
	  {accept = 0; boltzmann_factor = 10000; }
    else if (energy_diff< 0)
	  {accept = 1;
	  boltzmann_factor = exp(-energy_diff/pkt);}
    else
	  {
        boltzmann_factor = exp(-energy_diff/pkt);
		if (random<boltzmann_factor)
		  {accept = 1;}
        else
            accept = 0;
		if ((accept == 1) && (energy_diff > 10))
			printf("unlikely acceptance\n"); 
	  }
	  if ((pdo_run_preset == YES) || (pgenerate_preset == YES))
	  {
    		preset_run[iteration].energy_diff = energy_diff;
		preset_run[iteration].boltzmann_factor = boltzmann_factor; 
		preset_run[iteration].accept = accept; 
	} 
    return accept;
}


//***********************************************************************/
void CopyParticleContents(  particle *part_to_copy,   particle *part_to_hold)
//copies the contents of a particle into another particle.
{
    int i;
	
    part_to_hold->crossings[X] = part_to_copy->crossings[X];
    part_to_hold->crossings[Y] = part_to_copy->crossings[Y];
    part_to_hold->crossings[Z] = part_to_copy->crossings[Z];
    part_to_hold->energy = part_to_copy->energy;
    part_to_hold->cluster = part_to_copy->cluster;
    part_to_hold->model_index = part_to_copy->model_index; 
   // part_to_hold->radius = part_to_copy->radius;
  //  part_to_hold->length = part_to_copy->length;
    part_to_hold->chain_length = part_to_copy->chain_length;
   // part_to_hold->chain1start = part_to_copy->chain1start; 
    //part_to_hold->chain2start = part_to_copy->chain2start; 
    //part_to_hold->ibindex = part_to_copy->ibindex; 
    //part_to_hold->num_chains = part_to_copy->num_chains; 
    //part_to_hold->cbend = part_to_copy->cbend; 
	//    vcopy(part_to_hold->head, part_to_copy->head);
    for (i = 0; i <part_to_copy->chain_length;i++)
	  {
		copy_bead(&part_to_copy->chain[i], &part_to_hold->chain[i]);
	  }
    //part_to_hold->num_angles = part_to_copy->num_angles; 
    //for (i = 0; i < part_to_copy->num_angles; i++)
    //{
    //	part_to_hold->cosangle[i] = part_to_copy->cosangle[i]; 
    //	part_to_hold->bondangles[i].b1 = part_to_copy->bondangles[i].b1; 
//	part_to_hold->bondangles[i].b2 = part_to_copy->bondangles[i].b2; 
//	part_to_hold->bondangles[i].b3 = part_to_copy->bondangles[i].b3;
//	part_to_hold->bondangles[i].cbend = part_to_copy->bondangles[i].cbend; 
//	part_to_hold->bondangles[i].equil = part_to_copy->bondangles[i].equil;  
  //  }
    //part_to_hold->num_nn = part_to_copy->num_nn;
}

void copy_bead(bead *part_to_copy, bead *part_to_hold)
{
    int i;
    vcopy(part_to_hold->position, part_to_copy->position);
    vcopy(part_to_hold->init, part_to_copy->init); 
    part_to_hold->ci.b = part_to_copy->ci.b;
    part_to_hold->ci.m = part_to_copy->ci.m;
    //part_to_hold->chainpos = part_to_copy->chainpos; 
	//part_to_hold->site_index = part_to_copy->site_index; 
	//part_to_hold->type = part_to_copy->type; 
    //part_to_hold->radius = part_to_copy->radius;
}
/************************************************************/
double CheckAcceptanceRate(int accept, double differential, int interval)
{
    //optimizes acceptance rate, changing move size(passed in as differential) accordingly
    extern   vector pbox_length;
    double acceptance_rate;
    double new_diff=0;
    acceptance_rate = GetAcceptanceRate(accept, interval);
    if ((acceptance_rate)<0.4)
	  {
        new_diff = differential*(0.95);
	  }
    else if ((acceptance_rate>0.6))
	  {
        if ( (4*differential<pbox_length[X])&& (4*differential<pbox_length[Y])&&
             (4*differential<pbox_length[Z]))
            new_diff = differential*(1.05);
        else
            new_diff = differential;
	  }
    else
        new_diff = differential;
    return new_diff;
}



/************************************************************/

double GetAcceptanceRate(int accept, int interval)
{
    double rate;
    if (interval> 0)
        rate = ((double)accept)/((double)interval);
    return rate;
}



/************************************************************/
double GetTotalAcceptanceRate(int accept, int i)
{
    //really no different from GetAcceptanceRate
    double rate;
    double accept2, i2;
    accept2 = (double)(accept);
    i2 = (double)(i+1);
    rate = accept2/i2;
    return rate;
}








