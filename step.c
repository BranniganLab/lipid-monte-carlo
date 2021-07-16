#include <stdlib.h>
#include "structures.h"
#include "accept.h"
/*****************************************************************/
void Thermalize(int temp_accept[], int total_accept[], double move_size[], int interval[])
/*for each degree of freedom, adjusts the move size according to the acceptance
rate*/

{
    int i;
    extern int pequilibrate;
	extern double pdomaindx, pdomaindr, domain_accepts, domain_attempts; 

    for (i = 0; i<(NUM_DOF); i++)
	  { if (interval[i] > 0)	/*as long as at least one move has been made*/
		{
		  if ((pequilibrate == YES))
			  move_size[i] = CheckAcceptanceRate(temp_accept[i],move_size[i], interval[i]);
		  if (move_size[i] < 0.05)
			  move_size[i] = 0.05;
		  total_accept[i] = total_accept[i] + temp_accept[i];
		  temp_accept[i] = 0;
		  interval[i] = 0;
        }
	  }
	pdomaindx = CheckAcceptanceRate(domain_accepts, pdomaindx, domain_attempts); 
	pdomaindr = CheckAcceptanceRate(domain_accepts, pdomaindr, domain_attempts); 
//printf("\n");

}

/********************************************************************************/
void ShuffleMoveClass(int move_class[])
/*this function just mixes up all the individual particle perturbations so they
don't happen each time. Every iteration, the move class (X,Y,Z, ROTATION or FLIP)
is shuffled, and the program perturbs each degree of freedom in the order they
appear in move_class[]. This function also determines whether a flip or rotation
will occur*/
{
    extern long seed;
    long *base= &seed;
    extern long pnum_particles;
    extern double pflip_moves;
    double random;
    int i, index, temp;
	
    for (i = 0; i<= ROTATION; i++)
	  {
        temp = move_class[i];
		
        /*pick an index at random, where temp will go*/
        index= (int)(ROTATION*ran3(base));
		
        /*if it's a flip, decide whether to switch to a rotation and vice versa*/
        if ((temp == FLIP )|| (temp == ROTATION) )
		  {random = (ran3(base)) * (pnum_particles-1) + (pflip_moves);
            if (random <pnum_particles)
                temp = ROTATION;
            else
                temp = FLIP; }
		
        /*swap the moves held in i and index*/
        move_class[i] = move_class[index];
        move_class[index] = temp;
	  }
	
}

double anneal(double current_temp, double goal_temp, double fluct, double interval)
{
    extern long seed;
    long *base= &seed;
    double temp_diff, dT, random;
    temp_diff = current_temp - goal_temp;
    dT = -temp_diff*interval/1000;
    random = (ran3(base) - 0.5) * fluct;
    dT = dT + random;
    return (current_temp + dT);
}

