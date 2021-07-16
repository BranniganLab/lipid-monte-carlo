/***************************************************************************
bounds.c  -  description
-------------------
    begin                : Mon Jan 20 2003
    copyright            : (C) 2003 by Grace Brannigan
    email                : gbrannig@physics.ucsb.edu
    ***************************************************************************/

/***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "structures.h"
#include "vectorops.h"
#include "bounds.h"

/******************************************************************/
int translate_periodic(particle *the_particle, mparticle  *the_type, vector translate, vector box_length)
{
    vector temp, wrap;
    int i;
    int ibindex; 
    ibindex = the_type->ibindex; 
    vadd(temp, the_particle->chain[ibindex].position, translate);
    wrap_around_vector(wrap, temp, box_length);
    vcopy(the_particle->chain[ibindex].position, temp);

    vsub(translate, translate, wrap);
   for (i = 0; i < ibindex; i++)
   {
   vadd(temp, the_particle->chain[i].position, translate);
   vcopy(the_particle->chain[i].position, temp);
   }
    for (i =ibindex+1; i < the_particle->chain_length; i++)
    {
        vadd(temp, the_particle->chain[i].position, translate);
        vcopy(the_particle->chain[i].position, temp);
    }

    //        get_tail(the_particle->tail, *the_particle, box_length);

}

/**************************************************************/
int wrap_around_vector(vector wrap, vector position, vector box_length)
{
    int shift;
    int i;
    vector temp_position;

    vcopy(temp_position, position);
    vclear(wrap);
    for (i = 0; i< DIMENSION; i++)
    {
        shift = 0;
        while (position[i] < 0)
      		{position[i] = position[i] + box_length[i];
            shift = shift - 1;}
        while ((position[i]) > box_length[i])
        {position[i] = position[i] - box_length[i];
            shift = shift + 1;}
        wrap[i] = shift * box_length[i];               }
    return IN_BOUNDS;

}

/**************************************************************/
int wrap_around(vector wrap, particle *the_particle, mparticle *the_type, vector box_length)
{
    int shift;
    int i, ibindex;
    vector temp_position;
    ibindex = the_type->ibindex; 
    vcopy(temp_position, the_particle->chain[1].position);
    vclear(wrap);
    for (i = 0; i< DIMENSION; i++)
    {
        shift = 0;
        while (the_particle->chain[ibindex].position[i] < 0)
      		{the_particle->chain[ibindex].position[i] = the_particle->chain[ibindex].position[i] + box_length[i];
            shift = shift - 1;}
        while ((the_particle->chain[ibindex].position[i]) > box_length[i])
        {the_particle->chain[ibindex].position[i] = the_particle->chain[ibindex].position[i] - box_length[i];
            shift = shift + 1;}
        wrap[i] = shift * box_length[i];               }
    return IN_BOUNDS;

}

/************************************************************/
int check_bounds(particle part, mparticle type, vector box_length)
{
	int ibindex; 
	ibindex = type.ibindex; 
	
        if (part.chain[ibindex].position[X] < 0)
        return OVERLAP;
        if (part.chain[ibindex].position[X] >box_length[X])
        return OVERLAP;
       if (part.chain[ibindex].position[Y] < 0)
        return OVERLAP;
        if (part.chain[ibindex].position[Y] >box_length[Y])
        return OVERLAP;
        if (part.chain[ibindex].position[Z] < 0)
        return OVERLAP;
        if (part.chain[ibindex].position[Z] >box_length[Z])
        return OVERLAP;
        return NO_OVERLAP;
}
/************************************************************/
void print_out_of_bounds(particle *the_membrane, mparticle *the_types, vector box_length)
{
    extern long pnum_particles;
    int i, ibindex;
    printf("*********OUT OF BOUNDS PARTICLES: ***********\n");
    for (i = 0; i< pnum_particles; i++)
    {
    	ibindex = the_types[the_membrane[i].model_index].ibindex; 
        if (the_membrane[i].chain[ibindex].position[X] < 0)
        {printf("particle %d: ", i);
            vprint(the_membrane[i].chain[ibindex].position);}
        else if (the_membrane[i].chain[ibindex].position[Y] < 0)
        {printf("particle %d: ", i);
            vprint(the_membrane[i].chain[ibindex].position);}
        else if (the_membrane[i].chain[ibindex].position[Z] < 0)
        {printf("particle %d: ", i);
            vprint(the_membrane[i].chain[ibindex].position);  }

    }
    printf("********************\n");
}