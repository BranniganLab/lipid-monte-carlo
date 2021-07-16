#include <stdlib.h>
#include <math.h>
#include "structures.h"
#include "vectorops.h"
#include "perturb.h"
#include "model.h"
#include "utils.h"
#include "bounds.h"
#include "domain.h"
#include "preset.h"
/******************************************************************************/
int GetPertIndex()
{
    //chooses a random particle index.
    int  test, pert_index; /*holds index of changed particle*/
    extern long seed;
    extern long pnum_particles;
    extern int pdo_run_preset, pgenerate_preset, iteration; 
    extern step *preset_run; 
    long *base= &seed;
	
    test = (int)((ran3(base)*(pnum_particles)));
    if (test<pnum_particles)
        pert_index = test;
    else
        pert_index = pnum_particles-1;
    if (pdo_run_preset == YES)
       pert_index = preset_run[iteration].m; 
	if (pgenerate_preset == YES) 
		preset_run[iteration].m = pert_index; 
    return pert_index;
}

int perturb_domain(domain *the_domain, particle *the_membrane, mparticle *the_types, vector translation, vector rotation, double delta, vector box_length)
{
	extern long seed;
	extern double pdomaindx, pdomaindr; 
	extern vector pdfixed; 
	extern int pdomain_fixed; 
	extern int pdomain_constrain; 
	extern int pdo_run_preset, pgenerate_preset, iteration; 
        extern step *preset_run; 
    long *base= &seed;
    double size;
	vector wrap; 
	int coord; 
	
	if (((pdomain_constrain == NO) || (pdomain_fixed == 0)) || (pdomain_fixed >MAX_DOMAINS))
	{
		size= ran3(base)-0.5;
		translation[X] = size*pdomaindx;
		size= ran3(base)-0.5;
		translation[Y] = size*pdomaindx;
		size= ran3(base)-0.5;
		translation[Z] = size*pdomaindx;
		if (pdo_run_preset == YES)
			vcopy(translation, preset_run[iteration].vector); 
		if (pgenerate_preset == YES)
			vcopy(preset_run[iteration].vector, translation);
	}
	else
	{
		size = ran3(base) - 0.5; 
		if (pdo_run_preset == YES)
			size = preset_run[iteration].scalar; 
		if (pgenerate_preset == YES)
			preset_run[iteration].scalar = size; 
		vsub_periodic(translation, pdfixed, the_domain->center, box_length); 
		vmult_scal(translation, translation, size); 
	}
	vadd(the_domain->center, the_domain->center, translation); 
	wrap_around_vector(wrap, the_domain->center, box_length); 
	translate_domain(the_domain, the_membrane, the_types, translation, box_length); 
	
	coord = (int)(3*ran3(base)); 
	vclear(rotation); 
	rotation[coord] = pdomaindr * ran_theta(base); 
	if (pdo_run_preset == YES)
		rotation[preset_run[iteration].coord] = preset_run[iteration].scalar2; 
	if (pgenerate_preset == YES)
		{preset_run[iteration].coord = coord;
		 preset_run[iteration].scalar2 = rotation[coord];
		 }
	rotate_domain(the_domain, the_membrane, rotation, box_length); 

}

int rotate_interface(particle *pert_part, mparticle *pert_part_type, double delta)
{extern int pdo_run_preset, pgenerate_preset,  iteration; 
        extern step *preset_run; 
    
    vector axis, p, r, n, newr, linear;
    double dphi;
    extern long seed;
    long *base= &seed;
    double random; 
    int fixedstart, rotatestart, rotatestop, ibindex, i;  
    ibindex = pert_part_type->ibindex; 
    random= ran3(base);
    if (pdo_run_preset == YES)
	random = preset_run[iteration].scalar; 
	if (pgenerate_preset == YES)
	preset_run[iteration].scalar = random; 
    if (random <0.5)
    	{fixedstart = pert_part_type->chain1start;
	rotatestart = pert_part_type->chain2start;
	rotatestop = pert_part_type->chain_length - 1; }  
    else
	{fixedstart = pert_part_type->chain2start; 
	rotatestart = pert_part_type->chain1start; 
	rotatestop = pert_part_type->chain2start - 1; }
	
    vsub(axis, pert_part->chain[ibindex - 1].position, pert_part->chain[fixedstart].position);
    vcopy(p, pert_part->chain[ibindex].position);
    vsub(r, p, pert_part->chain[ibindex-1].position);
    vunit(n, axis);
    dphi = delta * ran_theta(base); 
    if (pdo_run_preset == YES)
	dphi = preset_run[iteration].scalar2;
	if (pgenerate_preset == YES)
	preset_run[iteration].scalar2 = dphi;
    vrotate(newr, r, n, dphi);
    vadd(pert_part->chain[ibindex].position, pert_part->chain[ibindex-1].position, newr);
    vsub(linear, pert_part->chain[ibindex].position, p); 
    for (i = rotatestart; i <=rotatestop; i++)
    {
    	vadd(pert_part->chain[i].position, pert_part->chain[i].position, linear); 
    }
    
    /*rotates chain with interface bead. probably not going to be accepted very often because it can \
    result in huge linear displacements, depending on the orientation of the chain.*/
 /* 
    for (i = rotatestart; i <= rotatestop; i++)
    	{
	    vcopy(p, pert_part->chain[i].position);
    	    vsub(r, p, pert_part->chain[ibindex-1].position);
	    vrotate(newr, r, n, dphi);
    	    vadd(pert_part->chain[i].position, pert_part->chain[ibindex-1].position, newr);
	}
*/
}

/******************************************************************************************/

int pivot_interface(particle *pert_part, mparticle *pert_part_type, double delta)
{
    extern long seed;
    int i;
    long *base= &seed;
    vector v, size; 
    double cosangle; 
    vector axis, p, r, n, newr; 
    double dphi; 
    int count;
    extern int pdo_run_preset, pgenerate_preset,  iteration; 
        extern step *preset_run; 
		
    vsub(axis, pert_part->chain[pert_part_type->chain1start].position, pert_part->chain[pert_part_type->chain2start].position);
    vcopy(p, pert_part->chain[pert_part_type->ibindex].position);
    vsub(r, p, pert_part->chain[pert_part_type->chain1start].position);
    vunit(n, axis);
    dphi = delta * ran_theta(base); 
	if (pdo_run_preset == YES)
	dphi = preset_run[iteration].scalar;
	if (pgenerate_preset == YES)
	preset_run[iteration].scalar = dphi;
    vrotate(newr, r, n, dphi);
    vadd(pert_part->chain[pert_part_type->ibindex].position, pert_part->chain[pert_part_type->chain1start].position, newr);
}




int rotate_bead(particle *pert_part, mparticle *pert_part_type, double delta, int bead)
{
    extern int pdo_run_preset, pgenerate_preset, iteration; 
    extern step *preset_run; 
    vector axis, p, r, n, newr;
    double dphi;
    int refbead; 
    extern long seed;
    long *base= &seed;
    refbead = bead - 1; 
    if (bead == pert_part_type->chain2start)
    	refbead = pert_part_type->ibindex; 
	
	vsub(axis, pert_part->chain[refbead].position, pert_part->chain[bead + 1].position);
	vcopy(p, pert_part->chain[bead].position);
    vsub(r, p, pert_part->chain[refbead].position);
    vunit(n, axis);
    dphi = delta * ran_theta(base); 
    if (pdo_run_preset == YES)
    	dphi = preset_run[iteration].scalar; 
	if (pgenerate_preset == YES)
    	preset_run[iteration].scalar = dphi; 
    vrotate(newr, r, n, dphi);
	vadd(pert_part->chain[bead].position, pert_part->chain[refbead].position, newr);
		
}

/******************************************************************************************/

int pivot_end(particle *pert_part, mparticle *pert_part_type, double delta, int end)
{
    extern long seed;
    extern int pdo_run_preset, pgenerate_preset, iteration; 
    extern step *preset_run; 
    int i;
    long *base= &seed;
    vector v, size; 
    double oldlength; 
    for (i = 0; i < DIMENSION; i++)
        size[i]= delta*(ran3(base)-0.5);
    
    if (pdo_run_preset == YES)
        vcopy(size, preset_run[iteration].vector); 
	if (pgenerate_preset == YES)
        vcopy(preset_run[iteration].vector, size); 

    	
    if (end == HEAD)
	  {
        vsub(v, pert_part->chain[0].position, pert_part->chain[1].position);
	oldlength = sqrt(vsquare(v)); 
        vadd(v, v, size);
        vresize(v,v, oldlength );
        vadd(pert_part->chain[0].position, pert_part->chain[1].position, v);
	  }
    else if ((end == TAIL2)|| (pert_part_type->num_chains ==1))
	  {
		vsub(v, pert_part->chain[pert_part->chain_length-1].position, pert_part->chain[pert_part->chain_length-2].position);
		oldlength = sqrt(vsquare(v)); 
        vadd(v, v, size);
        vresize(v, v,oldlength);
        vadd(pert_part->chain[pert_part->chain_length-1].position, pert_part->chain[pert_part->chain_length-2].position, v);
	  }
   else if (end == TAIL)
	  {
		vsub(v, pert_part->chain[pert_part_type->chain2start-1].position, pert_part->chain[pert_part_type->chain2start-2].position);
				oldlength = sqrt(vsquare(v)); 
        vadd(v, v, size);
        vresize(v, v,oldlength);
        vadd(pert_part->chain[pert_part_type->chain2start-1].position, pert_part->chain[pert_part_type->chain2start-2].position, v);
	  }
}



/******************************************************************************/
int perturb_chain(particle *pert_part,mparticle *pert_part_type, int bead_index, double delta, vector box_length)
{
    /*change one bond angle in a chain*/
    int i;
    int nbeads;
    int count, chain1stop; 
    double cosangle; 
    vector empty_vect; 
	nbeads = pert_part->chain_length; 
    if ((bead_index == pert_part_type->ibindex) && (pert_part_type->num_chains == 2))
    {
    	if (bead_index == 0)
		pivot_interface(pert_part, pert_part_type, delta); 
	else
		return; 
		//rotate_interface(pert_part, delta); 
    }
    else if (bead_index == 0)
        pivot_end(pert_part, pert_part_type, delta, HEAD);
    else if (pert_part_type->num_chains == 1)
    	{if (bead_index == (nbeads - 1))
	 	pivot_end(pert_part, pert_part_type, delta,TAIL);
	else
		rotate_bead(pert_part, pert_part_type, delta, bead_index); 
	}
    else if (pert_part_type->num_chains == 2)
    {
    	if (bead_index == (nbeads - 1))
		pivot_end(pert_part, pert_part_type,  delta, TAIL2); 
	else if (bead_index == (pert_part_type->chain2start - 1))
		pivot_end(pert_part, pert_part_type, delta, TAIL); 
	else
		rotate_bead(pert_part, pert_part_type, delta, bead_index); 
	
    }
    
    vclear(empty_vect); 
	if (bead_index == 1)
		translate_periodic(pert_part, pert_part_type,empty_vect, box_length); 
    count = 0;
    for (i = 0; i < pert_part_type->num_angles; i++)
    	{cosangle = get_angle(*pert_part, pert_part_type->bondangles[i].b1,  pert_part_type->bondangles[i].b2,pert_part_type->bondangles[i].b3);
	pert_part->cosangle[i] = cosangle; 
	}
/*chain1stop = pert_part->chain2start; 
if (chain1stop > pert_part->chain_length)
	chain1stop = pert_part->chain_length; 
    for (i = 0; i < chain1stop - 2;i++)
	  {
		cosangle = get_angle(*pert_part, i, i+ 1, i+ 2);
		pert_part->cosangle[count] = cosangle;
		count++; 
	  }
    if (pert_part->chain2start < nbeads)
    {
    cosangle = get_angle(*pert_part, 0, pert_part->chain2start, pert_part->chain2start + 1);
    pert_part->cosangle[count] = cosangle; 
    count++;
    for (i = pert_part->chain2start; i < pert_part->chain_length - 2; i++)
    {
	cosangle = get_angle(*pert_part, i, i+ 1, i+ 2);
	pert_part->cosangle[count] = cosangle;
	count++; 
	}
    }	*/		
}
/*******************************************************************************/
int perturb_bond(particle *pert_part, mparticle *pert_type, int bond_to_perturb, double delta, vector box_length)
{       
	extern int pdo_run_preset, pgenerate_preset, iteration; 
	extern step *preset_run; 
	extern long seed; 
	long *base= &seed;
	int b1, b2; 
	double size; 
	vector bond, bondperturb; 
	if (pert_type->bonds[bond_to_perturb].constrained == YES)
		return NO; 
	size = delta*(ran3(base)-0.5);
	if (pdo_run_preset == YES)
		size = preset_run[iteration].scalar;
	if (pgenerate_preset == YES)
		preset_run[iteration].scalar = size;

	b1 = pert_type->bonds[bond_to_perturb].b1;
	b2 = pert_type->bonds[bond_to_perturb].b2;
	vsub_periodic(bond, pert_part->chain[b1].position, pert_part->chain[b2].position, box_length); 
	vmult_scal(bondperturb, bond, size);
	vadd(pert_part->chain[b1].position, pert_part->chain[b1].position, bondperturb); 
	vmult_scal(bondperturb, bondperturb, -1); 
	vadd(pert_part->chain[b2].position, pert_part->chain[b2].position, bondperturb); 
	return YES; 
}
/*******************************************************************************/
void uniform_perturb_position(particle *pert_part, mparticle *pert_type, vector translation, double delta,
							  vector box_length )
{
    //generates a particle translation
    extern int pdo_run_preset, pgenerate_preset,  iteration; 
    extern step *preset_run;
    extern long seed;
    long *base= &seed;
    double size;
	

    size= ran3(base)-0.5;
    translation[X] = size*delta;
    size= ran3(base)-0.5;
    translation[Y] = size*delta;
    size= ran3(base)-0.5;
    translation[Z] = size*delta;
    if (pdo_run_preset == YES)
	vcopy(translation, preset_run[iteration].vector);
	if (pgenerate_preset == YES)
	vcopy(preset_run[iteration].vector, translation);
    translate_periodic(pert_part, pert_type, translation, box_length);
}
/*******************************************************************************/
void uniform_perturb_z(particle *pert_part,mparticle *pert_part_type, vector translation, double delta,
							  vector box_length )
{
    //generates a particle translation
    extern int pdo_run_preset, pgenerate_preset, iteration; 
    extern step *preset_run;
    extern long seed;
    long *base= &seed;
    double size;
    size= ran3(base)-0.5;
    translation[Z] = size*delta;
    if (pdo_run_preset == YES)
    	translation[Z] = preset_run[iteration].scalar; 
	if (pgenerate_preset == YES)
    	preset_run[iteration].scalar = translation[Z]; 
    translate_periodic(pert_part,pert_part_type, translation, box_length);
	
}

/*******************************************************************************/
void uniform_perturb_radius(particle *pert_part, mparticle *pert_type, vector translation, double delta,
							  vector box_length, uconstraint the_constraint)
{
    //generates a particle translation
	extern int pdo_run_preset, pgenerate_preset, iteration; 
    extern step *preset_run;
    extern long seed;
    long *base= &seed;
    double size;
    size= delta*(ran3(base)-0.5);
	the_constraint.arrow[Z] = 0;
    if (pdo_run_preset == YES)
    	size = preset_run[iteration].scalar;  
	if (pgenerate_preset == YES)
    	preset_run[iteration].scalar = size;  
    vmult_scal(translation, the_constraint.arrow, -size); 
    translate_periodic(pert_part, pert_type, translation, box_length);
	size = delta*(ran3(base) - 0.5); 
if (pdo_run_preset == YES)
    	size = preset_run[iteration].scalar2; 
	if (pgenerate_preset == YES)
    	preset_run[iteration].scalar2 = size;  
	vclear(translation); translation[Z] = size; 
	translate_periodic(pert_part, pert_type, translation, box_length); 
}
/*******************************************************************************/
int perturb_particle(particle *pert_part, mparticle *pert_part_type, chainindex ci, vector translation, double delta, vector box_length, int code)
{
	int i, pert_index, bead_index, constraint_index; 
	extern int pnumber_user_constraints; 
	extern uconstraint *puser_constraints; 
	
	pert_index = ci.m; 
	bead_index = ci.b; 
	
	constraint_index = -1; 
	for (i = 0; i < pnumber_user_constraints; i++)
	{
		if (puser_constraints[i].ci.m==pert_index)
			if ((puser_constraints[i].ci.b==bead_index) || ((code == PERTURB_WHOLE_CHAIN)||(puser_constraints[i].type == FIXED_POINT)))
				constraint_index = i; 
		}
	if (constraint_index == -1)
	{ 
		if (code == PERTURB_WHOLE_CHAIN)
			uniform_perturb_position(pert_part, pert_part_type, translation, delta, box_length); 
		else if (code == PERTURB_BEAD)
			perturb_chain(pert_part, pert_part_type, bead_index, delta, box_length);
		else if (code == PERTURB_BOND)
			return perturb_bond(pert_part, pert_part_type, bead_index, delta,  box_length);  
		return YES; 
	}
	
		if (code == PERTURB_BEAD)
		return NO; 

	if (puser_constraints[constraint_index].type == FIXED_ANGLE)
	{
		uniform_perturb_radius(pert_part, pert_part_type, translation, delta, box_length, puser_constraints[constraint_index]); 
		}
	else
	{
		uniform_perturb_z(pert_part, pert_part_type, translation,delta,box_length);  
	}
	return YES; 	
}
