#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "structures.h"
#include "vectorops.h"
#include "preset.h"



int clear_preset_step(step *preset_step)
{
	vclear(preset_step->vector);
	preset_step->scalar = 0; 
	preset_step->scalar2 = 0; 
	preset_step->acceptance_ran = 0; 
	preset_step->energy_diff = 0; 
	preset_step->boltzmann_factor = 0; 
	preset_step->domain_move_type = 0; 
	preset_step->regular_move_type = 0; 
	preset_step->bond_index = 0; 
	preset_step->coord = 0; 
	preset_step->m = 0; 
	preset_step->b = 0; 
	preset_step->accept = 0; 
	preset_step->particle_or_volume = 0; 
	preset_step->acceptance_ran = 0; 
}

int write_out_preset_steps(step *preset_step, int num_presets)
{
FILE *fp; 
int  i, p_or_v, rmove, dmove, m, b, bond_index, coord; 
double scalar, scalar2, x, y, z; 
fp = fopen("preset","w"); 
if (fp == NULL)
{
	printf("can't open preset file\n"); 
	return 0; 
}
fprintf(fp, "num_presets %d\n", num_presets); 
fprintf(fp, "step p_or_v rmove dmove m b bond coord scalar scalar2 x y z ran\n"); 
for (i = 0; i < num_presets; i++)
{
	fprintf(fp,"%d ",i); 
	fprintf(fp,"%d ",preset_step[i].particle_or_volume); 
	fprintf(fp,"%d ",preset_step[i].regular_move_type); 
	fprintf(fp,"%d ",preset_step[i].domain_move_type); 
	fprintf(fp,"%d ",preset_step[i].m); 
	fprintf(fp,"%d ",preset_step[i].b); 
	fprintf(fp,"%d ",preset_step[i].bond_index); 
	fprintf(fp,"%d ",preset_step[i].coord); 
	fprintf(fp,"%22.18lf ",preset_step[i].scalar); 
	fprintf(fp,"%22.18lf ",preset_step[i].scalar2); 
	fprintf(fp,"%22.18lf ",preset_step[i].vector[X]); 
	fprintf(fp,"%22.18lf ",preset_step[i].vector[Y]); 
	fprintf(fp,"%22.18lf ",preset_step[i].vector[Z]); 
	fprintf(fp,"%22.18lf\n",preset_step[i].acceptance_ran); 
	}
fclose(fp); 
}


int read_in_preset_steps(step **preset_step)
{
FILE *fp; 
char temp[201]; 
int num_presets, i, p_or_v, rmove, dmove, m, b, bond_index, coord; 
double scalar, scalar2, x, y, z, aran; 
printf("got to read presets\n"); 
fp = fopen("preset","r"); 
if (fp == NULL)
{
	printf("can't open preset file\n"); 
	return 0; 
}
fscanf(fp, "%*s %d\n", &num_presets); 
 printf("going to read %d presets.\n", num_presets); 
*preset_step = (step *)(malloc(sizeof(step) * num_presets)); 
printf("allocated room for preset.\n"); 
fgets(temp, sizeof(temp),fp);
puts(temp); 
for (i = 0; i < num_presets; i++)
{
	clear_preset_step(&((*preset_step)[i])); 
	fscanf(fp,"%*d %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf\n", &p_or_v, &rmove, &dmove, &m, &b, &bond_index,  &coord, &scalar, &scalar2, &x, &y, &z, &aran); 
	(*preset_step)[i].particle_or_volume = p_or_v; 	
//	printf("p_or_v: %d\n", p_or_v); 
	(*preset_step)[i].regular_move_type = rmove; 
	(*preset_step)[i].domain_move_type = dmove; 
	(*preset_step)[i].m = m; 
	(*preset_step)[i].b = b; 
	(*preset_step)[i].bond_index = bond_index; 
	(*preset_step)[i].coord = coord; 
	(*preset_step)[i].scalar = scalar; 
	(*preset_step)[i].scalar2 = scalar2; 
	(*preset_step)[i].vector[X] = x; 
	(*preset_step)[i].vector[Y] = y; 
	(*preset_step)[i].vector[Z] = z; 
	(*preset_step)[i].acceptance_ran = aran; 
}
fclose(fp); 
printf("finished reading presets\n"); 
return num_presets; 
}

int output_preset_results(step preset_step)
{
	extern int iteration; 
	FILE *fp; 
	fp = fopen("preset_out","a"); 
	fprintf(fp, "%d %d %10.6lf %10.6lf\n", iteration, preset_step.accept, preset_step.energy_diff, preset_step.boltzmann_factor); 
	fclose(fp); 
}