/***************************************************************************
mpref.c  -  description
-------------------
    begin                : Mon Nov 25 2002
    copyright            : (C) 2002 by Grace Brannigan
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
/*pref.c*/
/*commented, tentatively released 9/1/02, some obsolete vars can be removed*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "structures.h"
#include "vectorops.h"
#include "init.h"
#include "model.h"

/*************************preference variables********************/

long pnum_steps, ptotal_steps, poutput_interval, pnum_particles, pmax_steps, ptest_interval, pnum_types;
vector  pbox_length;
vector pdfixed; 
double paspect, ptension,pcluster_moves, ppatch_moves, pvol_moves, pflip_moves,pswap_moves, pcluster_length, ppatch_length,  pkt, pgoal_kt, pmove_size[NUM_DOF], ppressure;
int pequilibrate, pgroup_move, pboundary;
char pinfile[30],poutfile[30],plogfile[30],pcrdfile[30],psffile[30],pdbfile[30],ptrjfile[30], prefix[30] ;
int pcore;
double pdiameter, plength, ptorsion;
double pdomaindx, pdomaindr, pumbrella_coeff,pumbrella_center; 
potential pmodel[10];
double domain_accepts, domain_attempts;
int pdomain_fixed, pdomain_constrain, pnumber_user_interactions,pnumber_user_constraints;
upotential *puser_potentials;  
uconstraint *puser_constraints; 
mtablesite pmodeltable[3][3]; 
int pdo_run_preset, pgenerate_preset; 

char *prefile = "mc.pref";

void write_sim_prefs()
{

    extern long pnum_steps, ptotal_steps, poutput_interval, pnum_particles, ptest_interval;
    extern  vector  pbox_length;
    extern double pcluster_moves, pvol_moves, pflip_moves, pcluster_length, pkt, pgoal_kt, pmove_size[], ptension;
    extern int pequilibrate, pboundary, preshaping_method, pdo_run_preset, pgenerate_preset;
    extern char prefix[30];
    FILE *fp;
    if ((fp = fopen("sim.pref","w"))==NULL)
    {puts("cannot open sim.pref for writing\n");
    }
    fprintf(fp, "number of steps %ld\n", pnum_steps);
    fprintf(fp, "total steps so far %ld\n", ptotal_steps);
    fprintf(fp, "current kT %lf\n", pkt);
    fprintf(fp, "goal kT %lf\n", pgoal_kt);
    fprintf(fp, "output interval %ld\n", poutput_interval);
    fprintf(fp, "acceptance test inverval %ld\n", ptest_interval);
    fprintf(fp, "thermalize?  (0=no,1=yes) %d\n", pequilibrate);
    fprintf(fp, "number of particles %ld\n", pnum_particles);
    fprintf(fp, "number of types %ld\n", pnum_types); 
    fprintf(fp, "boxlength.x %lf\n", pbox_length[X]);
    fprintf(fp, "boxlength.y %lf\n", pbox_length[Y]);
    fprintf(fp, "boxlength.z %lf\n", pbox_length[Z]);
    fprintf(fp, "tangential pressure %lf\n", ptension);
    fprintf(fp,"volume changes/MC step %lf\n", pvol_moves);
    fprintf(fp,"cluster attempts/MC step %lf\n", pcluster_moves);
    fprintf(fp,"patch attempts/MC step %lf\n", ppatch_moves);
    fprintf(fp,"flip tries/MC step %lf\n", pflip_moves);
    fprintf(fp, "cluster length scale %lf\n", pcluster_length);
    fprintf(fp, "patch length scale %lf\n", ppatch_length);
    fprintf(fp, "boundary (0=none,1=meanfield,2=periodic,3=plane) %d\n", pboundary);
    fprintf(fp, "dX %lf\n", pmove_size[X]);
    fprintf(fp, "dY %lf\n", pmove_size[Y]);
    fprintf(fp, "dZ %lf\n", pmove_size[Z]);
    fprintf(fp, "dR %lf\n", pmove_size[ROTATION]);
    fprintf(fp, "dA[X] %lf\n", pmove_size[VOLUME+X]);
    fprintf(fp, "dA[Y] %lf\n", pmove_size[VOLUME+Y]);
    fprintf(fp, "dC %lf\n", pmove_size[CLUSTER]);
    fprintf(fp, "dP %lf\n", pmove_size[PATCH]);
    fprintf(fp, "%s %20s\n", "fileprefix ",prefix);
    fprintf(fp, "%s %d\n", "preset_run?", pdo_run_preset); 
	fprintf(fp, "%s %d\n", "preset_generate?", pgenerate_preset); 
    fprintf(fp, "reshaping method: (0=possibly_constrained,1=definitely_constrained) %d\n", preshaping_method); 
    fclose(fp);
}

void write_model_prefs()
/*writes the current "preferences" out- some of these change during the course
of the simulation(i.e., move size, num_steps, box_size, etc...)- overwrites the
original pref file*/
{

    FILE *fp;

    extern int pcore;
    extern double pdiameter, plength;
    extern potential pmodel[];
    extern double ptorsion;
    int i;

    if ((fp = fopen("model.pref","w"))==NULL)
    {puts("cannot open model.pref for writing\n");
    }
    else{
        fprintf(fp, "sphere_diameter %lf\n", pdiameter);
        fprintf(fp, "length %lf\n", plength);
        fprintf(fp, "core(0=hard,1=soft) %d\n", pcore);
        fprintf(fp, "number of alignment terms %d\n", pmodel[HINTERFACE].num_terms);
        fprintf(fp, "range of alignment terms %lf\n", pmodel[HINTERFACE].range);
        for (i = 0; i<pmodel[HINTERFACE].num_terms; i++)
            {fprintf(fp,"%lf %lf\n", (pmodel[HINTERFACE].coefficients[i][0]), (pmodel[HINTERFACE].coefficients[i][1]));
               }
        fprintf(fp, "number of attraction terms %d\n", pmodel[HTAIL].num_terms);
        fprintf(fp, "range of attraction terms %lf\n", pmodel[HTAIL].range);

        for (i = 0; i<pmodel[HTAIL].num_terms; i++)
            {fprintf(fp,"%lf %lf\n", (pmodel[HTAIL].coefficients[i][0]), (pmodel[HTAIL].coefficients[i][1]));    }
        fprintf(fp, "number of softcore terms %d\n", pmodel[HSOFTCORE].num_terms);
        fprintf(fp, "range of softcore terms %lf\n", pmodel[HSOFTCORE].range);
        for (i = 0; i<pmodel[HSOFTCORE].num_terms; i++)
            fprintf(fp,"%lf %lf\n", (pmodel[HSOFTCORE].coefficients[i][0]), (pmodel[HSOFTCORE].coefficients[i][1]));
        fprintf(fp,"torsion %lf\n", ptorsion);
        fclose(fp);
    }
}

void get_sim_prefs()
/*reads in the preferences from "prefile"*/
{
    extern long pnum_steps, ptotal_steps, poutput_interval, pnum_particles, ptest_interval;
    extern  vector  pbox_length;
    extern double pcluster_moves, pvol_moves, pflip_moves, pcluster_length, ppatch_length,pkt, pgoal_kt, pmove_size[], ptension, paspect;
    extern char pinfile[30],poutfile[30],plogfile[30],pcrdfile[30],psffile[30],pdbfile[30],ptrjfile[30], prefix[30] ;
    extern int pequilibrate, pboundary, preshaping_method, pdo_run_preset, pgenerate_preset;
    FILE *fp;

    if ((fp = fopen("sim.pref","r"))==NULL)
    {puts("cannot open preferences\n");
    }
    else
    {
        init_potential(pmodel[HSOFTCORE]);
        fscanf(fp,"%*s %*s %*s %ld\n",&pnum_steps);
        fscanf(fp,"%*s %*s %*s %*s %ld\n", &ptotal_steps);
        fscanf(fp,"%*s %*s %lf\n", &pkt);
        fscanf(fp,"%*s %*s %lf\n", &pgoal_kt);
        fscanf(fp,"%*s %*s %ld\n", &poutput_interval);
        fscanf(fp,"%*s %*s %*s %ld\n", &ptest_interval);
        fscanf(fp,"%*s %*s %d\n", &pequilibrate);
        fscanf(fp,"%*s %*s%*s %ld\n", &pnum_particles);
	fscanf(fp,"%*s %*s%*s %ld\n", &pnum_types); 
        fscanf(fp,"%*s %lf\n", &(pbox_length[X]));
        fscanf(fp,"%*s %lf\n",&(pbox_length[Y]));
        paspect = pbox_length[Y]/pbox_length[X];
        fscanf(fp,"%*s %lf\n", &(pbox_length[Z]));
        fscanf(fp,"%*s %*s%lf\n", &ptension);
        fscanf(fp,"%*s %*s %*s%lf\n", &pvol_moves);
        fscanf(fp,"%*s %*s %*s%lf\n", &pcluster_moves);
        fscanf(fp,"%*s %*s %*s%lf\n", &ppatch_moves);
        fscanf(fp,"%*s %*s %*s%lf\n", &pflip_moves);
        fscanf(fp,"%*s %*s %*s%lf\n", &pcluster_length);

        fscanf(fp,"%*s %*s %*s%lf\n", &ppatch_length);
        fscanf(fp,"%*s %*s%d\n", &pboundary);
        fscanf(fp,"%*s %lf\n", &(pmove_size[X]));
        fscanf(fp,"%*s %lf\n", &(pmove_size[Y]));
        fscanf(fp,"%*s %lf\n", &(pmove_size[Z]));
        fscanf(fp,"%*s %lf\n", &(pmove_size[ROTATION]));
        fscanf(fp,"%*s %lf\n", &(pmove_size[VOLUME+X]));
        fscanf(fp,"%*s %lf\n", &(pmove_size[VOLUME+Y]));
        fscanf(fp,"%*s %lf\n", &(pmove_size[CLUSTER]));
        fscanf(fp,"%*s %lf\n", &(pmove_size[PATCH]));
        fscanf(fp,"%*s %s ",prefix);
        strcpy(pinfile, prefix); strcat(pinfile, ".in");
        strcpy(poutfile, prefix); strcat(poutfile, ".out");
        strcpy(plogfile, prefix); strcat(plogfile, ".log");
        strcpy(ptrjfile, prefix); strcat(ptrjfile, ".trj");
        strcpy(pcrdfile, prefix); strcat(pcrdfile, ".crd");
        strcpy(psffile, prefix); strcat(psffile, ".psf");
        strcpy(pdbfile, prefix); strcat(pdbfile, ".pdb");
		fscanf(fp, "%*s %d\n", &pdo_run_preset); 
		fscanf(fp, "%*s %d\n", &pgenerate_preset); 
		if ((pdo_run_preset == YES) && (pgenerate_preset == YES))
		{
			printf("Cannot simultaneously run and generate preset, not generating. \n"); 
			pgenerate_preset = NO; 
		}
	if (feof(fp))
		printf("Reshaping method not set. Default is POSSIBLY_CONSTRAINED.\n"); 
	else
		fscanf(fp, "%*s %*s %*s %d ", &preshaping_method); 
    }
    fclose(fp);



}


void get_model_prefs()
/*reads in the preferences from "prefile"*/
{
    FILE *fp;
    extern int pcore;
    extern double pdiameter;
    extern potential pmodel[];
    extern double ptorsion;

    int i;

    if ((fp = fopen("model.pref","r"))==NULL)
    {puts("cannot open preferences\n");
    }
    else
    {
        fscanf(fp,"%*s %lf\n", &pdiameter);
        fscanf(fp,"%*s %lf\n",&plength);
        fscanf(fp,"%*s %d\n", &pcore);

        fscanf(fp,"%*s %*s %*s %*s %d\n", &(pmodel[HINTERFACE].num_terms));
        fscanf(fp,"%*s %*s %*s %*s %lf\n", &(pmodel[HINTERFACE].range));
        for (i = 0; i<pmodel[HINTERFACE].num_terms; i++)
           { fscanf(fp,"%lf %lf\n", &(pmodel[HINTERFACE].coefficients[i][0]), &(pmodel[HINTERFACE].coefficients[i][1]));
                 pmodel[HINTERFACE].coefficients[i][2] = poly_energy_noshift(pow(pmodel[HINTERFACE].range,2),pmodel[HINTERFACE],pdiameter);
}
        fscanf(fp,"%*s %*s %*s %*s %d\n", &(pmodel[HTAIL].num_terms));
        fscanf(fp,"%*s %*s %*s %*s %lf\n", &(pmodel[HTAIL].range));
        for (i = 0; i<pmodel[HTAIL].num_terms; i++)
        {fscanf(fp,"%lf %lf\n", &(pmodel[HTAIL].coefficients[i][0]), &(pmodel[HTAIL].coefficients[i][1]));
           pmodel[HTAIL].coefficients[i][2] = poly_energy_noshift(pow(pmodel[HTAIL].range,2),pmodel[HTAIL],pdiameter);

        }


        fscanf(fp,"%*s %*s %*s %*s %d\n", &(pmodel[HSOFTCORE].num_terms));
        fscanf(fp,"%*s %*s %*s %*s %lf\n", &(pmodel[HSOFTCORE].range));
        for (i = 0; i<pmodel[HSOFTCORE].num_terms; i++)
        {fscanf(fp,"%lf %lf\n", &(pmodel[HSOFTCORE].coefficients[i][0]), &(pmodel[HSOFTCORE].coefficients[i][1]));
           pmodel[HSOFTCORE].coefficients[i][2] = poly_energy_noshift(pow(pmodel[HSOFTCORE].range,2),pmodel[HSOFTCORE],pdiameter);}

        fscanf(fp,"%*s %lf\n", &(ptorsion));

        fclose(fp);
    }
	make_model_table();

}

void get_domain_prefs()
{
	extern double pdomaindx, pdomaindr, pumbrella_coeff, pumbrella_center, pmove_size[]; 
	FILE *fp; 
	if ((fp = fopen("domain.pref","r"))==NULL)
		{
			printf("No domain preferences read.\n"); 
			pdomaindx = pmove_size[X]; 
			pdomaindr = pmove_size[ROTATION]; 
			pumbrella_coeff = 0.0; 
			pumbrella_center = 0.0; 
			return; 
		}
	fscanf(fp,"%*s %lf\n", &pdomaindx); 
	fscanf(fp,"%*s %lf\n", &pdomaindr); 
	fscanf(fp,"%*s %*s %lf\n", &pumbrella_coeff);
	fscanf(fp,"%*s %*s %lf\n", &pumbrella_center);  
	fscanf(fp,"%*s %*s %d\n", &pdomain_fixed); 
	fscanf(fp,"%*s %*s %d\n", &pdomain_constrain); 
	fclose(fp);
}

void write_domain_prefs()
{
	extern double pdomaindx, pdomaindr, pumbrella_coeff, pumbrella_center;
	FILE *fp;  
	if ((fp = fopen("domain.pref","w"))==NULL)
		{
			printf("Cannot open domain preferences for writing.\n"); 
			return; 
		}
	fprintf(fp,"DX: %lf\n", pdomaindx);
	fprintf(fp,"DR: %lf\n", pdomaindr);
	fprintf(fp,"UMBRELLA COEFFICIENT: %lf\n", pumbrella_coeff);
	fprintf(fp,"UMBRELLA CENTER: %lf\n", pumbrella_center);
	fprintf(fp,"FIXED DOMAIN: %d\n", pdomain_fixed);
	fprintf(fp,"CONSTRAIN ORIENTATION: %d\n", pdomain_constrain);  
	fclose(fp);
}

void get_user_prefs()
{
	FILE *fp; 
	int mol1, bead1, mol2, bead2, type; 
	double coefficient, exponent, r0, range; 
	int i; 
	double a,b,c,d,e,f; 
	if ((fp = fopen("user.pref","r"))==NULL)
		{
			printf("No user preferences read.\n");
			pnumber_user_interactions = 0; 
			puser_potentials = NULL; 
			return; 
		}
	fscanf(fp,"%*s %d\n", &pnumber_user_interactions);
	printf("number user interactions specified: %d\n", pnumber_user_interactions); 
	if (pnumber_user_interactions == 0)
		{puser_potentials = NULL;}
	else 
	{
		puser_potentials = (upotential *)(malloc(pnumber_user_interactions * sizeof(upotential))); 
		fscanf(fp, "%*s %*s %*s %*s %*s %*s %*s %*s %*s \n"); 
		for (i = 0; i < pnumber_user_interactions; i++)
		{
			fscanf(fp, "%d %d %d %d %d %lf %lf %lf %lf\n", &mol1, &bead1, &mol2, &bead2, &type, &coefficient, &exponent, &r0, &range); 
			puser_potentials[i].ci1.m = mol1; 
			puser_potentials[i].ci2.m = mol2; 
			puser_potentials[i].ci1.b = bead1; 
			puser_potentials[i].ci2.b = bead2; 
			puser_potentials[i].type = type; 
			puser_potentials[i].coefficient = coefficient;
			puser_potentials[i].exponent = exponent; 
			puser_potentials[i].r0 = r0; 
			puser_potentials[i].range = range; 
		}
	}		
	
	fscanf(fp,"%*s %d\n", &pnumber_user_constraints);
	if (pnumber_user_constraints == 0)
		{puser_constraints = NULL;return; }
	
	puser_constraints = (uconstraint *)(malloc(pnumber_user_constraints * sizeof(uconstraint))); 
	fscanf(fp, "%*s %*s %*s %*s %*s %*s %*s %*s %*s \n"); 
	for (i = 0; i < pnumber_user_constraints;i++)
	{
		fscanf(fp, "%d %d %d %lf %lf %lf %lf %lf %lf\n", &mol1, &bead1, &type, &a, &b, &c, &d, &e, &f); 
		puser_constraints[i].ci.m = mol1; 
		puser_constraints[i].ci.b = bead1; 
		puser_constraints[i].type = type;
		printf("read one %d type\n", type); 
		puser_constraints[i].point[X] = a; 
		puser_constraints[i].point[Y] = b; 
		puser_constraints[i].point[Z] = c; 			
		if (type = FIXED_ANGLE)
		{
			puser_constraints[i].arrow[X] = d; 
			puser_constraints[i].arrow[Y] = e; 
			puser_constraints[i].arrow[Z] = f;
			vunit(puser_constraints[i].arrow, puser_constraints[i].arrow);  			
		}		
	}
	fclose(fp);
}

int update_prefs()
{
    extern long pnum_steps, ptotal_steps;
    extern double pkt;
    long temp_num_steps, temp_total_steps;
    FILE *fp;
    double temp_kt;
    char *updatefile = "up.pref";

    if ((fp = fopen(updatefile,"r"))==NULL)
        return NO_UPDATE;
    else
    {
        fclose(fp);
        temp_num_steps=pnum_steps;
        temp_total_steps = ptotal_steps;
        temp_kt = pkt;
        get_sim_prefs();
        get_model_prefs();
        remove(updatefile);
        pnum_steps = temp_num_steps;
        ptotal_steps = temp_total_steps;
      		pkt = temp_kt;
                return UPDATE;
    }

}
