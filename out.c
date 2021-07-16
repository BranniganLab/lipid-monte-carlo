#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "structures.h"
#include "out.h"
#include "vectorops.h"
#include "pref.h"
#include "domain.h"
#include "constraints.h"

/****************************************************************/
void write_version()
{
	extern long ptotal_steps;
	FILE *fp;
 	fp = fopen("run_info","a");
  	fprintf(fp,"step %ld\n", ptotal_steps);
	fprintf(fp, "version:5.1 read in topoLogy, separate 'types' stored instead of topology being stored in each molecule");
	fprintf(fp, "version:5.0 multiple architectures"); 
	fprintf(fp, "version:4.5 ghosts + fourier moves again"); 
	fprintf(fp, "version:4.4 allows user constraints\n"); 
	fprintf(fp, "version: 4.3 reads in user preferences, reads bead type from mc.model, reads domain preferences for umbrella sampling\n"); 
   	fprintf(fp, "version: 4.2 allows more than one inclusion \n");
   	fprintf(fp, "version: 4.1 allows inclusions. \n");
     fprintf(fp, "version:4.0 overhauled, lattice structure\n");
   	fprintf(fp, "version: 3.0 each molecule has own bending potential\n");
   	fprintf(fp,"version: 2.6 faster algorithm for neighbor list update. Fixed small bug for systems without cell lists\n");
   	fprintf(fp, "version:2.5.2 uses proper boltzmann weight for non zero tension\n");
  	fprintf(fp,"version: 2.5.1\n");
   // if (USE_CELL_LIST == YES)
	//	fprintf(fp, "neighbor list: %d cell list: %d \n", USE_LIST, NUM_CELLS);
//    else
	//	fprintf(fp, "neighbor list: %d cell list: %d \n", USE_LIST,USE_CELL_LIST);
	fprintf(fp, "does not currently rebalance tree.\n");
    fprintf(fp,"2.5.1 checks for abort move\n"); 
	fprintf(fp, "2.5 uses fourier modes, long range jump capability\n");
	fprintf(fp,"2.4 uses chain list\n");
   	fprintf(fp,"2.3 has center-tail attractions\n");
	
    fclose(fp);
}
/*****************************************************************/
void write_interval(particle *the_membrane, mparticle *the_types,  double total_energy, int total_accept[],
                    double move_size[], long i, vector box_length)
{
    extern char poutfile[30], plogfile[30];
    extern long pnum_particles;
    extern long ptotal_steps;
    extern vector pbox_length;
    FILE *fp;
    double total_points;
    int num_clusters;
    int j;

	//    if (pnum_particles < 500)
	//   num_clusters =   smake_cluster(the_membrane, pbox_length) ;
	//  else
	// num_clusters = 0;
	
	//    printf("num clusters: %d\n", num_clusters);
    //sort_clusters(*the_membrane, num_clusters);
    write_trajectory("trj.out", the_membrane, the_types,  total_energy,  move_size,  i) ;
    WriteXYZ(poutfile, the_membrane, the_types);
    WriteLog(plogfile, the_membrane, total_energy, total_accept, move_size, i);
    write_frame("mc.pdb", the_membrane,the_types, "w");
    write_simple_psf("mc.psf", the_membrane, the_types);

    if (ptotal_steps == 0)
	  {write_crd("mc.crd", the_membrane,"w");
        printf("overwriting mc.crd\n");
        fp = fopen("duda","w");
        fclose(fp);
        fp = fopen("edist", "w");
        fclose(fp);
		fp = fopen("mc.user", "w"); 
		fclose(fp); 
	  }
    else
        write_crd("mc.crd", the_membrane,"a");
	print_constraints(the_membrane, the_types, box_length); 
   	write_domains("mc.domains", the_membrane);
    write_sim_prefs();
    write_model_prefs();
   // write_domain_prefs();
	write_user_distances(the_membrane, box_length); 
	
	//    if (use_inclusions == YES)
	//  {
	
    //    write_isim_prefs();
    //    write_imodel_prefs();      }
    printf("Iteration Energy at step %ld: %12.6f\n", i,total_energy);
}

/********************************************************/
void write_domains(char *outfile, particle *the_membrane)
{
	int i, j;
	extern long pnum_particles, ptotal_steps;
	double radius;
	extern double pumbrella_coeff, pumbrella_center;
	double total = 0;
 	extern vector pbox_length;
	domain *domain_list[MAX_DOMAINS];
 	FILE *fp;

    	for (i = 0; i < MAX_DOMAINS; i++)
     		domain_list[i] = NULL;
    if (ptotal_steps < 2)
        fp = fopen(outfile,"w");
    else
        fp = fopen(outfile,"a");

     if (fp == NULL)
     	{
          printf("error. \n"); return; }
    for (i = 0; i<pnum_particles; i++)
		if (the_membrane[i].domain_index != 0)
  			{if (the_membrane[i].domainptr != NULL)
				domain_list[the_membrane[i].domain_index] = the_membrane[i].domainptr;
    			else
       			printf("non zero domain index has null domain pointer!\n");
          	}
// printf("here.\n");
   	for (i = 0; i < MAX_DOMAINS; i++)
    	for (j= i+1; j < MAX_DOMAINS; j++)
			if ((domain_list[i]!=NULL) && (domain_list[j]!=NULL))
				{
					radius = sqrt(vdistance_periodic(domain_list[i]->center, domain_list[j]->center,pbox_length));
					fprintf(fp, "%lf\n", radius);      				
				}
     fclose(fp);
    	
}



/*****************************************************************/
void write_trajectory(char *outfile,particle *the_membrane,mparticle *the_types,
                      double total_energy, double move_size[], long step)
/*append current ensemble to a file containing ensembles in XYZ format*/
/*specify step number, energy, box size, move sizes, to enable reconstruction*/
{
    extern double pvol_moves,pcluster_moves;
    extern long pnum_particles, ptotal_steps;
    extern vector pbox_length;
	
    FILE *fp;
    int i;  int j;
    if (ptotal_steps == 1)
        fp = fopen(outfile,"w");
    else
        fp = fopen(outfile,"a");
    if (fp == NULL)
	  {printf("cannot open trajectory file:\n");}
    else{
        fprintf(fp, "%8ld  %8.6lf  ",(long) ((double)ptotal_steps/((double)pnum_particles)),
                total_energy);
        fprintf(fp, "%8.6lf %8.6lf %8.6lf\n", pbox_length[X], pbox_length[Y], pbox_length[Z]);
        fprintf(fp,"%6.4lf %6.4lf %6.4lf %6.4lf %6.4lf %6.4lf\n", move_size[0],move_size[1],
                move_size[2],move_size[ROTATION],move_size[VOLUME],move_size[CLUSTER]);
        fprintf(fp, "%6.4lf %6.4lf\n", pvol_moves, pcluster_moves);
        for (i = 0; i<pnum_particles; i++)
		  {
            fprintf(fp,"%10.5f  %10.5f %10.5f ",the_membrane[i].chain[1].position[X],
                    the_membrane[i].chain[1].position[Y], the_membrane[i].chain[1].position[Z]);
            fprintf(fp,"%10.5f  %10.5f %10.5f ",the_membrane[i].chain[0].position[X],
                    the_membrane[i].chain[0].position[Y], the_membrane[i].chain[0].position[Z]);
            for (j = 2; j < the_membrane[i].chain_length; j++)
                fprintf(fp,"%10.5f  %10.5f %10.5f ",the_membrane[i].chain[j].position[X],
                        the_membrane[i].chain[j].position[Y], the_membrane[i].chain[j].position[Z]);
            fprintf(fp,"%10.5f", the_types[the_membrane[i].model_index].chain[1].radius * 2);
            fprintf(fp,"%10.5f ", the_types[the_membrane[i].model_index].length);
            fprintf(fp,"\n");
		  }
		
        fprintf(fp, "\n");
        fclose(fp);
    }
}


/**************************************************************************************/
void WriteXYZ(char *outfile,   particle *the_membrane, mparticle *the_types)
/*write out the current ensemble to a file that
can then be used as input for another run. */
{
    extern long pnum_particles;
    FILE *fp;
    int i,j;
    fp = fopen(outfile,"w");
    if (fp == NULL)
        printf("cannot open outfile:%s\n", outfile);
    else
	  {
        for (i = 0; i<pnum_particles; i++)
		  {
			
            fprintf(fp,"%10.5f  %10.5f %10.5f ",the_membrane[i].chain[1].position[X],
                    the_membrane[i].chain[1].position[Y], the_membrane[i].chain[1].position[Z]);
            fprintf(fp,"%10.5f  %10.5f %10.5f ",the_membrane[i].chain[0].position[X],
                    the_membrane[i].chain[0].position[Y], the_membrane[i].chain[0].position[Z]);
            for (j = 2; j < the_membrane[i].chain_length; j++)
                fprintf(fp,"%10.5f  %10.5f  %10.5f ",the_membrane[i].chain[j].position[X],
                        the_membrane[i].chain[j].position[Y], the_membrane[i].chain[j].position[Z]);
            fprintf(fp,"%10.5f ",the_types[the_membrane[i].model_index].chain[1].radius * 2.0);
            fprintf(fp,"%10.5f ",the_types[the_membrane[i].model_index].length);
            fprintf(fp,"\n");
		  }
        fclose(fp);
	  }
}
/***********************************************************************************/
void WriteLog(char *outfile,   particle *the_membrane, double total_energy,
              int total_accept[], double move_size[], long step)
/*write out the simulation's progress, some runtime vars*/
{
    extern long pnum_particles,ptotal_steps;
    extern vector pbox_length;
    extern double pkt;
    FILE *fp;
    int i;
	
	
	
    if (ptotal_steps == 1)
        fp = fopen(outfile,"w");
    else
        fp = fopen(outfile,"a");
    if (fp == NULL)
	  {printf("cannot open logfile:%s\n", outfile);}
    else{
        fprintf(fp, "Step: %8ld Energy %8.6lf   Length: %8.6lf  KT: %8.6lf",(long) ((double)ptotal_steps/
                                                                                    ((double)pnum_particles)), total_energy, pbox_length[X], pkt );
        fprintf(fp,"dX:%6.4lf ", move_size[0]);
        fprintf(fp,"dY:%6.4lf ", move_size[1]);
        fprintf(fp,"dZ:%6.4lf ", move_size[2]);
        fprintf(fp,"dR:%6.4lf ", move_size[ROTATION]);
        fprintf(fp,"dA:%6.4lf ", move_size[VOLUME]);
        fprintf(fp, "\n");
        fclose(fp);
    }
}



/***************************************************************/
void write_crd(char *outfile,   particle *the_membrane, char *code)
/*append to a crd file for viewing in vmd*/
{
    extern long pnum_particles;
    extern double ihydro_fraction;
    extern vector pbox_length;
    extern int *swaps;
    extern double *initlengths;
    FILE *fp;
    int i, count,index,j;
    double x, y, z, a, b, c,l, tempdouble;
    vector head, tail, center,temp_vect;
	
	

	
    if (strcmp(code,"a") == 0)
	  { fp = fopen(outfile,"r");
        if (fp != NULL)
            fclose(fp);
        else
            code ="w"; }
	
	
    if ((fp = fopen(outfile,code))==NULL)
	  {puts("cannot open crdfile\n ");
	  }
    else{
        if (strcmp(code,"w") == 0)
		  {fprintf(fp,"TITLE\n");}
		
        count = 0;
        for (i = 0; i<pnum_particles; i++)
		  {
            index = i;
            x = the_membrane[index].chain[1].position[X];
            y = the_membrane[index].chain[1].position[Y];
            z = the_membrane[index].chain[1].position[Z];
			
            vcopy(head, the_membrane[index].chain[0].position);
			
			
            /****************tail *********************/
            for (j =0; j <  the_membrane[i].chain_length; j++)
			  { vcopy(tail, the_membrane[index].chain[j].position);
                fprintf(fp,"%8.3f",tail[X]);
                count++;
                if (count  ==  8)
				  {fprintf(fp,"\n"); count = 0; }
				
                fprintf(fp,"%8.3f",tail[Y]);
                count++;
                if (count  ==  8)
				  {fprintf(fp,"\n"); count = 0; }
				
                fprintf(fp,"%8.3f",tail[Z]);
                count++;
                if (count  ==  8)
				  {fprintf(fp,"\n"); count = 0; } }
			
		  }
        fclose(fp);
    }
}

/***************************************************************/
void write_frame(char *outfile,   particle *the_membrane, mparticle *the_types,  char *code)
/*write a single pdb file for vmd viewing, containing
the current ensemble*/
{
    extern long pnum_particles;
    extern vector pbox_length;
	
    FILE *fp;
    int i,j;
    double x, y, z, a, b, c;
    int chainlength;
    vector head, tail, center ;
	int count = 1; 
	
	
    code = "w";
    if ((fp = fopen(outfile,code))==NULL)
	  {puts("cannot open pdbfile\n ");
	  }
    else{
        if (code == "w") {
			
            fprintf(fp,"COMPND    Monte Carlo Simulation of Lipids\n");
            fprintf(fp,"AUTHOR    Created by Grace Brannigan at UCSB\n"); }
        for (i = 0; i<pnum_particles; i++)
		  {
            x = the_membrane[i].chain[1].position[X];
            y = the_membrane[i].chain[1].position[Y];
            z = the_membrane[i].chain[1].position[Z] ;
            chainlength = the_membrane[i].chain_length;
            vcopy(head, the_membrane[i].chain[0].position);
            for (j = 0; j <the_membrane[i].chain_length; j++)
			  {vcopy(tail, the_membrane[i].chain[j].position);
			  if (the_membrane[i].chain[j].typeptr->type == TAIL)
                fprintf(fp,"HETATM %4d  C        %4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %2d\n",
                        count,i+1,tail[X], tail[Y], tail[Z], 0.0, 0.0,
                        (int)the_types[the_membrane[i].model_index].cbend + 10);
			  else if (the_membrane[i].chain[j].typeptr->type == INTERFACE)
                fprintf(fp,"HETATM %4d  G        %4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %2d\n",
                        count,i+1,tail[X], tail[Y], tail[Z], 0.0, 0.0,
                        (int)the_types[the_membrane[i].model_index].cbend + 10);
			  else
				fprintf(fp,"HETATM %4d  N        %4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %2d\n",
                        count,i+1,tail[X], tail[Y], tail[Z], 0.0, 0.0,
                        (int)the_types[the_membrane[i].model_index].cbend + 10);
				count++; }
			
		  }
        fclose(fp);
    }
}

/*****************************************************************/
void write_simple_psf(char *outfile,   particle *the_membrane, mparticle *the_types)
/*write the psf file for vmd- here the particles can not cross
over the boundary*/
{
    extern long pnum_particles;
    char *type;
    int count=1;
    int chainlength;
    int chain1start, chain2start, chain1stop, ibindex, ib, numchains; 
    FILE *fp;
    int i, j;
	int total_particles = 0; 
	int total_bonds = 0; 
	int count2 = 1; 
	
	for (i = 0; i < pnum_particles; i++)
		{
			total_particles = total_particles + the_membrane[i].chain_length;
			total_bonds = total_bonds + the_membrane[i].chain_length - 1; 
		} 
   // chainlength = the_membrane[0].chain_length;
    if ((fp = fopen(outfile,"w"))==NULL)
	  {puts("cannot open psffile\n ");}
    else
	  {
        fprintf(fp, "PSF\n");
        fprintf(fp, "%8d %s\n", 1,"!NTITLE:");
        fprintf(fp, " REMARKS FILENAME = 'mc.psf'\n");
        fprintf(fp, " REMARKS conect information for mc.pdb\n");
        fprintf(fp, "\n%8ld %s\n", total_particles,"!NATOM:");
        for (i = 0; i<pnum_particles; i++)
		  {
		  	
			chainlength = the_membrane[i].chain_length;
		  for (j = 0; j<(chainlength); j++)
			{
			  if (the_membrane[i].chain[j].typeptr->type == HEAD)
				  type = "N";
			  else if (the_membrane[i].chain[j].typeptr->type == INTERFACE)
				  type = "G";
			  else
				  type = "C";
			  /***/
			  fprintf(fp, "%8d             %4s    %13.10f%13.10f\n", count,
					  type, the_types[the_membrane[i].model_index].cbend*10 + 10,(double)
					  the_types[the_membrane[i].model_index].cbend + 10.0
					  );
				count++; 
			  //hack on the +10, it doesn't seem to read charges less than ten,
			  //probably an outputting problem.
			}
		  }
		count = 0; 
		count2 = 1; 
        fprintf(fp, "\n%8ld %s\n",total_bonds,"!NBONDS:bonds");
        for (i= 0; i<pnum_particles; i++)
		  {
		  	chain1start = the_types[the_membrane[i].model_index].chain1start; 
			chain2start = the_types[the_membrane[i].model_index].chain2start; 

			ibindex = the_types[the_membrane[i].model_index].ibindex; 
			chainlength = the_membrane[i].chain_length;
			if (chain2start < chainlength)
			{
				numchains = 2; 
				chain1stop = chain2start; }
			else
			{
				numchains = 1; 
				chain1stop = chainlength; 
			}	
	    		for (j = 0; j < chain1stop-1; j++)
	    		{
				fprintf(fp, "%8d%8d", count2, count2+1);
				if (j == ibindex) ib = count2; 
				count2++; 
				count++;  
				if (count == 4)
					{fprintf(fp,"\n"); count = 0;  }
			}
	    		if (numchains == 2)
	    		{
				count2++; 
				fprintf(fp, "%8d%8d", ib, count2);
				count++;  
				if (count == 4)
					{fprintf(fp,"\n"); count = 0;  }
				for (j = chain1stop; j < chainlength-1; j++)
	    			{
					fprintf(fp, "%8d%8d", count2, count2+1);
					count2++; 
					count++;  
					if (count == 4)
						{fprintf(fp,"\n"); count = 0;  }
				}
			}
			count2++; 
		  }
        fprintf(fp, "\n%8d %s\n", 0, "!NTHETA:angles");
        fprintf(fp, "\n%8d %s\n", 0, "!NPHI:dihedrals");
        fprintf(fp, "\n%8d %s\n", 0, "!NIMPHI:impropers");
		
        fclose(fp);
	  }
}
/*************************************************/
int write_user_distances(particle *the_membrane, vector box_length)
{
	FILE *fp; 
	int i; 
	vector radius; 
	extern int pnumber_user_interactions; 
	extern upotential *puser_potentials; 
	fp = fopen("mc.user","a"); 
	
	for (i = 0; i < pnumber_user_interactions; i++)
		{
			//printf("ci1: %d, ci2:%d\n", puser_potentials[i].ci1.m, puser_potentials[i].ci2.m); 
			vsub_periodic(radius,the_membrane[puser_potentials[i].ci1.m].chain[puser_potentials[i].ci1.b].position, 
									the_membrane[puser_potentials[i].ci2.m].chain[puser_potentials[i].ci2.b].position, box_length);
			radius[Z] = 0; 
			fprintf(fp, "%lf ", sqrt(vsquare(radius))); 
		}
	fprintf(fp, "\n"); 
	fclose(fp); 
}
