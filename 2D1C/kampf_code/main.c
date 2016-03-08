/**********************************************************
 * Program that generates capillary pattern according to the 
 * statistics of several deseases
 * this is after phyica A Vol 369 (2006) pp 599-611
 * for details read this paper
 */
 
 
 
/* Standard includes */
# include<stdio.h>
# include<stdlib.h>
# include<math.h>
# include<time.h>

/* eigene includes */

# include "simstruct.h"




int main(void)
{
	int c1,c2;
	
	/* declaration of the File output pointers */
	FILE *capillarydistribution;
	FILE *qualitycheck;  
	FILE *fileInitialization;
	
	
	simstruct simdata;
	simstruct *simdatap;
	
	simdatap = &simdata;
	
	
	
	/* the following is done to "clear the according files from their 
	 * content. These files are opened with the append argument by functions in this
	 * programm and thus it has to be ensured that there is no content in these files
	 * if they exist from preceding runnings of the program since they are thought to
	 * used in the statistical analysis of the generated capillary distribution
	 */
	 
	/* capillary distribution file */
	fileInitialization = fopen("capillarydistributions.dat","w");
	fclose(fileInitialization);

	/* distance distribution file */
	fileInitialization = fopen("distancedistribution.dat","w");
	fclose(fileInitialization);
	
	/* pair correlation function file */
	fileInitialization = fopen("pair_correlation_function.dat","w");
	fclose(fileInitialization);
	
	/* additional parameters for simulation check */
	fileInitialization = fopen("simulation_check_parameters.dat","w");
	fclose(fileInitialization);

		
	
	/* initializing the structure */
	initialize_structure(simdatap);
	
		
	/* initializing the capillary distribution */
	initialize_capillary_distribution(simdatap);

	/* calculation of the energy of the current
	 * capillary distribution */
	calculate_variable_energy_part(simdatap);

	
	/* Loop ueber die gewuenschten Konfigurationen
	 */
	for(c1 = 0; c1<simdatap->number_of_configurations; c1++)
	{
		simdatap->currConfig = c1;
		for(c2 = 0; c2 < simdatap->number_of_intermed_steps; c2++)
		{
			calculate_next_intermed_step(simdatap);		
		}


		/* Printing the data of the current capillary distribution
		 * done this way, thus it is possible to observe intermediate results
		 * and the file only open when data are written
		 */

		/* initialization of the file output pointers */
		capillarydistribution = fopen("capillarydistributions.dat", "a");
		qualitycheck = fopen("simulation_check_parameters.dat", "a");
	
		for(c2 = 0; c2 < simdatap->capillaryNumber; c2++)
		{
			fprintf(capillarydistribution,"%lf  %lf   %d   %d\n",  simdatap->x[c2], simdatap->y[c2], c1+1, c2);
		}

		fprintf(capillarydistribution,"\n");
	
		/* printing additional parameters to judge on the quality of the simulation data */
		fprintf(qualitycheck,"%d %lf  %lf  %lf  %e  %d\n", c1,  simdatap->old_energy, simdatap->cur_energy, 
		simdatap->cur_energy - simdatap->old_energy, 
		fmin(1.0,exp(-4.0*(simdatap->cur_energy - simdatap->old_energy)*simdatap->coupling_parameter)), simdatap->currCapillary);


		
		/* closing the outputfiles */
		fclose(capillarydistribution);	
		fclose(qualitycheck);

		/* Analysing the current capillary distribution for statistical purpose */
		nearest_neighbour_dist(simdatap);
		pair_correlation_function(simdatap);
		
		printf("%d configs of %d\n",c1, simdatap->number_of_configurations);	
	}
	
	
	
	return 0;
} 
