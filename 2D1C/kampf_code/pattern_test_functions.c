/***********************************************************
 * in this file the functions to calculate the distribution
 * nearest neighbour distances as well as the distance correlation
 * function are analysed
 */
 
 
/* Standard includes */
# include<stdio.h>
# include<stdlib.h>
# include<math.h>
# include<time.h>

/* eigene includes */

# include "simstruct.h"



/* This function creates a histogram on the distribution of the distance to the nearest neigbbour
 * The resolution of the hystogram is chosen by distanceResolution which gives the number of subsections
 * for the natural distance. The natural distance is given by the distance between two capillaries.
 */

int nearest_neighbour_dist(simstruct *simdatap)
{
	double mindist[DEF_VAL_CAPILLARYNUMBER];
	double cur_distance;
	int c1,c2;
	double largest_possible_distance;
	double normalizing_distance;		/* radius of the mean support or dephasing or whatever it is called volume */
	int index_number;
	FILE *ofp_distHyst;
 	double *distanceHyst;
	int index;
	double normalizing_constant;	
 	
	/* opening the outputfile, the flag a is for appending */
	ofp_distHyst = fopen("distancedistribution.dat","a");
	
	
	
	largest_possible_distance = sqrt(simdatap->lattice_vector_ax*simdatap->lattice_vector_ax + simdatap->lattice_vector_by*simdatap->lattice_vector_by);
	
	/* loop over the particle the distance is calculated to */
	for(c1 = 0; c1 < simdatap->	capillaryNumber; c1++)
	{
		/* initializing the array with redicules large values (larger than possible)
		 */	
		mindist[c1] = largest_possible_distance;

		simdatap->currCapillary = c1;
		
		for(c2 = 0; c2 < simdatap->capillaryNumber; c2++)
		{
			/* test if both capillaries are the same */
			if(c2 == c1)
			{
				/* distance to itself is zero and neglected */
			}
			else
			{
				simdatap->second_capillary_for_distance_calculation = c2;
			
				cur_distance = sqrt(calculate_distance_squared(simdatap));
				
				if(cur_distance < mindist[c1])
				{
					mindist[c1] = cur_distance;
				}
				else{};
			}
		
		}
	}
	
	/* vector mindist now contains the information about the clostest distance of the next capillary */
	normalizing_distance = sqrt(simdatap->simArea/simdatap->capillaryNumber/M_PI);
	
	/* calculation of the number of vector indizes needed */
	index_number = (int)(ceil(largest_possible_distance/normalizing_distance * simdatap->distanceResolution ));

	distanceHyst = (double *) malloc(index_number * sizeof(double));
	
	/* test if the memory is available */
	if(distanceHyst == NULL)
	{
		/* not thus exiting the programm */
		printf("Memory could not be reserved for the hystogramm!\n");
		exit(0);
	}
	else{};

	/* initializing the hystogram */
	for(c1 = 0; c1 < index_number; c1++)
	{
		distanceHyst[c1] = 0.0;
	}
	
	/* creating the hystogram */
	for(c1 = 0; c1 < simdatap->capillaryNumber; c1++)
	{
		index = (int)round(mindist[c1]/largest_possible_distance*index_number);

		if((index < index_number) && (index >= 0))
		{
			distanceHyst[index] += 1.0;
		}
		else
		{
			printf("ohoh, das sollte nicht passieren!!!!!!!!!!\n");
			printf("index = %d; capillary = %d; index_number = %d; min. distance = %lf; max. poss. distance = %lf\n",index, c1, index_number,mindist[c1],largest_possible_distance);
			exit(0);
		};
	}
	
	/* determining the normalization constant */
	normalizing_constant = 0.0;
	for(c1 = 0; c1 < index_number; c1++)
	{
		normalizing_constant += distanceHyst[c1];
	}
	
	/* normalizing the histogram */
	for(c1 = 0; c1 < index_number; c1++)
	{
		distanceHyst[c1] /= normalizing_constant;
	}
	

		
	/* plotting the hystogramm in the hystogram file */
	for(c1 = 0; c1 < index_number; c1++)
	{
		fprintf(ofp_distHyst," %lf   %lf  %d\n",(c1 + 0.5)/simdatap->distanceResolution ,distanceHyst[c1], simdatap->currConfig);
	}

	fprintf(ofp_distHyst,"\n");
	
	/* closing the filepointer */
	fclose(ofp_distHyst);
 	
	/* Freeing the memory of the hystogramm */
	free(distanceHyst);
	
	return 0;
}









/* function that calculates the pair correlation function which is nothing more than
 * the distribution of the distances between different capillaries. Thus it is 
 * given in the form of a hystogram.
 * The resolution of the hystogram is choosen by distanceResolution which gives the number of subsections
 * for the natural distance. The natural distance is given by the distance between two capillaries.
 */



int pair_correlation_function(simstruct *simdatap)
{

	FILE *ofp_pair_correlation;
	double *pair_correlation_function;
	
	int c1, c2;
	int index_number;
	double largest_possible_distance;
	double normalizing_distance;
	double cur_distance;
	int index;
	double normalizing_constant;
	


	/* opening the outputfile, the flag a is for appending */
	ofp_pair_correlation = fopen("pair_correlation_function.dat","a");

	/* calculating the size of the hystogram for the pair correlation function
	 */
	largest_possible_distance = sqrt(simdatap->lattice_vector_ax*simdatap->lattice_vector_ax + simdatap->lattice_vector_by*simdatap->lattice_vector_by);

	/* vector mindist now contains the information about the clostest distance of the next capillary */
	normalizing_distance = sqrt(simdatap->simArea/simdatap->capillaryNumber/M_PI);
	
	/* calculation of the number of vector indizes needed */
	index_number = (int)(ceil(largest_possible_distance/normalizing_distance * simdatap->distanceResolution ));
	

	/* allocating the memory for the hystogramm */
	pair_correlation_function = (double *) malloc(index_number * sizeof(double));
	
	
	/* loop over the current capillary
	 */
	for(c1 = 0; c1 < simdatap->capillaryNumber; c1++)
	{
		simdatap->currCapillary = c1;
		
		/* loop over the second capillary 
		 * This loop is only over the capillaries left to avoid 
		 * the distance zero if c1 == c2
		 * Furthermore this way, doublecounting distances due to the 
		 * symmetry of the distance behavior ( dist(a,b) == distance(b,z) is avoided
		 */
		for(c2 = c1 + 1; c2 < simdatap->capillaryNumber; c2++)
		{
			simdatap->second_capillary_for_distance_calculation = c2;
			
			cur_distance = sqrt(calculate_distance_squared(simdatap));
			
			index = (int)round(cur_distance/largest_possible_distance*index_number);

			if((index < index_number) && (index >= 0))
			{
				pair_correlation_function[index] += 1.0;
			}
			else
			{
				printf("ohoh, das sollte nicht passieren!!!!!!!!!!\n");
				printf("index = %d; capillary = %d; index_number = %d; min. distance = %lf; max. poss. distance = %lf\n",index, c1, index_number, cur_distance, largest_possible_distance);
				exit(0);
			};
		}
	}
	
	/* determining the normalization constant */
	normalizing_constant = simdatap->capillaryNumber * simdatap->capillaryNumber * M_PI;  /* the two from the volume is missing since in the computation the 
																						   * off the pair correlation function the symmetry in the expression
																						   * is used und thus only the half of the functions is actually evaluated
																						   * this is compensated by neglecting the two in the volume
																						   */
	normalizing_constant *= normalizing_distance * normalizing_distance;
	normalizing_constant /= simdatap->simArea * simdatap->distanceResolution * simdatap->distanceResolution;
	
	
	/* normalizing the histogram */
	for(c1 = 0; c1 < index_number; c1++)
	{
		pair_correlation_function[c1] /= (normalizing_constant*c1);
	}
	

	
	
	/* printing the hystogram in the file */
	for(c1 = 0; c1 < index_number; c1++)
	{	
		fprintf(ofp_pair_correlation,"%lf   %lf   %d\n",(c1 + 0.5)/simdatap->distanceResolution ,pair_correlation_function[c1], simdatap->currConfig);
	}
	fprintf(ofp_pair_correlation,"\n");


	/* closing the filepointer */
	fclose(ofp_pair_correlation);

	/* freeing allocated memory */
	free(pair_correlation_function);
	
	return 0;
}
