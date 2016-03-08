/**************************************************
 * Functions to initialize the capillary distribution 
 * and the resulting interaction energy of the particles
 * calculated for each capillary
 * the idea for the later one
 */





/*************************************************
 * function that creates a random capillary distribution
 * in the hexagonal simulation area 
 * here it is used as the initial distribution of capillaries !!!!!
 */
 
 
/* Standard includes */
# include<stdio.h>
# include<stdlib.h>
# include<math.h>
# include<time.h>

/* eigene includes */

# include "simstruct.h"

/* #define M_PI 3.14159265
 */

int initialize_capillary_distribution(simstruct *simdatap)
{
	int c1;
	double sum_of_circle_areas;
	
	
	/* preliminary check if the desired capillary distribution is 
	 * prinicpally to dense or not. The highest Volume fraction possible
	 * is Pi/2/sqrt(3) which is about 90.7 percent
	 */
	 sum_of_circle_areas = M_PI*simdatap->rCapillary*simdatap->rCapillary*simdatap->capillaryNumber;
	 
	 if(sum_of_circle_areas > simdatap->max_vol_fraction*simdatap->simArea)
	 {
	 	printf("With a volumefraction of %lf and above you are close to optimal packing\n"
				"The current volume fraction is %lf \n"
				"This is very unlikely to achive by random numbers, try it with a regular lattice!\n", simdatap->max_vol_fraction, sum_of_circle_areas/simdatap->simArea);
	 	exit(0);
	 }
	 
	 
	
	for(c1 = 0; c1 < simdatap->capillaryNumber; c1++)
	{
		simdatap->currCapillary = c1;
		
		do
		{
			/* initializing the random position in a box enclosing the hexagon 
			 */
			simdatap->x[c1] = simdatap->square_basis_length*(ran2_NRIC(&(simdatap->randseed)) - 0.5);
			simdatap->y[c1] = simdatap->square_basis_length*(ran2_NRIC(&(simdatap->randseed)) - 0.5);
		
			ensure_in_square(simdatap);
		
		}while(initialize_capillary_overlap_test(simdatap) != 1);
	}

	return 0;
}






