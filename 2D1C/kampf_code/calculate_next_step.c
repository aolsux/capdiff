/********************************************************************
 * Function to calculate the next intermediate step to create a new
 * capillary distribution according to the model described in
 * Physica A 369 599-611 (2006)
 */


/* Standard includes */
# include<stdio.h>
# include<stdlib.h>
# include<math.h>
# include<time.h>

/* eigene includes */

# include "simstruct.h"


int calculate_next_intermed_step(simstruct *simdatap)
{
	int c1;
	double acceptance_threshold; 
	double randwert1,randwert2;
	double cur_x;
	double cur_y;
	

	/* saving the current capillary distribution */
	for(c1 = 0; c1 < simdatap->capillaryNumber; c1++)
	{
		simdatap->old_x[c1] = simdatap->x[c1];
		simdatap->old_y[c1] = simdatap->y[c1];
	}

	/* saving the current energy state */
	simdatap->old_energy = simdatap->cur_energy;
	
		
	/* loop is done until a new valid configuration is reached
	 */
	do
	{
		/* restoring the last accepted capillary distribution 
		 */
		for(c1 = 0; c1 < simdatap->capillaryNumber; c1++)
		{
			simdatap->x[c1] = simdatap->old_x[c1];
			simdatap->y[c1] = simdatap->old_y[c1];
		}
				

		/* This is the whole difference!!!!!!
		 * instead of changing all capillaries only a fraction is chosen 
		 * by a random number and only this is varied
	 	*/
		
	
		for(c1 = 0; c1 < simdatap->fraction_to_change; c1++)
		{
			/* do the loop until a valid new capillary distribution is found 
			 */
			do
			{
				/* restoring the old capillary distribution */
				simdatap->x[simdatap->currCapillary] = simdatap->old_x[simdatap->currCapillary];
				simdatap->y[simdatap->currCapillary] = simdatap->old_y[simdatap->currCapillary];
			
				/* choose the capillary by a random number */
				randwert1 = ran2_NRIC(&(simdatap->randseed));
				simdatap->currCapillary = (int)floor(simdatap->capillaryNumber*randwert1);
				cur_x = simdatap->old_x[simdatap->currCapillary];
				cur_y = simdatap->old_y[simdatap->currCapillary];

				randwert1 = ran2_NRIC(&(simdatap->randseed));
				simdatap->x[simdatap->currCapillary] = cur_x + (2.0*simdatap->max_x_deviation*randwert1 - simdatap->max_x_deviation);

				randwert2 = ran2_NRIC(&(simdatap->randseed));
				simdatap->y[simdatap->currCapillary] = cur_y + (2.0*simdatap->max_y_deviation*randwert2 - simdatap->max_y_deviation);
					
				/* ensures that the capillary on the new position is inside the rectancle */
				
				ensure_in_square(simdatap);

			}while(test_overlap_of_current_capillary(simdatap) != 1); /* one as result of the function means that ther is no overlap */

		}
	
	
		/* calculating the energy of the configuration
		 */
		calculate_variable_energy_part(simdatap);
		
		 
		/* calculate the acceptance threshhold 
		 * this includes the renormalization of the complete energy change to the 
		 * energychange per particle which is independent from the changed particle number (on average)
		 * and thus the quantity that should be relevant. In any other case the threshhold is dependent
		 * on the number of changed particles or it would result in different gamma values for a different
		 * number of particles!s
		 */
		acceptance_threshold = fmin(1.0,exp(-(simdatap->cur_energy - simdatap->old_energy)*simdatap->coupling_parameter/((double)(simdatap->fraction_to_change))));
	}while(acceptance_threshold < ran2_NRIC(&(simdatap->randseed)));  /* if the caluclated wall is smaller than the random value the thermal disorder is to large to be acceptable */ 

	return 0;
}




