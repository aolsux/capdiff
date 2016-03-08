/**********************************************************
 * function that initializes the simulation structure with
 * resonable default values
 */
 


/* Standard includes */
# include<stdio.h>
# include<stdlib.h>
# include<math.h>
# include<time.h>

/* eigene includes */

# include "simstruct.h"

 
int initialize_structure(simstruct *simdatap)
{
	int c1;

	/* initializing the random seed by the system time
	 * according to NRIC the seed should be negative
	 * thus the negative of the system time is used
	 */
	 
	simdatap->randseed = -((long)time(NULL));



	
	/* initializing the independent scalar variables
	 */
	
	simdatap->capillaryNumber = DEF_VAL_CAPILLARYNUMBER;	
	simdatap->simArea = DEF_VAL_SIMAREA;
	simdatap->rCapillary = DEF_VAL_RCAPILLARY;
	simdatap->max_vol_fraction = DEF_VAL_MAX_VOL_FRAC;
	simdatap->number_of_configurations = DEF_VAL_NUM_OF_CONFIG;
	simdatap->number_of_intermed_steps = DEF_VAL_NUM_OF_INT_STEPS;
	simdatap->max_x_deviation = DEF_VAL_MAX_Y_DEVIATION;
	simdatap->max_y_deviation = DEF_VAL_MAX_Y_DEVIATION;
	simdatap->damping_const_square = DEF_VAL_DAMP_CONST_SQ;
	simdatap->min_neighoring_x_number = DEF_VAL_MIN_NEIGH_X_NUM;
	simdatap->max_neighoring_x_number = DEF_VAL_MAX_NEIGH_X_NUM;
	simdatap->min_neighoring_y_number = DEF_VAL_MIN_NEIGH_Y_NUM;
	simdatap->nax_neighoring_y_number = DEF_VAL_MAX_NEIGH_Y_NUM;
	simdatap->coupling_parameter = DEF_VAL_COUPLING_PRARAM;
	simdatap->distanceResolution = DEF_VAL_DIST_RES_IN_HYST;
	simdatap->min_neighoring_kx_number = DEF_VAL_MIN_NEIGH_KX_NUM;
	simdatap->max_neighoring_kx_number = DEF_VAL_MAX_NEIGH_KX_NUM;
	simdatap->min_neighoring_ky_number = DEF_VAL_MIN_NEIGH_KY_NUM;
	simdatap->nax_neighoring_ky_number = DEF_VAL_MAX_NEIGH_KY_NUM;
	simdatap->fraction_to_change = (int)fmax(1.0,DEF_VAL_CAPILLARYNUMBER*DEF_VAL_FRAC_TO_CHANGE);	/* ensures that at least one particle is moved!*/
	
	




	
	
	/* initializing the dependent scalar variables
	 */
	 
	simdatap->square_basis_length = sqrt(simdatap->simArea);
	
	
	/* initializing the lattice vectors */
	simdatap->lattice_vector_ax = simdatap->square_basis_length;
	simdatap->lattice_vector_ay = 0.0;
	
	simdatap->lattice_vector_bx = 0.0;
	simdatap->lattice_vector_by = simdatap->square_basis_length;
	
	simdatap->lattice_inv_det = simdatap->lattice_vector_ax*simdatap->lattice_vector_by - simdatap->lattice_vector_ay*simdatap->lattice_vector_bx;
	
	
	
	/* dependent fields and variables
	 */
	 
	for(c1 = 0; c1 < DEF_VAL_CAPILLARYNUMBER; c1++)
	{
		simdatap->x[c1] = 0.0;
		simdatap->y[c1] = 0.0;

		simdatap->old_x[c1] = 0.0;
		simdatap->old_y[c1] = 0.0;
	}




	return 0;
}
