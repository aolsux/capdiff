/******************************************************
 * This file contains the header for the structure which contains
 * all necessary information on the Simulation data for the creation
 * of a two dimensional capillary distribution by a 2d- one component plasma
 * This is after phyica A Vol 369 (2006) pp 599-611
 *
 * it provides the independt variables
 *		capillaryNumber 			= Number of capillaries in the specified area which is
 *										periodically continued
 *		simArea						= Size of the Hexagonal Area which is continued. The Area is given 
 * 										in the units of rCapillary squared
 *		rCapillary					= radius of a single capillary (necessary to prevent to close spacing
 *										of the capillaries
 * 		max_vol_fraction			= maximal allowed volume fraction for savety that it is not a to dense 
 *										packing of the capillaries (should not be necessary is just for savety)
 *		number_of_configuration		= number of the configurations to be calculated for the given coupling parameter
 *
 *		number_of_intermed_steps	= number of intermediate step to be performed to get a new capillary distribution 
 *										that is not to closely related to the capillary distribution before
 *		coupling_parameter			= is the parameter that determines the strength of disorder in the capillary distribution
 *
 *		max_x_deviation				= is the length of the maximum deviation of a capillary form its current position 
 *										in x-direction for one reordering step
 *		max_y_deviation				= is the length of the maximum deviation of a capillary form its current position 
 *										in y-direction for one reordering step
 *
 *
 * it provides dependent areas
 *		x = vector with the x-value of the position of the capillaries
 * 		y = vector with the y-value of the position of the capillaries	
 */



/* predefined values of the independent parameters
 */
#define DEF_VAL_CAPILLARYNUMBER 	1000
#define DEF_VAL_SIMAREA				400.0
#define DEF_VAL_RCAPILLARY			0.2	
#define DEF_VAL_MAX_VOL_FRAC		0.8				/* theoretical maximum is about .9069.... */
#define DEF_VAL_NUM_OF_INT_STEPS	20000
#define DEF_VAL_NUM_OF_CONFIG		20
#define DEF_VAL_COUPLING_PRARAM		1.5	
#define DEF_VAL_MAX_X_DEVIATION		0.1*sqrt(DEF_VAL_SIMAREA/M_PI/DEF_VAL_CAPILLARYNUMBER)		/* This value is not neccessary it is just used this way and can be reset to other arbitrary values */
#define DEF_VAL_MAX_Y_DEVIATION		DEF_VAL_MAX_X_DEVIATION		/* This value is not neccessary it is just used this way and can be reset to other arbitrary values */
#define DEF_VAL_DAMP_CONST_SQ		36.0/DEF_VAL_SIMAREA	/* taken from the paper assumed to have good convergence */	
#define DEF_VAL_MIN_NEIGH_X_NUM		-1
#define DEF_VAL_MAX_NEIGH_X_NUM		1
#define DEF_VAL_MIN_NEIGH_Y_NUM		-1
#define DEF_VAL_MAX_NEIGH_Y_NUM		1
#define DEF_VAL_MIN_NEIGH_KX_NUM	-15
#define DEF_VAL_MAX_NEIGH_KX_NUM	15
#define DEF_VAL_MIN_NEIGH_KY_NUM	-15
#define DEF_VAL_MAX_NEIGH_KY_NUM	15
#define DEF_VAL_DIST_RES_IN_HYST	4
#define DEF_VAL_FRAC_TO_CHANGE		1/DEF_VAL_CAPILLARYNUMBER	/* fractions of capillaries to change bevor calculating the next energy step */


/* definition of the basic structure which contains all necessary parameters for the
 * running of the simulation and has to be shared by the functions of the programm
 */
typedef struct
{
	/* independent variables 
	 */
	int capillaryNumber;			/* initialized with DEF_VAL_CAPILLARYNUMBER */
	double simArea;					/* initialized with DEF_VAL_SIMAREA */
	double rCapillary;				/* initialized with DEF_VAL_RCAPILLARY */
	double max_vol_fraction;		/* initialized with DEF_VAL_MAX_VOL_FRAC */
	int number_of_configurations;	/* initialized with DEF_VAL_NUM_OF_CONFIGS */
	int number_of_intermed_steps;	/* initialized with DEF_VAL_NUM_OF_INT_STEPS */
	double coupling_parameter;		/* initialized with DEF_VAL_COUPLING_PRARAM */
	double max_x_deviation;			/* initialized with DEF_VAL_MAX_X_DEVIATION */
	double max_y_deviation;			/* initialized with	DEF_VAL_MAX_Y_DEVIATION */
	double damping_const_square;	/* initialized with DEF_VAL_DAMP_CONST_SQ */
	int min_neighoring_x_number;	/* initialized with DEF_VAL_MIN_NEIGH_X_NUM */
	int max_neighoring_x_number;	/* initialized with DEF_VAL_MAX_NEIGH_X_NUM */
	int min_neighoring_y_number;	/* initialized with DEF_VAL_MIN_NEIGH_Y_NUM */
	int nax_neighoring_y_number;	/* initialized with DEF_VAL_MAX_NEIGH_Y_NUM */
	int min_neighoring_kx_number;	/* initialized with DEF_VAL_MIN_NEIGH_KX_NUM */
	int max_neighoring_kx_number;	/* initialized with DEF_VAL_MAX_NEIGH_KX_NUM */
	int min_neighoring_ky_number;	/* initialized with DEF_VAL_MIN_NEIGH_KY_NUM */
	int nax_neighoring_ky_number;	/* initialized with DEF_VAL_MAX_NEIGH_KY_NUM */
	int distanceResolution;			/* initialized with DEF_VAL_DIST_RES_IN_HYST */
	int fraction_to_change;			/* initialized with DEF_VAL_FRAC_TO_CHANGE */
	
	/* dependent scalars
	 */
	double square_basis_length;
	double lattice_vector_ax;
	double lattice_vector_ay;
	double lattice_vector_bx;
	double lattice_vector_by;
	double lattice_inv_det;

	
	/* dependent fields
	 */
	 
	double x[DEF_VAL_CAPILLARYNUMBER];
	double y[DEF_VAL_CAPILLARYNUMBER];

	double old_x[DEF_VAL_CAPILLARYNUMBER];
	double old_y[DEF_VAL_CAPILLARYNUMBER];
	

	/* auxilliary parameters neccesary for parameter exchange during 
	 * the running of the programm
	 */
	long randseed;						/* seed of the random number generator used */
	int currCapillary;					/* gives back the current capillary */
	int currConfig;						/* gives back the current confignumber */
	int second_capillary_for_distance_calculation;	/* is the index of the capillary the distance to the current capillary is calculated */
	double cur_energy;					/* current energy of the capillary distribution */
	double old_energy;					/* old energy of the capillary distribution */
	double other_x;						/* x-postion of the other capillary the distance is calculated at */
	double other_y;						/* y-postion of the other capillary the distance is calculated at */

} simstruct;






/* Prototypen der verwendeten Funktionen 
 */
 
int initialize_capillary_distribution(simstruct *simdatap);
int initialize_structure(simstruct *simdatap);	
int ensure_in_square(simstruct *simdatap);	
int test_capillary_overlap_all(simstruct *simdatap);	
int initialize_capillary_overlap_test(simstruct *simdatap);
int test_overlap_of_current_capillary(simstruct *simdatap);
int calculate_next_intermed_step(simstruct *simdatap);
int calculate_variable_energy_part(simstruct *simdatap);
int nearest_neighbour_dist(simstruct *simdatap);
int pair_correlation_function(simstruct *simdatap);
int capillary_cur_old_distribution_comparison(simstruct *simdatap);

double calculate_distance_squared(simstruct *simdatap);

double ran2_NRIC(long *idum);
double gasdev_NRIC(long *idum);
double expint1_NRIC(double x);



/* Prototyps of the test functions
 */
 

