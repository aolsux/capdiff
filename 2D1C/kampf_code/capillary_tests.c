/************************************************
 * File contains the functions for the variouse 
 * capillary tests
 */
 
 
/* Standard includes */
# include<stdio.h>
# include<stdlib.h>
# include<math.h>
# include<time.h>

/* eigene includes */

# include "simstruct.h"


/*#define M_PI 3.14159265
 */
 
/* Function that ensures that the capillary is inside the 
 * square
 */
 
int ensure_in_square(simstruct *simdatap)
{
	double aux1, aux2;
	double gamma, lambda;
	
	/* The idea behind this function is that there is a periodizity of the centers of each square
	 * each vector can be displayed as x' = x + n*a + m*b where b and a are non colinear vectors to 
	 * from the middle of the square to the middle of two nearest neighboor squares, n and m are 
	 * integer, and x = lambda*a + gamma*b with gamma, lambda element R [0,0.5)
	 * Damit kann man ein Gleichungssystem fuer die lambda+n, gamma+m aufbauen. 
	 * Abtrennen des ganzahligen Teils liefert dann lambda sowie gamma und damit auch x
	 * for squares the lattice vectors are quite easy: |a| = a_x und |b| = b_y the others components are zero
	 * Thus the calculation of the lambda and gamma is pretty easy since the vectors are orthogonal
	 */
	
	/* determining aux1 = lambda + n and aux2 = gamma + m 
	 * by solving the linear equation system
	 */
	aux1 = (simdatap->x[simdatap->currCapillary])/simdatap->lattice_vector_ax;
	aux2 = (simdatap->y[simdatap->currCapillary])/simdatap->lattice_vector_by;
	
	/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 *!!! THE ROUND FUNCTION EXPLICITLY REQUIRES THE -std=c99  !!!!
	 *!!! FLAG TO BE SET. OTHERWISE THERE WILL BE A WARNING    !!!!
	 *!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 */
	
	/* setting the cooridinates mirrored in the mean square
	 * first part reducing the area by truncating the integer multiple of the
	 * latice vectors */
	lambda = aux1 - round(aux1);
	gamma = aux2 - round(aux2);
	
	
	simdatap->x[simdatap->currCapillary] = lambda*simdatap->lattice_vector_ax;
	simdatap->y[simdatap->currCapillary] = gamma*simdatap->lattice_vector_by;
	
	return 0;
}








/* Function that checks if the capillaries overlap
 * each other
 * returns 1 if not and 0 if they overlap
 */

int test_capillary_overlap_all(simstruct *simdatap)
{
	int c1,c2;
	double distancesquared;
	double mindistancesquared;
	
	mindistancesquared = 4.0*simdatap->rCapillary*simdatap->rCapillary;

	for(c1 = 0; c1 < simdatap->capillaryNumber - 1 ; c1++)   /* last capillary cannot be tested against still missing capillaries thus only capillaryNumber-1 */
	{
		for(c2 = c1 + 1; c2 < simdatap->capillaryNumber; c2++)
		{
			simdatap->currCapillary = c1;
			simdatap->second_capillary_for_distance_calculation = c2;
			distancesquared = calculate_distance_squared(simdatap);

			if(distancesquared < mindistancesquared)
			{
				/* The distance between capillaries is closer than the double of the radius of both capillaries
				 * thus the capillaries overlap. Return to calling function
				 */
				return 0;
			}
			else{};
		}
	}
	
	/* reached this point -> no overlap of the capillaries
	 */
	return 1;
}






/* Function to check if the current capillary overlaps with the former in the initialization
 * returns 1 if not and 0 if they overlap
 */
int initialize_capillary_overlap_test(simstruct *simdatap)
{
	int c1;
	double distancesquared;
	double mindistancesquared;
	
	mindistancesquared = 4.0*simdatap->rCapillary*simdatap->rCapillary;
	
	/* in this case only the already existing capillaries have to be tested
	 */
	for(c1 = 0; c1 < simdatap->currCapillary; c1++)
	{
		simdatap->second_capillary_for_distance_calculation = c1;
		distancesquared = calculate_distance_squared(simdatap);
			
		if(distancesquared < mindistancesquared)
		{
			/* The distance between capillaries is closer than the double of the radius of both capillaries
			 * thus the capillaries overlap. Return to calling function
			 */
			return 0;
		}
		else{};
	
	}

	/* reached this point -> no overlap of the capillaries
	 */
	return 1;
}






/* Function to check if the current capillary overlaps with any of the other capillaries
 * returns 1 if not and 0 if they overlap
  */

int test_overlap_of_current_capillary(simstruct *simdatap)
{
	int c1;
	double distancesquared;
	double mindistancesquared;
	
	mindistancesquared = 4.0*simdatap->rCapillary*simdatap->rCapillary;

	/* testing the capillaries with smaller index */
	for(c1 = 0; c1 < simdatap->currCapillary; c1++)
	{
		simdatap->second_capillary_for_distance_calculation = c1;
		distancesquared = calculate_distance_squared(simdatap);
		
/*		printf("%d  %d  %lf  %lf  %lf   %lf  %lf  %lf\n", simdatap->currCapillary, c1, simdatap->x[simdatap->currCapillary],simdatap->y[simdatap->currCapillary], simdatap->x[c1], simdatap->y[c1],distancesquared, mindistancesquared);
 */			
		if(distancesquared < mindistancesquared)
		{
			/* The distance between capillaries is closer than the double of the radius of both capillaries
			 * thus the capillaries overlap. Return to calling function
			 */
			return 0;
		}
		else{};
	}
	
	/* testing the capillaries with higher index */
	for(c1 = simdatap->currCapillary + 1; c1 < simdatap->capillaryNumber; c1++)
	{
		simdatap->second_capillary_for_distance_calculation = c1;
		distancesquared = calculate_distance_squared(simdatap);
			
/*		printf("%d  %d  %lf  %lf  %lf   %lf  %lf  %lf\n", simdatap->currCapillary, c1, simdatap->x[simdatap->currCapillary],simdatap->y[simdatap->currCapillary], simdatap->x[c1], simdatap->y[c1],distancesquared, mindistancesquared);
 */			
		if(distancesquared < mindistancesquared)
		{
			/* The distance between capillaries is closer than the double of the radius of both capillaries
			 * thus the capillaries overlap. Return to calling function
			 */
			return 0;
		}
		else{};
	}
	
	/* reached this point -> no overlap of the capillaries
	 */
	return 1;
}






/* Function that calculates the distance between the current capillary and a second capillary
 * it returns the distance as a double value
 */
 
 double calculate_distance_squared(simstruct *simdatap)
 {
 	double cur_distance_squared;
	int c1,c2;
	double cur_x;
	double cur_y;
	double other_x;
	double other_y;
	
	c1 = simdatap->second_capillary_for_distance_calculation;
	c2 = simdatap->currCapillary;
	
	cur_x = simdatap->x[c2];
	cur_y = simdatap->y[c2];
	
	other_x = simdatap->x[c1];
	other_y = simdatap->y[c1];
 
	/* the idea is to check what is the closest periodically continued 
	 * copy of the second capillary
	 * since all capillary ar inside the mean volume this is easy by checking
	 * if the difference in each direction is larger as basis_length/2.0 of the square
	 * or not and to test if the difference is negative or positiv
	 * Thus one can deside if the second capillary inside the main square has to be used
	 * or one in the neighbouring squares!
	 */
	
	
	/* find the closest x of the second capillary */
	if((cur_x - other_x) > simdatap->square_basis_length/2.0)
	{
		other_x += simdatap->square_basis_length;
	}
	else if((other_x - cur_x) > simdatap->square_basis_length/2.0)
	{
		other_x -= simdatap->square_basis_length;
	}
	else{}
	
	/* find the closest y of the second capillary */
	if((cur_y - other_y) > simdatap->square_basis_length/2.0)
	{
		other_y += simdatap->square_basis_length;
	}
	else if((other_y - cur_y) > simdatap->square_basis_length/2.0)
	{
		other_y -= simdatap->square_basis_length;
	}
	else{}
	
	cur_distance_squared = (other_x - cur_x)*(other_x - cur_x) + (other_y - cur_y)*(other_y - cur_y);
	
	simdatap->other_x = other_x;
	simdatap->other_y = other_y;
	 	
 	return cur_distance_squared;
 }



/* Function that compares the current and old distribution and checks how
 * many capillary possitions are different
 */

int capillary_cur_old_distribution_comparison(simstruct *simdatap)
{
	int c1;
	
	/* difference is only allowed in the current capillary */
	
	for(c1 = 0; c1 < simdatap->capillaryNumber; c1++)
	{
		if((simdatap->old_x[c1] == simdatap->x[c1]) && (simdatap->old_y[c1] == simdatap->y[c1]))
		{
			/* everything is ok, nothing to concern */
		}
		else if( c1 == simdatap->currCapillary)
		{
			/* still everthing is ok since that is the capillary that is allowed to be different */
		}
		else
		{
			printf("Fucking crap this is not allowed to be happen!!\n");
			printf("Only the capillary %d is allowed to change but also the capillary %d is different\n", simdatap->currCapillary, c1);
			exit(0);
		
		}
	
	
	
	
	}

	return 0;
} 
