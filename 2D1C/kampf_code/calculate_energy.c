/*************************************************************** 
 * Function to calculate the energy of the system
 * However this is not the energy as in Equation 11 in Physica A 
 * 369 (2006) 599-611 but the normalized energy which is the same
 * only devided by q^2 which is present in all terms
 * however since only the energy difference between two states has
 * to be calcualted the term that does not depent explicitly on r_ij
 * is neglected
 */



/* Standard includes */
# include<stdio.h>
# include<stdlib.h>
# include<math.h>
# include<time.h>

/* eigene includes */

# include "simstruct.h"


/* local prototype functions that are only of interest in the 
 * context of this function */
 
double calculate_direct_space_part(simstruct *simdatap);
double calculate_reciprocal_space_part(simstruct *simdatap);





 
int calculate_variable_energy_part(simstruct *simdatap)
{
	double direct_space_energy;
	double reciprocal_space_energy;
	
	int save_cur_capillary_index;
	
	
	/* saving the current capillary index since it is changed during the
	 * calculation of the energy
	 */
	save_cur_capillary_index = simdatap->currCapillary;
	
	/* initializing the energy parts */
	direct_space_energy = 0.0;
	reciprocal_space_energy =0.0;
	
	/* calculating the direct lattice part of the energy */
	direct_space_energy = calculate_direct_space_part(simdatap);
	
	/* calculating the reciprocal lattice part of the energy */
	reciprocal_space_energy = calculate_reciprocal_space_part(simdatap);


	/* writing the energy */
	simdatap->cur_energy = reciprocal_space_energy + direct_space_energy;

	/* writing back the current index of the capillary thus it is
	 * the same as in the beginning
	 */
	simdatap->currCapillary = save_cur_capillary_index;

	return 0;
}
 



/* function to calculate the direct space part of the interaction energy in the system
 * here the factor of q^2/4 is excluded, since it is in both terms (direct and reciprocal 
 * term)
 */
 
double calculate_direct_space_part(simstruct *simdatap)
{
	double direct_space_part;
	double current_distance_sq;
	int c1, c2, c3, c4;
	double cur_x, cur_y, other_x, other_y,lattice_vector_ax, lattice_vector_by, damping_const_square;
	int min_neighoring_x_number;
	int max_neighoring_x_number;
	int min_neighoring_y_number;
	int nax_neighoring_y_number;
	
	lattice_vector_ax = simdatap->lattice_vector_ax;
	lattice_vector_by = simdatap->lattice_vector_by;
	damping_const_square = simdatap->damping_const_square;
	
	min_neighoring_x_number = simdatap->min_neighoring_x_number;
	max_neighoring_x_number = simdatap->max_neighoring_x_number;
	min_neighoring_y_number = simdatap->min_neighoring_y_number;
	nax_neighoring_y_number = simdatap->nax_neighoring_y_number;

	
	direct_space_part = 0.0;

	/* loop over the first particle
	 */
	for(c4 = 0; c4 < simdatap->capillaryNumber; c4++)
	{
		simdatap->currCapillary = c4;
	
		/* loop over the second particles first part 
		 * capillary index smaller than index of the current capillary
		 */
	 
		for(c1 = 0; c1 < simdatap->currCapillary; c1++)
		{
			
			simdatap->second_capillary_for_distance_calculation = c1;
			current_distance_sq = calculate_distance_squared(simdatap);
		
			/* loop over all necessary neighboors in the real space
			 * includes also the original simulation volume
			 */
			for(c2 = min_neighoring_x_number; c2 < max_neighoring_x_number + 1; c2++)
			{
				for(c3 = min_neighoring_y_number; c3 < nax_neighoring_y_number + 1; c3++)
				{	
					cur_x = simdatap->x[simdatap->currCapillary];
					cur_y = simdatap->y[simdatap->currCapillary];
				
					other_x = simdatap->other_x;
					other_y = simdatap->other_y;
				
					/* calculating the current distance */
					current_distance_sq = (other_x + c2*lattice_vector_ax - cur_x)*(other_x + c2*lattice_vector_ax - cur_x);
					current_distance_sq += (other_y + c3*lattice_vector_by - cur_y)*(other_y + c3*lattice_vector_by - cur_y);
					
					direct_space_part += expint1_NRIC(damping_const_square * current_distance_sq);
				}
			}
		}
	
	

		/* loop over the second particles second part 
		 * interaction with the simillar capillary in the neighboring cells
		 */
		/* loop over all necessary neighboors in the real space
		 * includes also the original simulation volume
		 */
		for(c2 = min_neighoring_x_number; c2 < max_neighoring_x_number + 1; c2++)
		{
			for(c3 = min_neighoring_y_number; c3 < nax_neighoring_y_number + 1; c3++)
			{	
				if((c2 == 0) && (c3 == 0))
				{
					/* nothing to do since this is the self interacting which is neglected */
				}
				else
				{
					cur_x = simdatap->x[simdatap->currCapillary];
					cur_y = simdatap->y[simdatap->currCapillary];
					
					other_x = simdatap->other_x;
					other_y = simdatap->other_y;
				
					/* calculating the current distance */
					current_distance_sq = (other_x + c2*lattice_vector_ax - cur_x)*(other_x + c2*lattice_vector_ax - cur_x);
					current_distance_sq += (other_y + c3*lattice_vector_by - cur_y)*(other_y + c3*lattice_vector_by - cur_y);
				
					direct_space_part += expint1_NRIC(damping_const_square * current_distance_sq);
				}
			}
		}
 


		/* loop over the second particles third part 
		 * capillary index larger than index of the current capillary
		 */
	 
		for(c1 = simdatap->currCapillary + 1 ; c1 < simdatap->capillaryNumber; c1++)
		{
		
			simdatap->second_capillary_for_distance_calculation = c1;
			current_distance_sq = calculate_distance_squared(simdatap);
			
			/* loop over all necessary neighboors in the real space
			 * includes also the original simulation volume
			 */
			for(c2 = min_neighoring_x_number; c2 < max_neighoring_x_number + 1; c2++)
			{
				for(c3 = min_neighoring_y_number; c3 < nax_neighoring_y_number + 1; c3++)
				{	
					cur_x = simdatap->x[simdatap->currCapillary];
					cur_y = simdatap->y[simdatap->currCapillary];
				
					other_x = simdatap->other_x;
					other_y = simdatap->other_y;
					
					/* calculating the current distance */
					current_distance_sq = (other_x + c2*lattice_vector_ax - cur_x)*(other_x + c2*lattice_vector_ax - cur_x);
					current_distance_sq += (other_y + c3*lattice_vector_by - cur_y)*(other_y + c3*lattice_vector_by - cur_y);
				
					direct_space_part += expint1_NRIC(damping_const_square * current_distance_sq);
				}
			}
		}
	}	

	return direct_space_part/4.0;
}



double calculate_reciprocal_space_part(simstruct *simdatap)
{
	double reciprocal_space_part;
	int c1, c2, c3;
	double cur_x, cur_y, cur_kax, cur_kby;
	double reciprocal_lattice_vector_ax, reciprocal_lattice_vector_by; 
	int min_neighoring_kx_number;
	int max_neighoring_kx_number;
	int min_neighoring_ky_number;
	int nax_neighoring_ky_number;
	double ksquare;
	double kfactor;
	double real_sum_part;
	double imag_sum_part;

	reciprocal_space_part = 0.0;
	reciprocal_lattice_vector_ax = 2.0*M_PI/simdatap->lattice_vector_ax;
	reciprocal_lattice_vector_by = 2.0*M_PI/simdatap->lattice_vector_by;

	min_neighoring_kx_number = simdatap->min_neighoring_kx_number;
	max_neighoring_kx_number = simdatap->max_neighoring_kx_number;
	min_neighoring_ky_number = simdatap->min_neighoring_ky_number;
	nax_neighoring_ky_number = simdatap->nax_neighoring_ky_number;


	
	for(c1 = min_neighoring_kx_number; c1 < max_neighoring_kx_number + 1; c1++)
	{
		for(c2 = min_neighoring_ky_number; c2 < nax_neighoring_ky_number + 1; c2++)
		{
			if((c1 == 0) && (c2 == 0))
			{
				/* nothing to do since this is a divergent self energy term which is neglected */
			}
			else
			{
				cur_kax = c1 * reciprocal_lattice_vector_ax;
				cur_kby = c2 * reciprocal_lattice_vector_by;
				
				ksquare = cur_kax*cur_kax + cur_kby*cur_kby;
				kfactor = exp(-ksquare/4.0/simdatap->damping_const_square);
				
				/* calculation of the structural factor equivalent */
				/* loop over the particles for adding the sum */
				
				real_sum_part = 0.0;
				imag_sum_part = 0.0;
				
				for(c3 = 0; c3 < simdatap->capillaryNumber; c3++)
				{
					cur_x = simdatap->x[c3];
					cur_y = simdatap->y[c3];
				
					real_sum_part += cos(cur_kax*cur_x + cur_kby*cur_y);
					imag_sum_part += sin(cur_kax*cur_x + cur_kby*cur_y);
				}
				
				/* calculating the complete additional factor for the according 
				 * reciprocal lattice vector
				 */
				 
				reciprocal_space_part += kfactor*(real_sum_part*real_sum_part + imag_sum_part*imag_sum_part);
			}
		}
	}

	/* multiplication of the additional normalization constant */
	
	reciprocal_space_part *= M_PI/simdatap->simArea;


	return reciprocal_space_part;
}



/*******************************************************
 * reimplementation of the Algorithm to calculate the
 * Exponential integral E_1(z) from the numerical recipes
 * This is no complete reimplementation since only the integrl
 * E1 is needed. It is reimplemented as double valued function
 */
 
 

#define MAXIT 1000                    	/*Maximum allowed number of iterations.*/
#define EULER 0.57721566490153286060651209008240243104215933593992           	/*Eulers constant */
#define FPMIN 1.0e-100                	/*Close to smallest representable floating-point number.*/
#define EPS 1.0e-7                   	/*Desired relative error, not smaller than the machine precision */


double expint1_NRIC(double x)
{
    int i;
    double a,b,c,d,del,fact,h,psi,ans;

    if (x <= 0.0)
	{
		printf("error in the evaluation of the exponential integral! used x-value is negative or zero!\n");
		exit(0);
	}
    else
	{
		if (x > 1.0)
		{                     /* using the Lentz's algorithm cf. Num Recipes in C */
        	b=x+1.0;
            c=1.0/FPMIN;
            d=1.0/b;
            h=d;
            for (i=1;i<=MAXIT;i++) 
			{
            	a = -i*(i);
                b += 2.0;
                d=1.0/(a*d+b);             /* Denominators cannot be zero.*/
                c=b+a/c;
                del=c*d;
                h *= del;
                if (fabs(del-1.0) < EPS)
				{
                	ans=h*exp(-x);
                    return ans;
                }
            }
            printf("continued fraction failed in expint1_NRIC! x = %lf\n", x);
			exit(0);
        }
		else
		{                           /* Evaluate series. */
        	ans = (-log(x)-EULER);          /* Set fïrst term.*/
            fact=1.0;
            
			for (i=1;i<=MAXIT;i++) 
			{
            	fact *= -x/i;
                if (i != 0)
				{
					del = -fact/(i);
				}
                else 
				{
                	psi = -EULER;          /* Compute psi(n). */
                    
					del=fact*(-log(x)+psi);
                }
                ans += del;
                
				if (fabs(del) < fabs(ans)*EPS)
				{
					return ans;
				}
				else{};
            }
            printf("series failed in expint1_NRIC! x = %lf\n", x);
			exit(0);
        }
    
	}
  	return ans;
}
 
