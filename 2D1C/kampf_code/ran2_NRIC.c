/******************************************************************
 * this isreimplementation of ran2 from the numerical recipes in C
 * it returns double values 
 * furthermore this file contains the method for generating gaussian
 * distributed random numbers provided by gasdev in NRiC
 */
 
 
/* Standard includes */
# include<stdio.h>
# include<stdlib.h>
# include<math.h>

/* eigene includes */

# include "simstruct.h"

/* predefines by NRIC */

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)



double ran2_NRIC(long *idum)
{
	int j;
  	long k;
  	static long idum2=123456789;
  	static long iy=0;
  	static long iv[NTAB];
  	float temp;
  	if (*idum <= 0) 
	{                      /* Initialize */
    	if (-(*idum) < 1)
		{
		 	*idum=1;         /* Be sure to prevent idum = 0 */
		}
		else
		{
		 	*idum = -(*idum);
		}
		
      	idum2=(*idum);
      	
		for (j=NTAB+7;j>=0;j--) 
		{          /* Load the shuﬄe table (after 8 warm-ups). */
          	k=(*idum)/IQ1;
          	*idum=IA1*(*idum-k*IQ1)-k*IR1;
        
		  	if (*idum < 0)
			{
				*idum += IM1;
			}
			else{}
          	if (j < NTAB)
			{
				iv[j] = *idum;
			}
			else{}
      	}
      	iy=iv[0];
  	}
  	
	k=(*idum)/IQ1;                         /* Start here when not initializing. */
  	
	*idum=IA1*(*idum-k*IQ1)-k*IR1;         /* Compute idum=(IA1*idum) % IM1 without */
  	
	if (*idum < 0) 
	{
		*idum += IM1;           /* overﬂows by Schrage’s method. */
	}
	else{};
	
  	k=idum2/IQ2;
  	idum2=IA2*(idum2-k*IQ2)-k*IR2;         /* Compute idum2=(IA2*idum) % IM2 likewise. */
  	
	if (idum2 < 0)
	{
		idum2 += IM2;
	}
	else{};
	
  	j=iy/NDIV;                             /* Will be in the range 0..NTAB-1. */
  	iy=iv[j]-idum2;                        /* Here idum is shuﬄed, idum and idum2 are */
  	iv[j] = *idum;                         /* combined to generate output. */
  	

	if (iy < 1) 
	{
		iy += IMM1;
	}
	else{};
  	if ((temp=AM*iy) > RNMX) 
	{
		return RNMX;  /* Because users do not expect endpoint values. */
	}
  	else
	{
		return temp;
	}
	
}




double gasdev_NRIC(long *idum)
{
    double ran2_NRIC(long *idum);
    static int iset=0;
    static double gset;
    double fac,rsq,v1,v2;
	
    if (*idum < 0) 
	{
		iset=0;                       /* Reinitialize.*/
	}
	else{};
    if (iset == 0) 
	{                             /* We don’t have an extra deviate handy, so */
        do 
		{
             v1=2.0*ran2_NRIC(idum)-1.0;              /* pick two uniform numbers in the square ex- */
             v2=2.0*ran2_NRIC(idum)-1.0;              /* tending from -1 to +1 in each direction,*/
             rsq=v1*v1+v2*v2;                    /* see if they are in the unit circle, */

      	} while (rsq >= 1.0 || rsq == 0.0);      /* and if they are not, try again. */
      
	  	fac=sqrt(-2.0*log(rsq)/rsq);
 
      	/* Now make the Box-Muller transformation to get two normal deviates. Return one and
       	 * save the other for next time.
	   	 */
      	gset=v1*fac;
      	iset=1;                               /* Set ﬂag. */
      	return v2*fac;
  	} 
	else 
	{                                  /* We have an extra deviate handy, */
      	iset=0;                               /* so unset the ﬂag, */
      	return gset;                          /* and return it. */
  	}
}
