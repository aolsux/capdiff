/*
 * Simulator.h
 *
 *  Created on: Aug 13, 2009
 *      Author: gstu0908
 */

#ifndef SIMULATOR_H_
#define SIMULATOR_H_

#include <vector>
#include <RandomNumberGenerator.h>
#include <Matrix2D.h>
#include <Rectangle.h>

#include "util/CollisionObject.h"
#include "util/CapillaryConfiguration.h"
#include "util/Describable.h"
#include "util/Compartement.h"

using namespace std;

template<typename Boundary_Type, typename Field_Type>
class Simulator: public Describable
{
	public:
		Simulator(Boundary_Type boundary, Field_Type field, double time, double diffusion, unsigned steps, double omega0);

		/*
		 * deletes allocated memory...
		 */
		virtual ~Simulator();

		/*
		 * does the random walk and returns the experienced dphi(phaseshifts) for each timestep
		 */
		double* runTrajectory() const;

		/*
		 * does the random walk and returns the experienced dphi(phase shifts) for each time step and the corresponding position
		 */
		Matrix2D runSpacialResolvedTrajectory() const;

		/*
		 * Sets the environment to the given capillary configuration and manages their boundary conditions
		 */
		void setCapillaryConfiguration(CapillaryConfiguration caps);

		/*
		 * Optimizes the given time settings:
		 * The mean step size is pushed up to the limit given by the smallest capillary radius.
		 */
		void optimizeTimeSettings();

		/*
		 * Calculates a set of n random walks from location A to B with average step length sigma.
		 * These Trajectories are then treated as a single step in the original timeStepSize.
		 * The return value contains A.x,A.y,B.x,B.y,error.
		 */
		Matrix2D getErrorMap(unsigned n, double sigma);

		static const double STEPLENGTH_RADIUS_RATIO = 4;

		/*Various getters.*/
			  unsigned 					getTimeSteps() 				const;
			  double					getTimeIntervall()			const;
			  double 					getTimeStepSize() 			const;
			  double 					getMeanStepSize() 			const;
			  double 					getOmega0() 				const;
			  double 					getDiffusionConsant() 		const;
			  double					getEnvironmentArea()		const;
			  Rectangle					getRectangle()				const;
			  double 					getVolumeFraction() 		const;
		const CapillaryConfiguration&	getCapillaryConfiguration() const;
		const Field_Type&				getField()					const;
		const Boundary_Type&			getBoundaryCondition()		const;
		const Compartement&				getCollisionTree()			const;

		inline Point getNextStepLocation(const Point &currentLocation) const;
	private:
		/*Hide default constructor.*/
		Simulator();
		/*Make Simulator non copy able by hiding copy constructor and assignment operator.*/
		Simulator(const Simulator &ref);
		Simulator &operator=(const Simulator &ref);

		/* Returns a valid next step and starting location.*/
		inline Point getRandomStartLocation() const;

		//Spatial settings
		CapillaryConfiguration capillaryConfiguration;
		Compartement *compartement;
		vector<CollisionObject const *> collisionObjects;
		Boundary_Type 	boundary;

		//Time Settings
		unsigned steps;				//[1]
		double timeIntervall;		//[s], Total simulation Duration
		double timeStepSize; 		//[s]

		//Magnetic Field
		double Bz;					//[T]
		double theta;				//[radiant] Tilting angle between magnetic field and capillary
		double deltaChi;			//[1]
		double gyromagneticRatio;	//[1/(sT)]
		double deltaomega0;			//[radiant/s]
		Field_Type 		field;

		//Diffusion
		double diffusionConstant;	//[m*m/s]
		double sigma;				//[m], mean step length of random walk (s = 2*D*dT)
		RandomNumberGenerator rng;
};

#include "Simulator.cpp"

#endif
