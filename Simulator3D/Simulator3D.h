/*
 * Simulator.h
 *
 *  Created on: Aug 13, 2009
 *      Author: gstu0908
 */

#ifndef SIMULATOR3D_H_
#define SIMULATOR3D_H_

#include <vector>
#include <RandomNumberGenerator.h>
#include <Matrix2D.h>
#include <Rectangle.h>

#include "CollisionObject3D.h"
#include "Configuration3D.h"
#include "Compartement3D.h"

using namespace std;

template<typename Boundary_Type, typename Field_Type>
class Simulator3D
{
	public:
		Simulator3D(Boundary_Type boundary, Field_Type field, double time, double diffusion, unsigned steps, double omega0);

		/*
		 * deletes allocated memory...
		 */
		virtual ~Simulator3D();

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
		void setConfiguration(SphereConfiguration caps);

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
			  Cuboid					getSimulationBox()			const;
			  double 					getVolumeFraction() 		const;
		const SphereConfiguration&		getCapillaryConfiguration() const;
		const Field_Type&				getField()					const;
		const Boundary_Type&			getBoundaryCondition()		const;
		const Compartement3D&			getCollisionTree()			const;

		inline Point3D getNextStepLocation(const Point3D &currentLocation) const;		/* Returns a valid next step location.*/
		inline Point3D getRandomStartLocation() 						 const;		/* Returns a valid location.*/

	private:
		/*Hide default constructor.*/
		Simulator3D();
		/*Make Simulator non copy able by hiding copy constructor and assignment operator.*/
		Simulator3D(const Simulator3D &ref);
		Simulator3D &operator=(const Simulator3D &ref);

		//Spatial settings
		SphereConfiguration capillaryConfiguration;
		Compartement3D *compartement;
		vector<CollisionObject3D const *> collisionObjects;
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

#include "Simulator3D.cpp"

#endif /*SIMULATOR3D_H_*/
