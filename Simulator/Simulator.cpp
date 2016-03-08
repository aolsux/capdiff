/*
 * Simulator.cpp
 *
 *  Created on: Aug 13, 2009
 *      Author: gstu0908
 */

#ifndef SIMULATOR_CPP_
#define SIMULATOR_CPP_

#include "Simulator.h"
#include "util/Field2D.h"
#include "util/LiveField2D.h"

#include <MyException.h>
#include <RandomNumberGenerator.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <float.h>

#include <Progress.h>

template<typename B,typename F>
Simulator<B,F>::Simulator(B b, F f, double _time, double _diffusion, unsigned _steps, double _omega0):
	//miscellaneous
	compartement		(NULL),
	boundary			(b),
	//time settings
	steps				(_steps),
	timeIntervall		(_time),
	timeStepSize		(timeIntervall / steps),
	//magnetic field
	Bz					(0),
	theta				(0),
	deltaChi			(0),
	gyromagneticRatio	(0),
	deltaomega0			(_omega0),
	field				(f),
	//diffusion settings
	diffusionConstant	(_diffusion),
	sigma				(sqrt(2.0 * diffusionConstant * timeStepSize)),
	rng					(1)
{
}

template<typename B,typename F>
Simulator<B,F>::~Simulator()
{
	delete compartement;
}

template<typename B,typename F>
void Simulator<B,F>::optimizeTimeSettings()
{
	double minradius = capillaryConfiguration.getMinimalRadius();

	//only increase step size
	if (sigma < minradius / STEPLENGTH_RADIUS_RATIO)return;
	//optimum already, nothing needs to change
	if (steps == 2.0 * diffusionConstant * STEPLENGTH_RADIUS_RATIO * STEPLENGTH_RADIUS_RATIO * timeIntervall / (minradius * minradius))return;

	cout << "Number of steps will get optimized for given parameters:"<<endl;
	cout << "\tr_min= " << minradius*1E6 << endl;
	cout << "\tD    = " << diffusionConstant*1E9 << endl;
	cout << "\tT    = " << timeIntervall <<endl;
	cout << "\ts    = " << steps << endl<< "\t->\t";
	steps = 2.0 * diffusionConstant * STEPLENGTH_RADIUS_RATIO * STEPLENGTH_RADIUS_RATIO * timeIntervall / (minradius * minradius);
	cout << steps << endl;
	cout << "\tdt   = " << timeStepSize*1E3 << "\t->\t";
	timeStepSize = timeIntervall / steps;
	cout << timeStepSize*1E3 << endl;
	cout << "\tsigma= " << sigma*1E6 << "\t->\t";
	sigma = sqrt(2.0 * diffusionConstant * timeStepSize);
	cout << sigma*1E6 << endl;
}

template<typename B,typename F>
void Simulator<B,F>::setCapillaryConfiguration(CapillaryConfiguration config)
{
	collisionObjects.clear();
	capillaryConfiguration = config;
	boundary.handleCapillaryCollisions(capillaryConfiguration);
	for (unsigned int i = 0; i < capillaryConfiguration.capillaries.size(); i++)
			collisionObjects.push_back(&(capillaryConfiguration.capillaries[i]));

	if (compartement != NULL) delete compartement;
	compartement = new Compartement(capillaryConfiguration.getSystemSize(), collisionObjects);
	//check if mean step size fits capillary radius
	if ((sigma > capillaryConfiguration.getMinimalRadius() / STEPLENGTH_RADIUS_RATIO) && (capillaryConfiguration.getMinimalRadius()!=0)){
		cerr << "Mean step size is too large for smallest capillary:" << endl;
		cerr << "\tr_min= " << capillaryConfiguration.getMinimalRadius()*1E6 << endl;
		cerr << "\tD    = " << diffusionConstant*1E9<< endl;
		cerr << "\tT    = " << timeIntervall <<endl;
		cerr << "\tdt   = " << timeStepSize <<endl;
		cerr << "\ts    = " << steps << endl;
		cerr << "\tsigma= " << sigma*1E6 << endl;
		optimizeTimeSettings();
	}
}

template<typename B,typename F>
double *Simulator<B,F>::runTrajectory() const
{
	double *currentTrajectory = new double[steps];
	Point currentLocation = getRandomStartLocation();
	currentTrajectory[0] = 0;
	for (unsigned i = 1; i < steps; i++) {//currentStep starts with 1
		//phase aufkumulieren
		currentTrajectory[i] = field.getValue(currentLocation) * timeStepSize * deltaomega0;
		//nachste position berechnen
		currentLocation = getNextStepLocation(currentLocation);
	}
	return currentTrajectory;
}

template<typename B,typename F>
Matrix2D Simulator<B,F>::runSpacialResolvedTrajectory() const
{
	Matrix2D matrix(steps,3);
	Point currentLocation = getRandomStartLocation();

	matrix(0,0)=0;
	matrix(0,1)=currentLocation.x;
	matrix(0,2)=currentLocation.y;

	for (unsigned i = 1; i < steps; i++) {//currentStep starts with 1
		currentLocation = getNextStepLocation(currentLocation);
		matrix(i,0) = field.getValue(currentLocation) * timeStepSize * deltaomega0;
		matrix(i,1) = currentLocation.x;
		matrix(i,2) = currentLocation.y;
	}

	return matrix;
}

template<typename B,typename F>
Point Simulator<B,F>::getRandomStartLocation() const
{
	Point p = boundary.getRandomInsideLocation(rng);
	while (compartement->getCollision(p) != NULL) {
		p = boundary.getRandomInsideLocation(rng);
	}
	return p;
}

template<typename B,typename F>
Point Simulator<B,F>::getNextStepLocation(const Point &currentLocation) const
{
	Point newLocation = currentLocation + rng.getNormalDistributedPoint(sigma);
	/*Check for Collision with boundary and Collision objects*/
	CollisionObject const *obj = (!boundary.isInside(newLocation)? &boundary : compartement->getCollision(newLocation));
	while(obj!=NULL){
		newLocation = obj->collide(newLocation);
		obj = (!boundary.isInside(newLocation)? &boundary : compartement->getCollision(newLocation));
	}
	return newLocation;
}


template<typename B,typename F>
Matrix2D Simulator<B,F>::getErrorMap(unsigned samples, double _sigma)
{
	Matrix2D errorMap(samples, 8);
	sigma = _sigma;
	double dt = sigma * sigma / (2 * diffusionConstant);
	LiveField2D liveField(capillaryConfiguration);
	//ProgressBar progress("Sampling " + util::num2str(samples) + " time steps", samples);
	for (unsigned n = 0; n < samples; n++) {
		Point start = getRandomStartLocation();
		Point p = start;
		double grid_field_phase_t = 0;
		double live_field_phase_t = 0;
		for (double t = 0; t < timeStepSize; t += dt) {
			p = this->getNextStepLocation(p);
			grid_field_phase_t += field.getValue(p) * dt * deltaomega0;
			live_field_phase_t += liveField.getValue(p) * dt * deltaomega0;
		}
		if ((start - p).abs() >= 10 * sqrt(2.0 * diffusionConstant * timeStepSize)) continue;
		Point mitte = (start + p) * 0.5;
		double grid_field_phase_T = field.getValue(mitte) * timeStepSize * deltaomega0;
		double live_field_phase_T = liveField.getValue(mitte) * timeStepSize * deltaomega0;
		errorMap(n, 0) = start.x;
		errorMap(n, 1) = start.y;
		errorMap(n, 2) = p.x;
		errorMap(n, 3) = p.y;
		errorMap(n, 4) = grid_field_phase_t;
		errorMap(n, 5) = grid_field_phase_T;
		errorMap(n, 6) = live_field_phase_t;
		errorMap(n, 7) = live_field_phase_T;
		//progress.progress();
	}
	sigma = sqrt(2.0 * diffusionConstant * timeStepSize);
	return errorMap;
}

template<typename B,typename F>
	double Simulator<B,F>::getVolumeFraction()									const{ return (capillaryConfiguration.getCapillaryArea() / boundary.getArea());}
template<typename B,typename F>
	double Simulator<B,F>::getEnvironmentArea()									const{ return boundary.getArea();}
template<typename B,typename F>
	Rectangle Simulator<B,F>::getRectangle()									const{ return boundary.getRectangle();}
template<typename B,typename F>
	double Simulator<B,F>::getTimeStepSize()									const{ return timeStepSize;}
template<typename B,typename F>
	double Simulator<B,F>::getTimeIntervall()									const{ return timeIntervall;}
template<typename B,typename F>
	double Simulator<B,F>::getDiffusionConsant()								const{ return diffusionConstant;}
template<typename B,typename F>
	const F& Simulator<B,F>::getField()											const{ return field;}
template<typename B,typename F>
	const B& Simulator<B,F>::getBoundaryCondition()								const{ return boundary;}
template<typename B,typename F>
	unsigned Simulator<B,F>::getTimeSteps()										const{ return steps;}
template<typename B,typename F>
	double Simulator<B,F>::getMeanStepSize()									const{ return sigma;}
template<typename B,typename F>
	double Simulator<B,F>::getOmega0()											const{ return deltaomega0;}
template<typename B,typename F>
	const CapillaryConfiguration& Simulator<B,F>::getCapillaryConfiguration()	const{ return capillaryConfiguration;}
template<typename B,typename F>
	const Compartement&	Simulator<B,F>::getCollisionTree()						const{ return *compartement;}

#endif
