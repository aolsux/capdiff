/*
 * SimulationEvaluator.cpp
 *
 *  Created on: Aug 18, 2009
 *      Author: gstu0908
 */

#ifndef SIMULATIONEVALUATOR3D_CPP_
#define SIMULATIONEVALUATOR3D_CPP_

#include "SimulationEvaluator3D.h"
#include <omp.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include <util.h>
#include <Progress.h>
#include <Matrix3D.h>
#include <MyException.h>

template<typename S>
SimulationEvaluator3D<S>::SimulationEvaluator3D(const S &_sim, unsigned _samples) :
	simulator		(_sim),
	real_signal		(1, simulator.getTimeSteps()),
	imag_signal		(1, simulator.getTimeSteps()),
	real_error		(1, simulator.getTimeSteps()),
	spatial_real	(1,1,1),
	spatial_imag	(1,1,1),
	spatial_error	(1,1,1),
	spatial_samples	(1,1,1),
	sequences		(1, simulator.getTimeSteps())
{
	#pragma omp parallel
		{
			if (omp_get_thread_num() == 0) nthreads = omp_get_num_threads();
		}
	for (unsigned t = 0; t < simulator.getTimeSteps(); t++)
		sequences(0, t) = 1;
	samples 		= _samples;
	sequence_count 	= 1;
}

template<typename S>
SimulationEvaluator3D<S>::~SimulationEvaluator3D()
{
	//nothing to do yet...
}

template<typename S>
void SimulationEvaluator3D<S>::evaluateTrajectories()
{
	ProgressBar progress("Sampling " + util::num2str(samples) + " Trajectories (" + util::num2str(nthreads) + "CPUs)", samples);
	omp_lock_t my_lock;
	omp_init_lock(&my_lock);
#pragma omp parallel for shared(progress, my_lock)
	for (unsigned i = 0; i < samples; i++) {
		//get random walk trajectory
		double *p = simulator.runTrajectory();
		for (unsigned s = 0; s < sequence_count; s++) {
			//reset accumulated phase
			double phi = 0;
			//calculate signal and errors
			for (unsigned t = 0; t < simulator.getTimeSteps(); t++) {
				phi += p[t] * sequences(s, t);
				double real = cos(phi);
				double imag = sin(phi);
				real_signal(s, t) += real;
				imag_signal(s, t) += imag;
				real_error(s, t) += real * real;
			}
		}
		delete[] p;
		omp_set_lock(&my_lock);
		progress.progress();
		omp_unset_lock(&my_lock);
	}
	omp_destroy_lock(&my_lock);
	for (unsigned s = 0; s < sequence_count; s++) {
		for (unsigned t = 0; t < simulator.getTimeSteps(); t++) {
			//verschiebungssatz + stichprobenvarianz
			double var = (real_error(s, t) - real_signal(s, t) * real_signal(s, t) / samples) / samples;
			//standartfehler aus varianz
			real_error(s, t)   = sqrt(var / (samples - 1));
			real_signal(s, t) /= samples;
			imag_signal(s, t) /= samples;
		}
	}
}


//template<typename S>
//void SimulationEvaluatorT<S>::getMagnetizationDistribution(double spaceres, unsigned timeres)
//{
//	double width 	= simulator.getRectangle().getWidth();
//	double height 	= simulator.getRectangle().getHeight();
//	unsigned size_x = width/spaceres+1;
//	unsigned size_y = height/spaceres+1;
//	unsigned size_t = simulator.getTimeSteps()/timeres+1;
//	Point ll=simulator.getRectangle().getLowerLeftCorner();
//
//	spatial_real	= Matrix3D(size_x,size_y,size_t);
//	spatial_imag	= Matrix3D(size_x,size_y,size_t);
//	spatial_error	= Matrix3D(size_x,size_y,size_t);
//	spatial_samples	= Matrix3D(size_x,size_y,size_t);
//
//	ProgressBar progress("Sampling Trajectories (" + util::num2str(nthreads) + "CPUs)", samples);
//	omp_lock_t my_lock;
//	omp_init_lock(&my_lock);
//	#pragma omp parallel for shared(progress, my_lock)
//		for (unsigned i = 0; i < samples; i++)
//		{
//			//get random walk trajectory
//			Matrix2D p = simulator.runSpacialResolvedTrajectory();
//
//			//reset accumulated phase
//			double phi = 0;
//			//calculate signal and errors
//			for (unsigned t = 0; t < simulator.getTimeSteps(); t++) {
//				phi += p(t,0);
//				if(t%timeres==0){
//					unsigned x = (p(t,1)-ll.x) / spaceres;
//					unsigned y = (p(t,2)-ll.y) / spaceres;
//
//					double real = cos(phi);
//					double imag = sin(phi);
//					omp_set_lock(&my_lock);
//					spatial_real   (x,y,t/timeres) += real;
//					spatial_imag   (x,y,t/timeres) += imag;
//					spatial_error  (x,y,t/timeres) += real * real;
//					spatial_samples(x,y,t/timeres) ++;
//					omp_unset_lock(&my_lock);
//				}
//			}
//			progress.progress();
//		}
//
//	progress.finish();
//	omp_destroy_lock(&my_lock);
//
//	for (unsigned x = 0; x < size_x; x++) {
//		for (unsigned y = 0; y < size_y; y++) {
//			for (unsigned t = 0; t < size_t; t++) {
//				//verschiebungssatz + stichprobenvarianz
//				double var = (spatial_error(x,y,t) - spatial_real(x,y,t) * spatial_real(x,y,t) / spatial_samples(x,y,t)) / spatial_samples(x,y,t);
//				//standartfehler aus varianz
//				spatial_error(x,y,t)  = sqrt(var / (spatial_samples(x,y,t) - 1));
//				spatial_real(x,y,t) /= spatial_samples(x,y,t);
//				spatial_imag(x,y,t) /= spatial_samples(x,y,t);
//			}
//		}
//	}
//}


//template<typename S>
//void SimulationEvaluatorT<S>::preparePeriodicMultiSpinEcho(unsigned number, unsigned periodResolution)
//{
//	sequence_count = number;
//	sequences.resize(sequence_count, simulator.getTimeSteps());
//	real_signal.resize(sequence_count, simulator.getTimeSteps());
//	imag_signal.resize(sequence_count, simulator.getTimeSteps());
//	real_error.resize(sequence_count, simulator.getTimeSteps());
//
//	for (unsigned s = 1; s < sequence_count; s++) {
//		unsigned period = s * periodResolution;
//		int currentRotationDirection = 1;
//		unsigned nextFlipTime = period / 2;
//		for (unsigned t = 0; t < simulator.getTimeSteps(); t++) {
//			if (t >= nextFlipTime) {
//				currentRotationDirection *= -1;
//				nextFlipTime += period;
//			}
//			sequences(s, t) = currentRotationDirection;
//		}
//	}
//}

/*
 * speichert eine kopie der uebergebenen sequenz und initialisiert die rotationsrichtungen.
 */
//template<typename S>
//void SimulationEvaluatorT<S>::addSpinEchoSequence(vector<double> sequence)
//{
//	if (sequence.size() == 0) throw MyException("invalid sequence");
//	sequence_count++;
//	sequences.resize(sequence_count, simulator.getTimeSteps());
//	real_signal.resize(sequence_count, simulator.getTimeSteps());
//	imag_signal.resize(sequence_count, simulator.getTimeSteps());
//	real_error.resize(sequence_count, simulator.getTimeSteps());
//
//	unsigned currentPulsIndex = 0;
//	int currentRotationDirection = 1;
//	double nextFlipTime = sequence[0];
//	for (unsigned t = 0; t < simulator.getTimeSteps(); t++) {
//		if (t * simulator.getTimeStepSize() >= nextFlipTime) {
//			if ((t * simulator.getTimeStepSize() - nextFlipTime) / simulator.getTimeStepSize() > 0.1) cout
//					<< "WARNING: Sequence puls position does not match timeresolution" << endl;
//			currentRotationDirection *= -1;
//			currentPulsIndex++;
//			double _nft = (sequence.size() > currentPulsIndex ? sequence[currentPulsIndex] : simulator.getTimeIntervall() * 10);
//			if (_nft - nextFlipTime < simulator.getTimeStepSize()) throw MyException("invalid sequence: sequence not in chronological order or to short pulse distance");
//			nextFlipTime = _nft;
//		}
//		sequences(sequence_count - 1, t) = currentRotationDirection;
//	}
//}

#endif /* SIMULATIONEVALUATORT_CPP_ */
