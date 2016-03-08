/*
 * SimulationEvaluator.h
 *
 *  Created on: Aug 18, 2009
 *      Author: gstu0908
 */

#ifndef SIMULATIONEVALUATOR3D_H_
#define SIMULATIONEVALUATOR3D_H_

#include <Matrix2D.h>
#include <Matrix3D.h>

using namespace std;


template<typename Simulator_Type>
class SimulationEvaluator3D
{
	public:
		SimulationEvaluator3D(const Simulator_Type &sim, unsigned samples);

		virtual ~SimulationEvaluator3D();

		/*Run random walks and save sequence data.*/
		void evaluateTrajectories();

		/*Samples the spatial distribution of magnetization by running LOTS of random walks....*/
//		void getMagnetizationDistribution(double spaceres, unsigned timeres);

		/*Add defined periodic multi spin echo sequences to evaluation process.*/
//		void preparePeriodicMultiSpinEcho(unsigned sequences, unsigned periodResolution);

		/*Add given sequence(s) to evaluation process.*/
//		void addSpinEchoSequence(vector<double> echos);
//		void addSpinEchoSequences(vector<vector<double> > sequences);

		/*Various getters.*/
		Simulator_Type const & getSimulator() 	const							{return simulator;}
		const Matrix2D &getRealSignal() 		const							{return real_signal;}
		const Matrix2D &getImagSignal() 		const							{return imag_signal;}
		const Matrix2D &getError()				const							{return real_error;}
		const Matrix2D &getSequences()			const							{return sequences;}
		unsigned 		getSamples()			const							{return samples;}
		unsigned 		getSequenceCount()		const							{return sequence_count;}
		const Matrix3D &getSpatialRealSignal() 	const							{return spatial_real;}
		const Matrix3D &getSpatialImagSignal() 	const							{return spatial_imag;}
		const Matrix3D &getSpatialError()		const							{return spatial_error;}
		const Matrix3D &getSpatialSamples()		const							{return spatial_samples;}

	private:
		const Simulator_Type &simulator;

		unsigned samples;

		/* Signal and error storage */
		unsigned sequence_count;
		Matrix2D real_signal, imag_signal;
		Matrix2D real_error;

		Matrix3D spatial_real, spatial_imag;
		Matrix3D spatial_error, spatial_samples;

		/* Sequence storage */
		Matrix2D sequences;

		unsigned nthreads;

};

#include "SimulationEvaluator3D.cpp"

#endif /* SIMULATIONEVALUATORT_H_ */
