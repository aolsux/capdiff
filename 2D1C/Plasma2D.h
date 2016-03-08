/*
 * Plasma2D.h
 *
 *  Created on: 23.11.2010
 *      Author: Martin Rueckl
 */

#ifndef PLASMA2D_H_
#define PLASMA2D_H_

#include <Point.h>
#include <vector>
#include <Matrix2D.h>
#include <Matrix3D.h>
#include <iostream>
#include <RandomNumberGenerator.h>
//#include <boost/threadpool/pool.hpp>

using namespace std;

class Charge
{
	public:
		Charge(double x, double y, double c) :
			location(x, y), charge(c)
		{
		}

		Point location;
		double charge;
};

class Plasma2D
{
		friend ostream& operator <<(ostream &os, const Plasma2D &obj);

	public:
		Plasma2D(vector<Charge> charges, double basis_length);
		Plasma2D(unsigned charge_count, double basis_length, double charge);

		virtual ~Plasma2D();

		void modify();
		void undo();

		double getEnergy() const;

		const Matrix2D& getDistanceMatrix() const;

		Matrix2D getPairCorrelation(double dr, double min = 0, double max = 20) const;

		unsigned getChargeCount()const{return charge_count;}

	private:
		double expint1_NRIC(double x) const;

		/* calculates the energy contibuted by particle i */
		void calculateParticleEnergy(unsigned i);

		/* Calculates total energy by summing over real and k-space matrices */
		double calculateEnergy();

		double sumKSpaceMatrices() const;
		double sumRealSpaceMatrix() const;

		void updateDistanceMatrix(unsigned j);
		bool collision(unsigned i,double distance);

	private:
		double energy;
		vector<Charge> charges;
		unsigned charge_count;
		double total_charge;
		double base_width;
		double base_height;
		double kSpace_base_width;
		double kSpace_base_height;
		int kSpace_cutoff;
		double eta_squared;/* "magical" convergence parameter taken from paper... */

		//store energy distributions of each particle
		Matrix2D realSpaceEnergyMatrix;/* obere dreiecksmatrix: enthält die realraum energieterme der partikel. gesammtmatrix enthält alles doppelt*/
		Matrix3D kSpaceEnergyMatrix_real;
		Matrix3D kSpaceEnergyMatrix_imag;

		//keep track of distance matrix for eventually doing collision detection
		Matrix2D distance_matrix;

		//keep track of changed particle
		unsigned modifiedParticle;
		Point oldLocation;
		double oldEnergy;

		RandomNumberGenerator rng;
};

#endif /* PLASMA2D_H_ */

