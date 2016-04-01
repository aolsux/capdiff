/*
 * RandomNumberGenerator.h
 *
 *  Created on: 15.10.2010
 *      Author: Martin Rueckl
 */

#ifndef RANDOMNUMBERGENERATOR_H_
#define RANDOMNUMBERGENERATOR_H_

#include <random>

#include <float.h>
#include <limits.h>
#include <iostream>
#include <time.h>
#include <omp.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include "Point.h"



using namespace std;

class RandomNumberGenerator {
private:
	RandomNumberGenerator(const RandomNumberGenerator &ref);//noncopyable
	RandomNumberGenerator& operator=(const RandomNumberGenerator &ref);//noncopyable

public:
	/* Construct random number generator. seeds differently for every thread.*/
	RandomNumberGenerator(unsigned long seed = 1);

	virtual ~RandomNumberGenerator();

	/* generator for uniformly distributed doubles out of [0,1]*/
	double getUniform() const;

	/* random number generator for two normal distributed (mean=m, variance=s*s) doubles.
	 * Use "Point getNormalDistributed(double s, double m = 0) for using both numbers.*/
	double getNormalDistributed(const double &s,const double &m = 0) const;

	/* random number generator for two normal distributed (mean=m, variance=s*s) doubles.
	 * Returns a Point containing the two random numbers.*/
	Point getNormalDistributedPoint(const double &s, const double &m = 0) const;

//	template<class Point_Type>
//	Point_Type getNormalDistributedPointT(const double &s, const double &m = 0) const;

	unsigned long getInitialSeed() const;

private:
	unsigned int rank;
	unsigned long initialSeed;

	mutable std::normal_distribution<double> norm;
	mutable std::uniform_real_distribution<double> uniform;
	mutable std::ranlux48 gen;
};

//template<>
//Point RandomNumberGenerator::getNormalDistributedPointT<Point>(const double &s, const double &m) const
//{
//	return Point(norm_distrib_rnd()*s+m, norm_distrib_rnd()*s+m);
//}
//
//template<>
//Point3D RandomNumberGenerator::getNormalDistributedPointT<Point3D>(const double &s, const double &m) const
//{
//	return Point3D(norm_distrib_rnd()*s+m, norm_distrib_rnd()*s+m, norm_distrib_rnd()*s+m);
//}

#endif /* RANDOMNUMBERGENERATOR_H_ */
