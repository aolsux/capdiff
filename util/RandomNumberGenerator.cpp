/*
 * RandomNumberGenerator.h
 *
 *  Created on: 15.10.2010
 *      Author: Martin Rueckl
 */

#include "RandomNumberGenerator.h"

/* auxiliary variables for the rng's thx @niclas*/
//unsigned long RandomNumberGenerator::x, RandomNumberGenerator::y, RandomNumberGenerator::z, RandomNumberGenerator::v,
//		RandomNumberGenerator::w;
//unsigned long RandomNumberGenerator::initialSeed = 0;




RandomNumberGenerator::RandomNumberGenerator(unsigned long seed):
	rank					(omp_get_thread_num()),
	initialSeed				((seed == 1) ? ((unsigned long) (time(NULL)))*(rank+1) : seed),
	gen						(initialSeed)
{
	cout<< "rng's rank:\t"<<rank<<"\tseed:"<<initialSeed<<endl;
}

RandomNumberGenerator::~RandomNumberGenerator(){}

unsigned long RandomNumberGenerator::getInitialSeed() const {
	return initialSeed;
}

double RandomNumberGenerator::getNormalDistributed(const double &s, const double &m) const {
	return norm(gen)*s+m;
}

Point RandomNumberGenerator::getNormalDistributedPoint(const double &s, const double &m) const {
	return Point(norm(gen)*s+m, norm(gen)*s+m);
}

double RandomNumberGenerator::getUniform() const {
	return uniform(gen);
}
