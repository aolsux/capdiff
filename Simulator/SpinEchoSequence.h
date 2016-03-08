/*
 * SpinEchoSequence.h
 *
 * in dieser klasse koennen spinechosequenzen initialisiert werden.
 * dasu werden einfach die zeitpunkte der pulse in den vektor eingefuegt.
 * timebase
 *
 *  Created on: Aug 19, 2009
 *      Author: gstu0908
 */

#ifndef SPINECHOSEQUENCE_H_
#define SPINECHOSEQUENCE_H_

#include <vector>

using namespace std;

class SpinEchoSequence: public vector<double> {
public:

	SpinEchoSequence(double _timeBase = 1.0) :
		timeBase(_timeBase) {
	}

	SpinEchoSequence(vector<double> _pulses, double _timeBase = 1.0) :
		vector<double> (_pulses), timeBase(_timeBase) {
	}

	double getTimeBase() const {
		return timeBase;
	}
private:
	double timeBase;

};

#endif /* SPINECHOSEQUENCE_H_ */
