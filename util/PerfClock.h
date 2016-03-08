/*
 * PerfClock.h
 *
 *  Created on: 19.10.2009
 *      Author: Administrator
 */

#ifndef PERFCLOCK_H_
#define PERFCLOCK_H_

#include <time.h>
#include <iostream>

/*
 *
 */
class PerfClock {
public:
	PerfClock();
	virtual ~PerfClock();
	static void tick() {
		start = (double) clock() / CLOCKS_PER_SEC;
	}

	static double tack() {
		stop = (double) clock() / CLOCKS_PER_SEC;
		return stop - start;
	}

	static void printTack() {
		std::cout << stop - start << "s" << std::endl;
	}

	static double dt() {
		return stop - start;
	}

protected:
	static double start;
	static double stop;

};

#endif /* PERFCLOCK_H_ */
