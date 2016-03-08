/*
 * LiveField2D.h
 *
 *  Created on: 15.10.2010
 *      Author: Martin Rueckl
 */

#ifndef LIVEFIELD2D_H_
#define LIVEFIELD2D_H_

#include "Field2D.h"
#include <Point.h>
#include "CapillaryConfiguration.h"

class LiveField2D: public Field2D
{
		friend ostream& operator <<(ostream &os, const LiveField2D &obj);

public:
	LiveField2D(CapillaryConfiguration config) {
		capillaryConfiguration = config;
		description = "Live-Field";
	}

	double getValue(const double &x, const double &y) const {
		Point p(x, y);
		return getValue(p);
	}

	double getValue(const Point &p) const {
		double dOmega = 0;
		//iterate over all capillaries in the system
		for (unsigned i = 0; i < capillaryConfiguration.capillaries.size(); i++) {
			double radius = capillaryConfiguration.capillaries[i].getRadius();
			Point dif = p - capillaryConfiguration.capillaries[i].getLocation();
			dOmega += radius * radius * (dif.x * dif.x - dif.y * dif.y) / (dif.absPow2() * dif.absPow2());
		}
		return dOmega;
	}

private:
	CapillaryConfiguration capillaryConfiguration;
};

ostream& operator <<(ostream &os, const LiveField2D &obj) {
	os << obj.capillaryConfiguration;
	return os;
}

#endif /* LIVEFIELD2D_H_ */
