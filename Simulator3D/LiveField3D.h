/*
 * LiveField2D.h
 *
 *  Created on: 15.10.2010
 *      Author: Martin Rueckl
 */

#ifndef LIVEFIELD3D_H_
#define LIVEFIELD3D_H_

#include <Point.h>
#include "Configuration3D.h"

class LiveField3D
{
	friend ostream& operator <<(ostream &os, const LiveField3D &obj);

public:
	LiveField3D(SphereConfiguration config)										{	capillaryConfiguration = config;	}

	double getValue(const double &x, const double &y, const double &z) const	{	return getValue(Point3D(x, y, z));	}

	string getDescription() const												{	return  "LiveField3D";	}

	double getValue(const Point3D &p) const
	{
		double dOmega = 0;
		//iterate over all elements in the system
		for (unsigned i = 0; i < capillaryConfiguration.elements.size(); i++) {
			double radius = capillaryConfiguration.elements[i].getRadius();
			Point3D dif = p - capillaryConfiguration.elements[i].getCenter();
			//R_c^3 (3cos^2(\theta) -1)/r^3
			dOmega += radius * radius * radius * (3.*(dif.x * dif.x + dif.y * dif.y)/dif.absPow2() - 1.) / (dif.absPow2() * dif.abs());
		}
		return dOmega;
	}

private:

	SphereConfiguration capillaryConfiguration;
};

ostream& operator <<(ostream &os, const LiveField3D &obj) {
	os << obj.capillaryConfiguration;
	return os;
}

#endif /* LIVEFIELD3D_H_ */
