/*
 * Capillary2D.h
 *
 *  Created on: Aug 12, 2009
 *      Author: gstu0908
 */

#ifndef CAPILLARY2D_H_
#define CAPILLARY2D_H_

#include "CollisionObject.h"
#include <Point.h>

using namespace std;

class Capillary2D: public CollisionObject {
	friend class CapillaryConfiguration;
public:
	Capillary2D(Point loc, double radius, double degrees = 0);

	bool isBoundaryCrossed(const Point &oldLoc, const Point &newLoc) const;
	Point collide(const Point &newLoc) const;

	double getArea() const;
	bool isInside(const Point &p) const;

	const Point &getLocation()	const {return location;}
	const double &getAngle()	const {return angle;}
	const Point &_getAngle()	const {return _angle;}
	const double &getRadius()	const {return radius;}
	Point getCenter()			const {return location;}

private:

	inline bool isDoubleCrossed(const Point &oldLoc, const Point &newLoc) const;

	Point location;

	//stores angle in degrees
	double angle;
	//stores angle as (cos(phi), sin(phi))
	Point _angle;
	double radius, radius_squared;

};

#endif /* CAPILLARY2D_H_ */
