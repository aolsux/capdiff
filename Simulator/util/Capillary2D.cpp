/*
 * Capillary2D.cpp
 *
 *  Created on: Aug 12, 2009
 *      Author: gstu0908
 */

#include "Capillary2D.h"
#include <math.h>

Capillary2D::Capillary2D(Point loc, double _radius, double degrees) {
	location = loc;
	this->radius = _radius;
	this->radius_squared = _radius * _radius;
	this->angle = degrees;
	this->_angle = Point(cos((double) degrees * M_PI / 180.), sin((double) degrees * M_PI / 180.));
}

bool Capillary2D::isBoundaryCrossed(const Point &oldLoc, const Point &newLoc) const {
	return isInside(newLoc) | isDoubleCrossed(oldLoc, newLoc);
}

Point Capillary2D::collide(const Point & newLoc) const {
	//Point s = location + (oldLoc - location).normed() * radius;//schnittpunkt
	//if ((s - newLoc).abs() >= radius) {
	//throw MeanStepSizeException("Stepped to far into Capillary!", (oldLoc - newLoc).abs());
	Point dif = newLoc - location;
	Point buf = location + dif.normed() * 2 * radius - dif;
	return buf;
}

double Capillary2D::getArea() const {
	return radius_squared * M_PI;
}

bool Capillary2D::isDoubleCrossed(const Point &oldLoc, const Point &newLoc) const {
	/*
	 * verbindungsgerade wird parametrisiert und abstand kreismittelpunkt zu param. geradenpunkt gleich kreisradius gesetzt
	 * aus diskriminante folgen dann anzahl schnittpunkte, aus der parametrisierung der/die schnittpunkt/e
	 */
	double a = oldLoc * oldLoc + newLoc * newLoc - 2 * (oldLoc * newLoc);
	double b = -2 * (location * newLoc) + 2 * (location * oldLoc) + 2 * (newLoc * oldLoc) - 2 * (oldLoc * oldLoc);
	double c = location * location - 2 * (location * oldLoc) + oldLoc * oldLoc - radius_squared;
	double d = b * b - 4 * a * c;//Diskriminante
	if (d <= 0)//ein oder kein schnittpunkt
		return false;
	d = sqrt(d);
	a *= 2;
	double t1 = (-b + d) / a;
	double t2 = (-b - d) / a;
	if (0 <= t1 && t1 <= 1 && 0 <= t2 && t2 <= 1)//beide in [0,1]->zwei schnittpunkte
		return true;
	return false;
}

bool Capillary2D::isInside(const Point &p) const {
	return (p - location).absPow2() < radius_squared;
}

