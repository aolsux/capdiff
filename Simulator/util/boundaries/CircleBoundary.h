/*
 * CircleBoundary.h
 *
 *  Created on: 19.10.2010
 *      Author: Martin Rueckl
 */

#ifndef CIRCLEBOUNDARY_H_
#define CIRCLEBOUNDARY_H_

#include "OuterBoundary.h"
#include <util.h>
#include <MyException.h>
#include <Rectangle.h>

class CircleBoundary: public OuterBoundary
{
	public:
		CircleBoundary(double _radius) :
			 radius(_radius), radius_squared(_radius * _radius)
		{
			description="Circle Boundary (radius="+util::num2str(radius)+")";
		}

		Point collide(const Point &newLoc) const
		{
			Point ret = newLoc;
			if(newLoc.abs()>radius)	ret = newLoc * (2 * radius / newLoc.abs() - 1);
			return ret;
		}

		Point getRandomInsideLocation(const RandomNumberGenerator &rnd) const
		{
			Point p;
			do {
				p.x = radius * (2. * rnd.getUniform() - 1);
				p.y = radius * (2. * rnd.getUniform() - 1);
			} while (!isInside(p));
			return p;
		}

		bool handleCapillaryCollisions(CapillaryConfiguration &config) const
		{
			for (unsigned u = 0; u < config.capillaries.size(); u++)
				if (config.capillaries[u].getLocation().abs() + config.capillaries[u].getRadius() > radius) throw MyException("Invalid Circle Geometry: Capillary Overlay with Boundary!");
			return true;
		}

		bool 		isBoundaryCrossed(const Point &oldLoc, const Point &newLoc) const	{return(isInside(oldLoc)&&!isInside(newLoc) | !isInside(oldLoc)&&isInside(newLoc));		}
		bool 		isInside(const Point &p) 									const	{return p.absPow2() < radius_squared;}
		double 		getArea() 													const	{return radius_squared * M_PI;}
		Point 		getCenter() 												const	{return Point(0,0);}
		Rectangle 	getRectangle() 												const	{return Rectangle(2*radius,2*radius);}


	private:
		double radius, radius_squared;
};

#endif /* CIRCLEBOUNDARY_H_ */
