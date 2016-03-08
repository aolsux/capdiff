/*
 * MixedBoundary.h
 *
 *  Created on: 17.02.2011
 *      Author: Martin Rueckl
 */

#include "OuterBoundary.h"
#include "CircleBoundary.h"
#include "RectangularBoundary.h"
#include <MyException.h>

#ifndef MIXEDBOUNDARY_H_
#define MIXEDBOUNDARY_H_

class MixedBoundary : public OuterBoundary
{
	public:
		MixedBoundary(double _radius) :
			radius				(_radius),
			radius_squared		(_radius*_radius)
		{
		}


		Point collide(const Point &newLoc) const
		{
			Point ret = newLoc;
			if (ret.absPow2() > radius_squared)		ret	  = newLoc * (2 * radius / newLoc.abs() - 1);
			if (ret.x < 0)							ret.x = -ret.x ;
			if (ret.y < 0)							ret.y = -ret.y;
			return ret;
		}

		double getArea() const
		{
			return M_PI * radius_squared/4.;
		}

		bool isInside(const Point &p) const
		{
			return (p.absPow2() < radius_squared) && (p.x > 0) && (p.y > 0);
		}

		Point getRandomInsideLocation(const RandomNumberGenerator &rnd) const
		{
			Point p;
			do {
				p.x = radius * rnd.getUniform();
				p.y = radius * rnd.getUniform();
			} while (!isInside(p));
			return p;
		}

		string printShapeInfo() const
		{
			return "MixedBoundary[1.Quadrant](radius=" + util::num2str(radius)  + ")";
		}

		bool handleCapillaryCollisions(CapillaryConfiguration &config) const
		{
			for (unsigned u = 0; u < config.capillaries.size(); u++)
				if (config.capillaries[u].getLocation().abs() + config.capillaries[u].getRadius() > radius) throw MyException("Invalid Circle Geometry: Capillary Overlay with Boundary!");
			return true;
		}

	private:
		double radius, radius_squared;
};

#endif /* MIXEDBOUNDARY_H_ */
