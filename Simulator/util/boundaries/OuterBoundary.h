/*
 * OuterReflectiveCirleBondary.h
 *
 *  Created on: Aug 14, 2009
 *      Author: gstu0908
 */

#ifndef OUTERBOUNDARY_H_
#define OUTERBOUNDARY_H_

#include <Point.h>
#include <Rectangle.h>
#include <MyException.h>
#include <RandomNumberGenerator.h>

#include "../Describable.h"
#include "../CapillaryConfiguration.h"
#include "../CollisionObject.h"

using namespace std;

class OuterBoundary: public Describable, public CollisionObject
{
	public:

		/*
		 * if given Point is outside the boundary returns mirrored Point.
		 * if no collision occurs, given Point is returned
		 */
		virtual Point collide(const Point &p) const =0;

		/*
		 * returns the area of the object
		 */
		virtual double getArea() const=0;

		/*
		 * return whether the selected point is inside the area of the object
		 */
		virtual bool isInside(const Point &p) const=0;

		virtual Point getRandomInsideLocation(const RandomNumberGenerator &rnd) const=0;

		/*
		 * Checks whether there are capillaries at the edge of boundary and adds mirrored capillaries if needed
		 */
		virtual bool handleCapillaryCollisions(CapillaryConfiguration &config) const =0;

		/*
		 * Returns the smallest Rectangle containing the whole Area.
		 */
		virtual Rectangle getRectangle() const =0;

};

#endif /* OUTERBOUNDARY_H_ */
