/*
 * CollisionObject.h
 *
 *  Created on: Aug 13, 2009
 *      Author: gstu0908
 */

#ifndef COLLISIONOBJECT_H_
#define COLLISIONOBJECT_H_

#include <Point.h>

class CollisionObject {
public:
	/*
	 *checks wether boundary is crossed
	 */
	virtual bool isBoundaryCrossed(const Point &oldLoc, const Point &newLoc) const=0;

	/*
	 * checks wether oldLoc and newLoc are seperated by the objects boundary.
	 * if so returns true and writes the the new location after collision into buffer.
	 * if no collision occures, buffer will contain newLoc
	 */
	virtual Point collide(const Point &newLoc) const=0;

	/*
	 * returns the area of the object
	 */
	virtual double getArea() const =0;

	/*
	 * return wether the selected point is inside the area of the object
	 */
	virtual bool isInside(const Point &p) const=0;

	virtual Point getCenter()const =0;
};

#endif /* COLLISIONOBJECT_H_ */
