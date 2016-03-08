/*
 * CollisionObject.h
 *
 *  Created on: Aug 13, 2009
 *      Author: gstu0908
 */

#ifndef COLLISIONOBJECTT_H_
#define COLLISIONOBJECTT_H_

#include <Point.h>
#include <Rectangle.h>

class CollisionObject3D {
public:

	/*checks whether boundary is crossed*/
	virtual bool isPlaneCrossed(const Point3D &p, const Point3D &norm) const=0;
	virtual bool isPlaneCrossed(const Point3D &p, const Point3D &v1, const Point3D &v2) const {return isPlaneCrossed(p,v1.cross(v2));}


	/*reflects point to corresponding location outside the object*/
	virtual Point3D collide(const Point3D &point) const=0;

	/* returns the area of the object */
	virtual double getArea() const =0;

	/* return whether the selected point is inside the area of the object */
	virtual bool isInside(const Point3D &p) const=0;

	virtual Point3D getCenter()const =0;

	/* center and radius define a circle which completely contains the object... */
	virtual double getRadius()const =0;

	virtual void translate(const Point3D &p) =0;
};

#endif /* COLLISIONOBJECTT_H_ */
