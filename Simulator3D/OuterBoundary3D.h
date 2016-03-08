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

#include "CollisionObject3D.h"
#include "Configuration3D.h"

using namespace std;

class OuterBoundary: public CollisionObject3D, public Cuboid
{

public:

	OuterBoundary(const Cuboid &ref):
		Cuboid(ref)
	{}

	/* If given Point is outside the boundary returns mirrored Point. If no collision occurs, given Point is returned */
	Point3D collide(const Point3D &newLoc) const
	{
		Point ret = newLoc - ll;
		if (ret.x <= 0) 				ret.x = ret.x + width;
		if (ret.x >= width) 			ret.x = ret.x - width;
		if (ret.y <= 0) 				ret.y = ret.y + height;
		if (ret.y >= height) 			ret.y = ret.y - height;
		if (ret.z <= 0) 				ret.z = ret.z + depth;
		if (ret.z >= depth) 			ret.z = ret.z - depth;
		return ret + ll;
	}

	virtual Point3D getRandomInsideLocation(const RandomNumberGenerator &rnd) const
	{
		Point3D p;
		p.x = rnd.getUniform() * width;
		p.y = rnd.getUniform() * height;
		p.z = rnd.getUniform() * depth;
		return p + llf;
	}

	bool isBoundaryCrossed(const Point3D &oldLoc, const Point3D &newLoc) const
	{
		return ((isInside(oldLoc)&&!isInside(newLoc)) | (!isInside(oldLoc)&&isInside(newLoc)));
	}

	/* Checks whether there are collisions with the boundary and adds mirrored elements if needed	 */
	template <class Element_Type>
	bool handleCollisions(Configuration3D<Element_Type> &config) const
	{
		//define facet planes
		Point3D norm[6] = {Point3D(1,0,0),Point3D(0,1,0),Point3D(0,0,1),Point3D(-1,0,0),Point3D(0,-1,0),Point3D(0,0,-1)};
		Point3D p[6]  	= {llf		   	 ,llf			,llf		   ,urb			  ,urb		     ,urb			};
		vector<Element_Type>toAdd;
		for (unsigned i = 0; i < config.getElementCount(); i++) {
			const Element_Type &elem = config.elements[i];
			const Point3D &l = elem.getCenter();
			double r = elem.getRadius();
			bool flg_collision=false;
			for(unsigned f=0;f<6;f++)//handle all facets
			{
				if(elem.isPlaneCrossed(p[f],norm[f])){
					Element_Type copy = elem;
					copy.translate(Point3D(norm[f].x*width, norm[f].y*height, norm[f].z*depth));
					toAdd.push_back(copy);
				}
			}

		}
		config.addCapillaries(toAdd);
		return toAdd.size() > 0;
	}

};

#endif /* OUTERBOUNDARY_H_ */
