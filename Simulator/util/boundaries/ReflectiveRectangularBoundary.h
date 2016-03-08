/*
 * ReflectiveRectangularBoundary.h
 *
 *  Created on: 19.10.2010
 *      Author: Martin Rueckl
 */

#ifndef REFLECTIVERECTANGULARBOUNDARY_H_
#define REFLECTIVERECTANGULARBOUNDARY_H_

#include "RectangularBoundary.h"
#include <util.h>

class ReflectiveRectangularBoundary: public RectangularBoundary
{
	public:
		ReflectiveRectangularBoundary(Point _ll, Point _ur):RectangularBoundary(_ll,_ur)
		{
			description = "Reflective " + description;
		}

		ReflectiveRectangularBoundary(double _width, double _height):RectangularBoundary(_width,_height)
		{
			description = "Cyclic " + description;
		}

		Point collide(const Point &newLoc) const
		{
			Point ret = newLoc - ll;
			if (ret.x <= 0) 				ret.x = -ret.x;
			if (ret.x >= width) 			ret.x = width - (ret.x - width);
			if (ret.y <= 0) 				ret.y = -ret.y;
			if (ret.y >= height) 			ret.y = height - (ret.y - height);
			return ret + ll;
		}
};

#endif /* RECTANGULARBOUNDARY_H_ */
