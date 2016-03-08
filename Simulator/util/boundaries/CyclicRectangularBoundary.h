/*
 * CyclicRectangularBoundary
 *
 *  Created on: 19.10.2010
 *      Author: Martin Rueckl
 */

#ifndef CYCLICRECTANGULARBOUNDARY_H_
#define CYCLICRECTANGULARBOUNDARY_H_

#include "RectangularBoundary.h"
#include <util.h>

class CyclicRectangularBoundary: public RectangularBoundary
{
	public:
		CyclicRectangularBoundary(Point _ll, Point _ur):RectangularBoundary(_ll,_ur)
		{
			description = "Cyclic " + description;
		}

		CyclicRectangularBoundary(double _width, double _height):RectangularBoundary(_width,_height)
		{
			description = "Cyclic " + description;
		}

		Point collide(const Point &newLoc) const
		{
			Point ret = newLoc - ll;
			if (ret.x <= 0) 				ret.x = ret.x + width;
			if (ret.x >= width) 			ret.x = ret.x - width;
			if (ret.y <= 0) 				ret.y = ret.y + height;
			if (ret.y >= height) 			ret.y = ret.y - height;
			return ret + ll;
		}
};

#endif /* CYCLICRECTANGULARBOUNDARY_H_ */
