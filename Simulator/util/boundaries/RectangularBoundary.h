/*
 * RectangularBoundary.h
 *
 *  Created on: 19.10.2010
 *      Author: Martin Rueckl
 */

#ifndef RECTANGULARBOUNDARY_H_
#define RECTANGULARBOUNDARY_H_

#include "OuterBoundary.h"
#include <util.h>

class RectangularBoundary: public OuterBoundary, public Rectangle
{
	public:
		/* Abstract Methods to be implemented by Reflective or Cyclic Boundaries */
		virtual Point collide(const Point &p) const =0;

		/* ======   Base Class Method implementations   ====== */
		RectangularBoundary(double _width, double _height) :
			Rectangle(_width,_height)
		{
			description="Rectangular Boundary (ll="+util::num2str(ll.x)+"/"+util::num2str(ll.y)+", width="+util::num2str(width)+", height="+util::num2str(height)+")";
		}

		RectangularBoundary(Point _ll, Point _ur) :
			Rectangle(_ll,_ur)
		{
			description="Rectangular Boundary (ll="+util::num2str(ll.x)+"/"+util::num2str(ll.y)+", width="+util::num2str(width)+", height="+util::num2str(height)+")";
		}

		virtual Point getRandomInsideLocation(const RandomNumberGenerator &rnd) const
		{
			Point p;
			p.x = rnd.getUniform() * width;
			p.y = rnd.getUniform() * height;
			return p + ll;
		}

		bool isBoundaryCrossed(const Point &oldLoc, const Point &newLoc) const
		{
			return ((Rectangle::isInside(oldLoc)&&!Rectangle::isInside(newLoc)) | (!Rectangle::isInside(oldLoc)&&Rectangle::isInside(newLoc)));
		}

		bool 		isInside(const Point &p)	const	{return Rectangle::isInside(p);}
		double 		getArea()					const	{return Rectangle::getArea();}
		Point		getCenter()					const	{return Rectangle::getCenter();}
		Rectangle	getRectangle()				const	{return Rectangle(*this);}

		bool handleCapillaryCollisions(CapillaryConfiguration &config) const
		{
			vector<Capillary2D> toAdd;

			string warning = "Warning: Capillary overlap while mirroring boundary capillaries. Assuming no mirroring is needed...\n";
			for (unsigned i = 0; i < config.getCapillaryCount(); i++) {
				const Capillary2D &cap = config.capillaries[i];
				const Point &l = cap.getLocation();
				double r = cap.getRadius();
				double angle = cap.getAngle();

				bool flg_tl = cap.isInside(getUpperLeftCorner());//Point(-width_half, height_half));
				bool flg_tr = cap.isInside(getUpperRightCorner());//Point(width_half, height_half));
				bool flg_ll = cap.isInside(getLowerLeftCorner());//Point(-width_half, -height_half));
				bool flg_lr = cap.isInside(getLowerRightCorner());//Point(width_half, -height_half));

				if (flg_tl | flg_tr | flg_ll | flg_lr) {//capillary is in one of the corners
					Point m1, m2, m3;

					if (flg_tl) {
						m1 = l + Point(width, 0);
						m2 = l + Point(width, -height);
						m3 = l + Point(0, -height);
					} else if (flg_tr) {
						m1 = l + Point(-width, 0);
						m2 = l + Point(-width, -height);
						m3 = l + Point(0, -height);
					} else if (flg_ll) {
						m1 = l + Point(width, 0);
						m2 = l + Point(width, height);
						m3 = l + Point(0, height);
					} else if (flg_lr) {
						m1 = l + Point(-width, 0);
						m2 = l + Point(-width, height);
						m3 = l + Point(0, height);
					}
					if (config.isInside(m1)) cout << warning;
					else toAdd.push_back(Capillary2D(m1, r, angle));
					if (config.isInside(m2)) cout << warning;
					else toAdd.push_back(Capillary2D(m2, r, angle));
					if (config.isInside(m3)) cout << warning;
					else toAdd.push_back(Capillary2D(m3, r, angle));

				} else {//capillary might be on one of the edges
					if ((l.x - r <= ll.x) && (ll.x <= l.x + r)) { //crossing with left boundary -> mirror to right side
						Point ml = l - Point(width, 0);
						if (config.isInside(ml)) cout << warning;
						else toAdd.push_back(Capillary2D(ml, r, angle));
					}
					if ((l.x - r <= ur.x) && (ur.x <= l.x + r)) {//crossing with right boundary -> mirror to left side
						Point ml = l + Point(width, 0);
						if (config.isInside(ml)) cout << warning;
						else toAdd.push_back(Capillary2D(ml, r, angle));
					}
					if ((l.y - r <= ll.y) && (ll.y <= l.y + r)) { //crossing with lower boundary -> mirror to upper side
						Point ml = l + Point(0, height);
						if (config.isInside(ml)) cout << warning;
						else toAdd.push_back(Capillary2D(ml, r, angle));
					}
					if ((l.y - r <= ur.y) && (ur.y <= l.y + r)) { //crossing with upper boundary -> mirror to lower side
						Point ml = l - Point(0, height);
						if (config.isInside(ml)) cout << warning;
						else toAdd.push_back(Capillary2D(ml, r, angle));
					}
				}
			}
			config.addCapillaries(toAdd);
			return toAdd.size() > 0;
		}
};

#endif /* RECTANGULARBOUNDARY_H_ */
