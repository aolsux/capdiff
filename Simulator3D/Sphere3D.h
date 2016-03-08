/*
 * Sphere.h
 *
 *  Created on: 11.04.2011
 *      Author: Martin Rueckl
 */

#ifndef SPHERE_H_
#define SPHERE_H_

#include "CollisionObject3D.h"
#include <MyException.h>
#include <Point.h>

class Sphere3D : public CollisionObject3D
{
	public:
		Sphere3D(Point3D _center, double _radius):
			center			(_center),
			radius			(_radius),
			radius_squared	(radius*radius)
		{}

		virtual Point3D collide(const Point3D &point) const
		{
			Point3D dif = point - center; //center->point
			Point3D buf = center + dif.normed() * 2 * radius - dif;
			return buf;
		}

		virtual bool	 isPlaneCrossed	(const Point3D &p, const Point3D &norm)	const {return (p-center)*norm<=radius;}
		virtual double	 getArea		()										const {return 4.*radius*radius*radius*M_PI/3.;}
		virtual bool	 isInside		(const Point3D &p) 						const {return (p-center).absPow2()<=radius_squared;}
		virtual Point3D  getCenter		()										const {return center;}
		virtual double	 getRadius		()										const {return radius;}
		virtual void	 translate		(const Point3D &p)							  {center=center+p;}

		virtual ~Sphere3D();

	private:
		Point3D center;
		double radius,radius_squared;
};

#endif /* SPHERE_H_ */
