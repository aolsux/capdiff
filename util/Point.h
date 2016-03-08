/*
 * Point.h
 *
 *  Created on: 15.10.2010
 *      Author: Martin Rueckl
 */

#ifndef POINT_H_
#define POINT_H_

#include <math.h>

class Point
{
		//friend double operator*(const Point &p1, const Point &p2);
		//friend double operator*(Point p1, Point p2);

	public:
		double x, y;

		Point():x(0), y(0){}
		Point(double _x, double _y):x(_x), y(_y){}

		Point operator+(const Point &p2)	const	{return Point(x + p2.x, y + p2.y);}
		Point operator-(const Point &p2)	const	{return Point(x - p2.x, y - p2.y);}
		Point operator*(const double &d)	const	{return Point(x * d, y * d);}
		Point operator*(const int &d) 		const	{return Point(x * d, y * d);}
		Point operator/(const double &d)	const	{return Point(x / d, y / d);}
		double operator*(const Point &ref)	const	{return x * ref.x + y * ref.y;}

		bool operator==(const Point &p2)	const	{return x == p2.x && y == p2.y;}
		bool operator!=(const Point &p2)	const	{return x != p2.x || y != p2.y;}

		Point &operator+=(const Point &p2)			{x += p2.x;	y += p2.y; return (*this);}
		Point &operator-=(const Point &p2)			{x -= p2.x;	y -= p2.y; return (*this);}

		Point normed()						const	{double b = abs();return Point(x / b, y / b);}
		double abs()						const	{return sqrt(x * x + y * y);}
		double absPow2()					const	{return x * x + y * y;}
};

class Point3D
{

	public:
		double x, y, z;

		Point3D():x(0), y(0), z(0){}
		Point3D(double _x, double _y, double _z):x(_x), y(_y), z(_z){}

		Point3D operator+(const Point3D &p2) const	{return Point3D(x + p2.x, y + p2.y, z + p2.z);}
		Point3D operator-(const Point3D &p2) const	{return Point3D(x - p2.x, y - p2.y, z - p2.z);}
		Point3D operator*(const double &d)	 const	{return Point3D(x * d, y * d, z * d);}
		Point3D operator*(const int &d)		 const	{return Point3D(x * d, y * d, z * d);}
		Point3D operator/(const double &d)	 const	{return Point3D(x / d, y / d, z / d);}
		double operator*(const Point3D &ref) const	{return x * ref.x + y * ref.y + z * ref.z;}

		Point3D &operator+=(const Point3D &p2)		{x += p2.x;	y += p2.y;	z += p2.z;	return (*this);	}
		Point3D &operator-=(const Point3D &p2)		{x -= p2.x;	y -= p2.y;	z -= p2.z;	return (*this);	}

		Point3D normed()					 const	{double b = abs();return Point3D(x / b, y / b, z / b);}
		double abs()						 const	{return sqrt(x * x + y * y + z * z);}
		double absPow2()					 const	{return x * x + y * y + z * z;}
		Point3D cross(const Point3D &ref)	 const	{return Point3D(y*ref.z - z*ref.y,
																	z*ref.x - x*ref.z,
																	x*ref.y - y*ref.x);}
};

#endif /* POINT_H_ */
