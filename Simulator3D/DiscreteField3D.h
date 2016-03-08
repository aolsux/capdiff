/*
 * Field2D.h
 *
 *  Created on: Aug 25, 2009
 *      Author: gstu0908
 */

#ifndef DISCRETEFIELD3D_H_
#define DISCRETEFIELD3D_H_

#include <util.h>
#include <Point.h>
#include <Matrix3D.h>
#include <Rectangle.h>
#include <MyException.h>

#include <string>
#include <iostream>

using namespace std;

/*
 * Scalar Field, internally represented by a three dimensional Matrix.
 * The y-dimension is in one memory block (See Matrix2D memory management).
 * For access of matrix elements use field(x,y,z) which is inherited from Matrix3D.
 * For fields values at arbitrary locations linear interpolation between matrix
 * indices is provided with field.getValue(x,y,z).
 */

class DiscreteField3D: public Cuboid, public Matrix3D
{
		friend istream& operator >>(istream &is, DiscreteField3D &obj);
		friend ostream& operator <<(ostream &os, const DiscreteField3D &obj);

	public:
		/*Create a field with 4 supporting points and size 1x1. the resolution of the field is then also 1.*/
		DiscreteField3D();
		DiscreteField3D(const DiscreteField3D &ref);
		/*Create a field with data from Matrix m and given size. the resolution of the field is determined by matrix size and geometric size.*/
		DiscreteField3D(const Matrix3D &ref, Point3D llf, Point3D urb);
		/*Create a not initialized field. resolution of the field is at least the target resolution*/
		DiscreteField3D(Point3D llf, Point3D urb, double target_resolution);

		virtual ~DiscreteField3D();

		/* Return value of the field by linear interpolating between grid points.*/
		double getValue(const double &x, const double &y, const double &z) const;
		double getValue(const Point3D &p) const;

		/* Return value of the field by zero order interpolating between grid points.*/
		double& zeroOrderInterpolation(const Point3D &p);
		/* Forward direct matrix element access.*/
		double& operator()(const unsigned &x, const unsigned &y, const unsigned &z)			{return Matrix3D::operator()(x,y,z);}
		double  operator()(const unsigned &x, const unsigned &y, const unsigned &z) const	{return Matrix3D::operator()(x,y,z);}

		/* Transforms given set of unsigned array indices to real space position. */
		Point3D indexToPoint(const unsigned &x, const unsigned &y, const unsigned &z) const	;

		double getXResolution()	const {return resX;}
		double getYResolution()	const {return resY;}
		double getZResolution()	const {return resZ;}

		DiscreteField3D &operator=(const DiscreteField3D &ref);

		void transform(const Cuboid &cub);

		string getDescription() const;

//		/*Various operators.*/
//		Matrix2D operator*(const double &d) 	const		{return PrecalculatedField2D(Matrix2D::operator*(d),ll,ur);}
//		Matrix2D &operator*=(const double &d)				{Matrix2D::operator*=(d);	return (*this);}
//		Matrix2D operator/(const double &d) 	const		{return PrecalculatedField2D(Matrix2D::operator/(d),ll,ur);}
//		Matrix2D &operator/=(const double &d)				{Matrix2D::operator/=(d);	return (*this);}


	private:
		double resX, resY, resZ;
};

#endif /* PRECALCULATEDFIELD2D_H_ */
