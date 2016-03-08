/*
 * Field2D.h
 *
 *  Created on: Aug 25, 2009
 *      Author: gstu0908
 */

#ifndef PRECALCULATEDFIELD2D_H_
#define PRECALCULATEDFIELD2D_H_

#include "MPICommunicator.h"
#include "Field2D.h"
//#include "IOManager.h"
#include <util.h>
#include <Point.h>
#include <Matrix2D.h>
#include <Rectangle.h>
#include <MyException.h>

#include <string>
#include <iostream>

using namespace std;

/*
 * Scalar Field, internally represented by a two dimensional Matrix.
 * The field is always assumed to be centered around the coordinate origin.
 * The y-dimension is in one memory block (See Matrix2D memory management).
 * For access of matrix elements use field(x,y) which is inherited from Matrix2D.
 * For fields values at arbitrary locations linear interpolation between matrix
 * indices is provided with field.getValue(x,y).
 * Maximal internal y index is matrix.columns.
 * Maximal internal x index is matrix.lines.
 */

class PrecalculatedField2D: public Field2D, public Rectangle, public Matrix2D
{
		friend istream& operator >>(istream &is, PrecalculatedField2D &obj);
		friend ostream& operator <<(ostream &os, const PrecalculatedField2D &obj);

	public:
		/*Create a field with 4 supporting points and size 1x1. the resolution of the field is then also 1.*/
		PrecalculatedField2D();
		PrecalculatedField2D(const PrecalculatedField2D &ref);
		/*Create a field with data from Matrix m and given size. the resolution of the field is determined by matrix size and geometric size.*/
		PrecalculatedField2D(const Matrix2D &ref, Point ll, Point ur);
		/*Create a not initialized field. resolution of the field is at least the target resolution*/
		PrecalculatedField2D(Point ll, Point ur, double target_resolution);

		virtual ~PrecalculatedField2D();

		/* Return value of the field by linear interpolating between grid points.*/
		double getValue(const double &x, const double &y) const;
		double getValue(const Point &p) const;

		/* Return value of the field by zero order interpolating between grid points.*/
		double& operator()(const Point &p);
		/* Forward direct matrix element access.*/
		double& operator()(const unsigned &x, const unsigned &y)		{return Matrix2D::operator()(x,y);}
		double  operator()(const unsigned &x, const unsigned &y) const	{return Matrix2D::operator()(x,y);}

		/* Transforms given set of unsigned array indices to real space position. */
		Point transform(const unsigned &x, const unsigned &y) const	;

		double getXResolution()	const {return resX;}
		double getYResolution()	const {return resY;}

		static PrecalculatedField2D readFile(string filename, Point ll, Point ur);

		PrecalculatedField2D &operator=(const PrecalculatedField2D &ref);

		/*Various operators.*/
		Matrix2D operator*(const double &d) 	const		{return PrecalculatedField2D(Matrix2D::operator*(d),ll,ur);}
		Matrix2D &operator*=(const double &d)				{Matrix2D::operator*=(d);	return (*this);}
		Matrix2D operator/(const double &d) 	const		{return PrecalculatedField2D(Matrix2D::operator/(d),ll,ur);}
		Matrix2D &operator/=(const double &d)				{Matrix2D::operator/=(d);	return (*this);}


	private:
		double resX, resY;

		void updateDescription();
};

#endif /* PRECALCULATEDFIELD2D_H_ */
