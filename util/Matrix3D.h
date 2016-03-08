/*
 * Field2D.h
 *
 *  Created on: Aug 25, 2009
 *      Author: gstu0908
 */

#ifndef MATRIX3D_H_
#define MATRIX3D_H_

#include "MyException.h"
#include <iostream>
#include "Matrix2D.h"

using namespace std;

class Matrix3D {
	friend istream& operator >>(istream &is, Matrix3D &obj);
	friend ostream& operator <<(ostream &os, const Matrix3D &obj);

public:
	Matrix3D(unsigned _slices, unsigned _lines, unsigned _columns);

	/* Copy Constructor needs to create local copy of data array to make a matrix element
	 * available as return type. If no local copy would be generated, statements like:
	 *
	 * Matrix3D bla()
	 * {
	 *	 	Matrix2D r(1,2,3);
	 * 		return r;
	 * }
	 *
	 * Would be illegal because r would go out of scope and data in r would be deleted.
	 * The pointer in the returned copy of r would still point to the deleted array.
	 */
	Matrix3D(const Matrix3D &ref);

	/* Access matrix data column and line*/
	double& operator()(unsigned s, unsigned l, unsigned c);
	double operator()(unsigned s, unsigned l, unsigned c) const ;

	Matrix3D &operator*=(const double &d);
	Matrix2D getSlice (unsigned index) const;
	Matrix2D getLine  (unsigned index) const;
	Matrix2D getColumn(unsigned index) const;

	virtual ~Matrix3D();

	double sum() const;
	Matrix3D &operator=(const Matrix3D &ref);
	unsigned columns, lines, slices;

private:
	/*
	 * Memory layout:
	 *   [(x0,y0),(x0,y1),(x0,y2),...,(x0,yi)]->
	 * ->[(x1,y0),(x1,y1),(x1,y2),...,(x1,yi)]->
	 * ->               ,...,				  ->
	 * ->[(xj,y0),(xj,y1),(xj,y2),...,(xj,yi)]
	 *
	 *then-> next slice!
	 * access to data should be performed via operator(x,y,z)!
	 */
	double *data;

};
#endif /* MATRIX3D_H_ */
