/*
 * Field2D.h
 *
 *  Created on: Aug 25, 2009
 *      Author: gstu0908
 */

#ifndef MATRIX4D_H_
#define MATRIX4D_H_

#include "MyException.h"
#include <string>
#include <iostream>
#include <iomanip>
#include "Matrix3D.h"

using namespace std;

class Matrix4D
{
		friend istream& operator >>(istream &is, Matrix4D &obj);
		friend ostream& operator <<(ostream &os, const Matrix4D &obj);

	public:

		Matrix4D(unsigned i1, unsigned i2, unsigned i3, unsigned i4);
		Matrix4D(const Matrix4D &ref);					/* Copy Constructor*/
		virtual ~Matrix4D();

		/* Access matrix data column and line.*/
		double& operator()(unsigned i1, unsigned i2, unsigned i3, unsigned i4);
		double operator()(unsigned i1, unsigned i2, unsigned i3, unsigned i4) const;

		Matrix3D getSubMatrix(unsigned dimension, unsigned index) const;

		/*Various operators.*/
		Matrix4D operator*(const double &d) const;
		Matrix4D &operator*=(const double &d);
		Matrix4D &operator=(const Matrix4D &ref);		/* Assigns given reference to *this. */
		double sum() const;								/* Sum all elements. */

		unsigned size() const;							/* Returns amount of data elements in Matrix. */

		unsigned i1,i2,i3,i4;

	protected:

		double *data;

};
#endif /* MATRIX4D_H_ */
