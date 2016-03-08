/*
 * Field2D.h
 *
 *  Created on: Aug 25, 2009
 *      Author: gstu0908
 */

#ifndef MATRIX2D_H_
#define MATRIX2D_H_

#include "MyException.h"
#include <string>
#include <iostream>
#include <iomanip>

using namespace std;


/* Implementation of 2D matrices.
 * Inner loops should run on second index(columns)because this one is in a single memory block.
 * Memory layout:
 *   [(x0,y0),(x0,y1),(x0,y2),...,(x0,yi)]->
 * ->[(x1,y0),(x1,y1),(x1,y2),...,(x1,yi)]->
 * ->               ,...,				  ->
 * ->[(xj,y0),(xj,y1),(xj,y2),...,(xj,yi)]
 *
 * access to data should be performed via operator(x,y)!*/

class Matrix2D
{
		friend istream& operator >>(istream &is, Matrix2D &obj);
		friend ostream& operator <<(ostream &os, const Matrix2D &obj);

	public:

		Matrix2D(unsigned _lines, unsigned _columns);
		Matrix2D(const Matrix2D &ref);					/* Copy Constructor*/
		virtual ~Matrix2D();

		/* Access matrix data column and line.*/
		double& operator()(unsigned l, unsigned c);
		double operator()(unsigned l, unsigned c) const;

		/*Various operators.*/
		Matrix2D operator*(const double &d) const;
		Matrix2D &operator*=(const double &d);
		Matrix2D &operator+=(const Matrix2D &ref);
		Matrix2D operator/(const double &d) const		{return (*this)*(1./d);}
		Matrix2D &operator/=(const double &d)			{return (*this)*=(1./d);}
		Matrix2D &operator=(const Matrix2D &ref);		/* Assigns given reference to *this.*/

		double sum()			const;					/* Sum all elements.*/
		double getMaximum()		const;					/* Returns the largest element of the matrix.*/
		double getMinimum()		const;					/* Returns the smallest element of the matrix.*/
		Matrix2D transpose()	const;					/* Returns the transposed of *this.*/

		unsigned getLines()		const					{return lines;}
		unsigned getColumns()	const					{return columns;}

		unsigned lines, columns;

	protected:

		double *data;

};
#endif /* MATRIX2D_H_ */
