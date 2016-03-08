/*
 * Matrix4D.cpp
 *
 *  Created on: 14.03.2011
 *      Author: Martin Rueckl
 */
#include "Matrix4D.h"

#include <iostream>
#include <sstream>
#include <stdio.h>
#include <fstream>
#include <algorithm>


Matrix4D::Matrix4D(unsigned _i1, unsigned _i2, unsigned _i3, unsigned _i4):
	i1(_i1), i2(_i2), i3(_i3), i4(_i4)
{
	data = new double[size()];
	for (unsigned i = 0; i < size(); i++)
		data[i] = 0;
}

/* Copy Constructor needs to create local copy of data array to make a matrix element
 * available as return type. If no local copy would be generated, statements like:
 *
 * Matrix2D bla()
 * {
 *	 	Matrix2D r(1,2);
 * 		return r;
 * }
 *
 * Would be illegal because r would go out of scope and data in r would be deleted.
 * The pointer in the returned copy of r would still point to the deleted array.
 */
Matrix4D::Matrix4D(const Matrix4D &ref) :
	i1(ref.i1),	i2(ref.i2),	i3(ref.i3),	i4(ref.i4)
{
	data = new double[size()];
	for (unsigned i = 0; i < size(); i++)
		data[i] = ref.data[i];
}

/*
 * access matrix data column and line
 */
double& Matrix4D::operator()(unsigned _i1, unsigned _i2, unsigned _i3, unsigned _i4)
{
	if (_i1 >= i1 || _i2 >= i2 || _i3 >= i3 || _i4 >= i4) {
		cerr << "i1=" << _i1 << " i2=" << _i2 << " i3=" << _i3 << " i4=" << _i4 << endl;
		throw MyException("Matrix subscript out of bounds");
	}
	return data[_i1*i2*i3*i4 + _i2*i3*i4 + _i3*i4 + _i4];
}

double Matrix4D::operator()(unsigned _i1, unsigned _i2, unsigned _i3, unsigned _i4) const
{
	if (_i1 >= i1 || _i2 >= i2 || _i3 >= i3 || _i4 >= i4) {
		cerr << "i1=" << _i1 << " i2=" << _i2 << " i3=" << _i3 << " i4=" << _i4 << endl;
		throw MyException("Matrix subscript out of bounds");
	}
	return data[_i1*i2*i3*i4 + _i2*i3*i4 + _i3*i4 + _i4];
}

Matrix4D &Matrix4D::operator*=(const double &d)
{
	for (unsigned i = 0; i < size(); i++)
		data[i] *= d;
	return *this;
}

Matrix4D Matrix4D::operator*(const double &d) const {
	Matrix4D m(*this);
	m*=d;
	return m;
}

Matrix4D::~Matrix4D()
{
	delete data;
}

double Matrix4D::sum() const
{
	double ans = 0;
	for (unsigned i = 0; i < size(); i++)
		ans += data[i];
	return ans;
}

/*
 * Assigns given reference to *this.
 */
Matrix4D &Matrix4D::operator=(const Matrix4D &ref)
{
	//check for self assignment
	if (this == &ref) return *this;

	delete data;
	data = new double[ref.size()];

	i1= ref.i1;i2 = ref.i2;i3 = ref.i3;i4 = ref.i4;

	for (unsigned i = 0; i < size(); i++)
		data[i] = ref.data[i];

	return *this;
}

unsigned Matrix4D::size() const
{
	return i1*i2*i3*i4;
}

Matrix3D getSubMatrix(unsigned dimension, unsigned index) const
{
	unsigned m3s1, m3s2, m3s3;
	switch(dimension){
		case 1:
			m3s1=i2;
			m3s2=i3;
			m3s3=i4;
			break;
		case 2:
			m3s1=i1;
			m3s2=i3;
			m3s3=i4;
			break;
		case 3:
			m3s1=i1;
			m3s2=i2;
			m3s3=i4;
			break;
		case 4:
			m3s1=i1;
			m3s2=i2;
			m3s3=i3;
			break;
		default:
			throw MyException("invalid dimension selection...");
	}
	Matrix3D ret(m3s1,m3s2,m3s3);
	double elem;
	for(unsigned c1=0;c1<m3s1;c1++)
		for(unsigned c2=0;c2<m3s2;c2++)
			for(unsigned c3=0;c3<m3s3;c3++)
			{
			switch(dimension){
				case 1:
					ret(c1,c2,c3)= (*this)(index,c1,c2,c3);
					break;
				case 2:
					ret(c1,c2,c3)= (*this)(c1,index,c2,c3);
					break;
				case 3:
					ret(c1,c2,c3)= (*this)(c1,c2,index,c3);
					break;
				case 4:
					ret(c1,c2,c3)= (*this)(c1,c2,c3,index);
					break;
				}
			}
	return ret;
}

/*
istream& operator >>(istream &is, Matrix4D &obj) {
	delete obj.data;
	obj.lines = count(istreambuf_iterator<char> (is), istreambuf_iterator<char> (), '\n') + 1;
	is.seekg(0);
	obj.columns = count(istreambuf_iterator<char> (is), istreambuf_iterator<char> (), ',') / obj.lines + 1;
	obj.data = new double[obj.lines*obj.columns];
	cout << "Reading Matrix2D, found " << obj.lines << " lines, " << obj.columns << " columns" << endl;
	is.seekg(0);
	for (unsigned l = 0; l < obj.lines; l++)
			for (unsigned c = 0; c < obj.columns; c++){
				is>>obj(l, c);
				is.ignore(1);
			}
	return is;
}*/

ostream& operator <<(ostream &os, const Matrix4D &obj) {
	os.precision(15);
	os.setf(ios::scientific, ios::floatfield);
	os.setf(ios::showpos);
	for (unsigned i = 0; i < size(); i++) {
		os << data[i] << ",";
	}
	return os;
}
