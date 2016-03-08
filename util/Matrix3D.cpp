/*
 * Matrix2D.cpp
 *
 *  Created on: 24.11.2010
 *      Author: Martin Rueckl
 */
#include "Matrix3D.h"

#include <iostream>
#include <stdio.h>
#include <algorithm>
#include "Matrix2D.h"

/* File Format:
	"{l0c0s0,l0c0s1}","{l0c1s0,l0c1s1}"
	"{l1c0s0,l1c0s1}","{l1c1s0,l1c1s1}"
 */
istream& operator >>(istream &is, Matrix3D &obj) {
	delete obj.data;
	obj.lines = count(istreambuf_iterator<char> (is), istreambuf_iterator<char> (), '\n') + 1;
	is.seekg(0);
	obj.columns = count(istreambuf_iterator<char> (is), istreambuf_iterator<char> (), '}') / obj.lines + 1;
	is.seekg(0);
	obj.slices = (count(istreambuf_iterator<char> (is), istreambuf_iterator<char> (), ',') / obj.lines - (obj.columns-1))/obj.columns;
	obj.data = new double[obj.lines*obj.columns*obj.slices];
	cout << "Reading Matrix3D, found " << obj.lines << " lines, " << obj.columns << " columns" << obj.slices << " slices" << endl;
	is.seekg(0);
	for (unsigned l = 0; l < obj.lines; l++)
	{
		for (unsigned c = 0; c < obj.columns; c++)
		{
			is.ignore(2);//Ignore '"{'
			for (unsigned s = 0; s < obj.slices; s++)
			{
				is >> obj(s, l, c);
				is.ignore(1);
			}
			is.ignore(2);//Ignore '}"'
			if(c+1<obj.columns)is.ignore(1);//Ignore ','
		}
		if(l+1<obj.lines)is.ignore(1);//Ignore '\n'
	}
	return is;
}
ostream& operator <<(ostream &os, const Matrix3D &obj) {
	os.precision(15);
	os.setf(ios::scientific, ios::floatfield);
	os.setf(ios::showpos);
	for (unsigned l = 0; l < obj.lines; l++) {
		for (unsigned c = 0; c < obj.columns; c++) {
			os << "\"{";
			for (unsigned s = 0; s < obj.slices; s++)
				os << obj(s, l, c) << (s+1<obj.slices?",":"");
			os << "\"}";
			os << (c+1<obj.columns?",":"");
		}
		os << (l+1<obj.lines?"\n":"");
	}
	return os;
}


Matrix3D::Matrix3D(unsigned _slices, unsigned _lines, unsigned _columns) :
	columns(_columns),
	lines(_lines),
	slices(_slices)
{
	data = new double[columns * lines * slices];
	for (unsigned i = 0; i < columns * lines * slices; i++)
		data[i] = 0;
}


Matrix3D::Matrix3D(const Matrix3D &ref) :
	columns(ref.columns),
	lines(ref.lines),
	slices(ref.slices)
{
	data = new double[columns * lines * slices];
	for (unsigned i = 0; i < columns * lines * slices; i++)
		data[i] = ref.data[i];
}


double& Matrix3D::operator()(unsigned s, unsigned l, unsigned c)
{
	if (c >= columns || l >= lines || s >= slices) {
		cerr << "s=" << s << " l=" << l << " c=" << c << endl;
		throw MyException("Matrix subscript out of bounds");
	}
	return data[columns * lines * s + columns * l + c];
}

double Matrix3D::operator()(unsigned s, unsigned l, unsigned c) const
{
	if (c >= columns || l >= lines || s >= slices) {
		cerr << "s=" << s << " l=" << l << " c=" << c << endl;
		throw MyException("Matrix subscript out of bounds");
	}
	return data[columns * lines * s + columns * l + c];
}

Matrix3D &Matrix3D::operator*=(const double &d)
{
	for (unsigned i = 0; i < columns * lines * slices; i++)
		data[i] *= d;
	return *this;
}

Matrix3D::~Matrix3D() {
	delete data;
}

double Matrix3D::sum() const {
	double ans = 0;
	for (unsigned i = 0; i < columns * lines * slices; i++)
		ans += data[i];
	return ans;
}

Matrix2D Matrix3D::getSlice(unsigned index) const
{
	Matrix2D m(lines,columns);
	for(unsigned l=0;l<lines;l++)
		for(unsigned c=0;c<columns;c++)
			m(l,c)=(*this)(index,l,c);
	return m;
}

Matrix2D Matrix3D::getLine(unsigned index) const
{
	Matrix2D m(slices,columns);
	for(unsigned s=0;s<slices;s++)
		for(unsigned c=0;c<columns;c++)
			m(s,c)=(*this)(s,index,c);
	return m;
}

Matrix2D Matrix3D::getColumn(unsigned index) const
{
	Matrix2D m(slices,lines);
	for(unsigned s=0;s<slices;s++)
		for(unsigned l=0;l<lines;l++)
			m(s,l)=(*this)(s,l,index);
	return m;
}

Matrix3D &Matrix3D::operator=(const Matrix3D &ref) {
	//check for self assignment
	if (this == &ref) return *this;

	delete data;
	data = new double[ref.columns * ref.lines * ref.slices];

	lines	= ref.lines;
	columns = ref.columns;
	slices  = ref.slices;
	for (unsigned i = 0; i < columns * lines * slices; i++)
		data[i] = ref.data[i];

	return *this;
}
