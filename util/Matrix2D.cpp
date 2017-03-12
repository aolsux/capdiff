/*
 * Matrix2D.cpp
 *
 *  Created on: 24.11.2010
 *      Author: Martin Rueckl
 */
#include "Matrix2D.h"

#include <iostream>
#include <sstream>
#include <stdio.h>
#include <fstream>
#include <algorithm>
#include <limits>


Matrix2D::Matrix2D(unsigned _lines, unsigned _columns):
	lines(_lines),
	columns(_columns)
{
	data = new double[columns * lines];
	for (unsigned i = 0; i < columns * lines; i++)
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
Matrix2D::Matrix2D(const Matrix2D &ref) :
	lines(ref.lines),
	columns(ref.columns)
{
	data = new double[columns * lines];
	for (unsigned i = 0; i < columns * lines; i++)
		data[i] = ref.data[i];
}

/*
 * access matrix data column and line
 */
double& Matrix2D::operator()(unsigned l, unsigned c)
{
	if (l >= lines || c >= columns) {
		cerr << "c=" << c << " l=" << l << endl;
		throw MyException("Matrix subscript out of bounds");
	}
	return data[columns * l + c];
}

double Matrix2D::operator()(unsigned l, unsigned c) const
{
	if (l >= lines || c >= columns) {
		cerr << "c=" << c << " l=" << l << endl;
		throw MyException("Matrix subscript out of bounds");
	}
	return data[columns * l + c];
}

Matrix2D &Matrix2D::operator*=(const double &d)
{
	for (unsigned i = 0; i < columns * lines; i++)
		data[i] *= d;
	return *this;
}

Matrix2D &Matrix2D::operator+=(const Matrix2D &ref){
	if(ref.lines!=lines||ref.columns!=columns)
		throw MyException("Cannot sum over matrices with different size");
	for (unsigned i = 0; i < columns * lines; i++)
			data[i] += ref.data[i];
	return *this;
}

Matrix2D Matrix2D::operator*(const double &d) const {
	Matrix2D m(*this);
	m*=d;
	return m;
}

Matrix2D::~Matrix2D()
{
	delete data;
}

double Matrix2D::sum() const
{
	double ans = 0;
	for (unsigned i = 0; i < columns * lines; i++)
		ans += data[i];
	return ans;
}

double Matrix2D::getMaximum() const
{
	double ans = numeric_limits<double>::min();
	for (unsigned i = 0; i < columns * lines; i++)
		ans = max(data[i], ans);
	return ans;
}
double Matrix2D::getMinimum() const
{
	double ans =  numeric_limits<double>::max();
	for (unsigned i = 0; i < columns * lines; i++)
		ans = min(data[i], ans);
	return ans;
}



Matrix2D Matrix2D::transpose() const
{
	Matrix2D m(columns,lines);
	for (unsigned l = 0; l <  lines; l++)
		for (unsigned c = 0; c <  columns; c++)
			m(c,l)=(*this)(l,c);
	return m;
}

/*
 * Assigns given reference to *this.
 */
Matrix2D &Matrix2D::operator=(const Matrix2D &ref)
{
	//check for self assignment
	if (this == &ref) return *this;

	delete data;
	data = new double[ref.columns * ref.lines];

	lines	= ref.lines;
	columns = ref.columns;
	for (unsigned i = 0; i < columns * lines; i++)
		data[i] = ref.data[i];

	return *this;
}

///*
// * Resizes the Matrix to the given values.
// * If this means the matrix gets smaller, then elements are dropped and lost.
// * For increasing size, new elements are initialized with zeroes.
// */
//Matrix2D &Matrix2D::resize(unsigned _lines, unsigned _columns)
//{
//	double *_data = new double[_columns * _lines];
//	for (unsigned i = 0; i < _lines * _columns; i++)
//		_data[i] = 0;
//	for (unsigned l = 0; l < min(lines, _lines); l++)
//		for (unsigned c = 0; c < min(columns, _columns); c++)
//			_data[_columns * l + c] = (*this)(l, c);
//	delete data;
//	data = _data;
//	columns = _columns;
//	lines = _lines;
//	return *this;
//}

istream& operator >>(istream &is, Matrix2D &obj) {
	delete obj.data;
	char c = '\n';
	is.seekg(-1,ios::end);
	is.get(c);
	is.seekg(0,ios::beg);
	obj.lines = count(istreambuf_iterator<char> (is.rdbuf()), istreambuf_iterator<char> (), '\n') + (c=='\n'?0:1);
	std::cout << "found " << obj.lines << " newlines" << std::endl;
	is.seekg(0);
	unsigned commas = count(istreambuf_iterator<char> (is.rdbuf()), istreambuf_iterator<char> (), ',');
	std::cout << "found " << commas << "commas" << std::endl;
	is.seekg(0);
	unsigned spaces= count(istreambuf_iterator<char> (is), istreambuf_iterator<char> (), ' ');
	std::cout << "found " << spaces << " spaces" << std::endl;
	obj.columns = max(commas,spaces) / obj.lines +1;
	char delim = (commas>spaces?',':' ');
	obj.data = new double[obj.lines*obj.columns];
	cout << "Reading Matrix2D, found " << obj.lines << " lines, " << obj.columns << " columns" << endl;
	cout << "Using columns delimiter '"<<delim<<"'. Found "<<(c=='\n'?"NO ":"")<<" newline at end of file."<<endl;
	is.seekg(0);
	for (unsigned l = 0; l < obj.lines; l++)
			for (unsigned c = 0; c < obj.columns; c++){	
				is>>obj(l, c);				
				is.ignore(1);
			}
	return is;
}

ostream& operator <<(ostream &os, const Matrix2D &obj) {
	os.precision(15);
	os.setf(ios::scientific, ios::floatfield);
	os.setf(ios::showpos);
	for (unsigned l = 0; l < obj.lines; l++) {
		for (unsigned c = 0; c < obj.columns; c++)
			os << obj(l, c) << (c+1<obj.columns?",":"");
		os << (l+1<obj.lines?"\n":"");
	}
	return os;
}
