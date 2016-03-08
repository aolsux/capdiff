/*
 * Field2D.cpp
 *
 *  Created on: Aug 25, 2009
 *      Author: gstu0908
 */

#include "PrecalculatedField2D.h"
#include <MyException.h>

#include <iostream>
#include <sstream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <algorithm>
#include <Matrix2D.h>
#include <boost/filesystem.hpp>

PrecalculatedField2D::PrecalculatedField2D() :
	Rectangle(1,1), Matrix2D(2,2)
{
	resX = width / ((double) lines - 1.);
	resY = height / ((double) columns - 1.);
	updateDescription();
}

PrecalculatedField2D::PrecalculatedField2D(const PrecalculatedField2D &ref) :
	Rectangle(ref), Matrix2D(ref)
{
	resX = ref.resX;
	resY = ref.resY;
	updateDescription();
}

PrecalculatedField2D::PrecalculatedField2D(const Matrix2D &ref, Point ll, Point ur) :
	Rectangle(ll, ur), Matrix2D(ref)
{
	resX = width / ((double) lines - 1.);
	resY = height / ((double) columns - 1.);
	updateDescription();
}

PrecalculatedField2D::PrecalculatedField2D( Point ll, Point ur, double _resolution) :
	Rectangle(ll,ur), Matrix2D(floor((ur.x-ll.x) / _resolution) + 1, floor((ur.y-ll.y) / _resolution) + 1)
{
	resX = width / ((double) lines - 1.);
	resY = height / ((double) columns - 1.);
	updateDescription();
}

PrecalculatedField2D &PrecalculatedField2D::operator=(const PrecalculatedField2D &ref)
{
	/*call assignment operators of base classes*/
	Matrix2D::operator =(ref);
	Rectangle::operator =(ref);

	resX = ref.resX;
	resY = ref.resY;
	return *this;
}

PrecalculatedField2D::~PrecalculatedField2D()
{
	//nothing to do....
}

void PrecalculatedField2D::updateDescription(){
	description=	string("Precalculated Field:")  +
					"\n\tlowerleft\t\t="		+util::num2str(ll.x)+"/"+util::num2str(ll.y)+
					"\n\twidth\t\t="		+util::num2str(width)+
					"\n\theight\t\t="		+util::num2str(height)+
					"\n\tx_resolution\t="	+util::num2str(resX)+
					"\n\ty_resolution\t="	+util::num2str(resY)+
					"\n\tlines="			+util::num2str(lines)+
					"\n\tcolumns\t="		+util::num2str(columns);
}

Point PrecalculatedField2D::transform(const unsigned &x, const unsigned &y) const
{
	return ll + Point(x * resX , y * resY);
}

double PrecalculatedField2D::getValue(const double &x, const double &y) const
{
	Point shifted = Point(x,y) - ll;
	//index of next lattice point
	unsigned x0 = floor(shifted.x / resX);
	unsigned y0 = floor(shifted.y / resY);

	//distance to next lattice point
	const double dx = shifted.x - x0 * resX;
	const double dy = shifted.y - y0 * resY;

	//linear interpolation
	double b1 = Matrix2D::operator()(x0	, y0	);
	double b2 = Matrix2D::operator()(x0 + 1, y0	) - Matrix2D::operator()(x0	, y0);
	double b3 = Matrix2D::operator()(x0	, y0 + 1) - Matrix2D::operator()(x0	, y0);
	double b4 = Matrix2D::operator()(x0	, y0	) - Matrix2D::operator()(x0 + 1, y0) - Matrix2D::operator()(x0, y0 + 1) + Matrix2D::operator()(x0 + 1, y0 + 1);
	return b1 + b2 * dx + b3 * dy + b4 * dx * dy;
}

double PrecalculatedField2D::getValue(const Point &p) const
{
	return getValue(p.x, p.y);
}

double& PrecalculatedField2D::operator()(const Point &p)
{
	Point shifted = p - ll;
	//index of next lattice point
	unsigned x0 = round(shifted.x / resX);
	unsigned y0 = round(shifted.y / resY);
	return Matrix2D::operator()(x0, y0);
}

ostream& operator <<(ostream &os, const PrecalculatedField2D &obj) {
	//os << static_cast<Matrix2D>(obj);
	os << (Matrix2D)obj;
	return os;
}

PrecalculatedField2D PrecalculatedField2D::readFile(string filename, Point ll, Point ur){
	if (!boost::filesystem::exists(filename)) throw MyException("The File " + filename + " does not exist.");
	cout << "Reading Field2D-File " << filename << "... (";
	double t1 = (double) clock() / CLOCKS_PER_SEC;
	ifstream stream(filename.c_str());
	Rectangle r(ll,ur);
	Matrix2D m(1,1);
	stream >> m;
	cout << m.lines << " lines, ";
	cout << m.columns << "columns ";
	double t2 = (double) clock() / CLOCKS_PER_SEC;
	cout << (t2 - t1) * 1000 << "ms)" << endl;
	/* width<->lines  height<->columns */
	/* This is for backwards compatibility if field file is stored with inverted dimensions, transpose matrix.*/
	if((r.getWidth()>r.getHeight()*1.1 && m.lines < m.columns) || (r.getHeight()>r.getWidth()*1.1 && m.columns < m.lines))
	{
		cout << "Assuming inverted matrix dimensions for " << filename <<" -> Transposing matrix" << endl;
		m=m.transpose();
	}
	return PrecalculatedField2D(m, ll, ur);
}
