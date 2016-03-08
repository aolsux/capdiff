/*
 * Field2D.cpp
 *
 *  Created on: Aug 25, 2009
 *      Author: gstu0908
 */

#include "DiscreteField3D.h"
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

DiscreteField3D::DiscreteField3D() :
	Cuboid  (1,1,1),
	Matrix3D(2,2,2)
{
	resX = width  / ((double) slices  - 1.);
	resY = height / ((double) lines   - 1.);
	resZ = depth  / ((double) columns - 1.);
}

DiscreteField3D::DiscreteField3D(const DiscreteField3D &ref) :
	Cuboid  (ref),
	Matrix3D(ref)
{
	resX = ref.resX;
	resY = ref.resY;
	resZ = ref.resZ;
}

DiscreteField3D::DiscreteField3D(const Matrix3D &ref, Point3D ll, Point3D ur) :
	Cuboid	 (ll, ur),
	Matrix3D(ref)
{
	resX = width  / ((double) slices  - 1.);
	resY = height / ((double) lines   - 1.);
	resZ = depth  / ((double) columns - 1.);
}

DiscreteField3D::DiscreteField3D( Point3D llf, Point3D urb, double _resolution) :
	Cuboid	(llf,urb),
	Matrix3D(floor((urb.x-llf.x) / _resolution) + 1,
			 floor((urb.y-llf.y) / _resolution) + 1,
			 floor((urb.z-llf.z) / _resolution) + 1)
{
	resX = width  / ((double) slices  - 1.);
	resY = height / ((double) lines   - 1.);
	resZ = depth  / ((double) columns - 1.);
}

DiscreteField3D &DiscreteField3D::operator=(const DiscreteField3D &ref)
{
	/*call assignment operators of base classes*/
	Matrix3D::operator =(ref);
	Cuboid::operator =(ref);

	resX = ref.resX;
	resY = ref.resY;
	resZ = ref.resZ;
	return *this;
}

DiscreteField3D::~DiscreteField3D()
{
	//nothing to do....
}

string DiscreteField3D::getDescription()const{
	return	string("DiscreteField3D:")  +
					"\n\tlowerleft\t\t="		+util::num2str(llf.x)+"/"+util::num2str(llf.y)+"/"+util::num2str(llf.z)+
					"\n\twidth\t\t="		+util::num2str(width)+
					"\n\theight\t\t="		+util::num2str(height)+
					"\n\tdepth\t\t="		+util::num2str(depth)+
					"\n\tx_resolution\t="	+util::num2str(resX)+
					"\n\ty_resolution\t="	+util::num2str(resY)+
					"\n\tz_resolution\t="	+util::num2str(resZ)+
					"\n\tlines="			+util::num2str(lines)+
					"\n\tcolumns\t="		+util::num2str(columns);
					"\n\tslices\t="			+util::num2str(slices);
}

Point3D DiscreteField3D::indexToPoint(const unsigned &x, const unsigned &y, const unsigned &z) const
{
	return llf + Point3D(x * resX , y * resY , z * resZ);
}

double DiscreteField3D::getValue(const double &x, const double &y, const double &z) const
{
	Point3D shifted = Point3D(x,y,z) - llf;
	//index of next lattice point
	unsigned x0 = floor(shifted.x / resX);
	unsigned y0 = floor(shifted.y / resY);
	unsigned z0 = floor(shifted.z / resZ);

	//distance to next lattice point
	const double dx = shifted.x - x0 * resX;
	const double dy = shifted.y - y0 * resY;
	const double dz = shifted.z - z0 * resZ;

	//bilinear interpolation in x/y-Plane for z=z0
	double b1 = Matrix3D::operator()(x0	 ,y0  ,z0);
	double b2 = Matrix3D::operator()(x0+1,y0  ,z0) - Matrix3D::operator()(x0,y0,z0);
	double b3 = Matrix3D::operator()(x0	 ,y0+1,z0) - Matrix3D::operator()(x0,y0,z0);
	double b4 = Matrix3D::operator()(x0	 ,y0  ,z0) - Matrix3D::operator()(x0+1,y0,z0) - Matrix3D::operator()(x0, y0+1, z0) + Matrix3D::operator()(x0+1, y0+1,z0);
	double t1 = b1 + b2 * dx + b3 * dy + b4 * dx * dy;

	//bilinear interpolation in x/y-Plane for z=z0+1
	double c1 = Matrix3D::operator()(x0	 ,y0  ,z0+1);
	double c2 = Matrix3D::operator()(x0+1,y0  ,z0+1) - Matrix3D::operator()(x0,y0,z0+1);
	double c3 = Matrix3D::operator()(x0	 ,y0+1,z0+1) - Matrix3D::operator()(x0,y0,z0+1);
	double c4 = Matrix3D::operator()(x0	 ,y0  ,z0+1) - Matrix3D::operator()(x0+1,y0,z0+1) - Matrix3D::operator()(x0, y0+1, z0+1) + Matrix3D::operator()(x0+1, y0+1,z0+1);
	double t2 = c1 + c2 * dx + c3 * dy + c4 * dx * dy;

	//linear interpolate in z-direction
	return (1.-dz/resZ)*t1+(dz/resZ)*t2;
}

double DiscreteField3D::getValue(const Point3D &p) const
{
	return getValue(p.x, p.y, p.z);
}

double& DiscreteField3D::zeroOrderInterpolation(const Point3D &p)
{
	Point3D shifted = p - llf;
	//index of next lattice point
	unsigned x0 = round(shifted.x / resX);
	unsigned y0 = round(shifted.y / resY);
	unsigned z0 = round(shifted.z / resZ);
	return Matrix3D::operator()(x0, y0, z0);
}

void DiscreteField3D::transform(const Cuboid &cub){
	Cuboid::operator =(cub);
	resX = width  / ((double) slices  - 1.);
	resY = height / ((double) lines   - 1.);
	resZ = depth  / ((double) columns - 1.);
}

ostream& operator <<(ostream &os, const DiscreteField3D &obj) {
	os << (Matrix3D)obj;
	return os;
}

istream& operator >>(istream &is, DiscreteField3D &obj) {
	Matrix3D *m = (Matrix3D*)(&obj);
	is >> *m;
	return is;
}
