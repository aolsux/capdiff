/*
 * FieldPreCalculator.cpp
 *
 *  Created on: Aug 21, 2009
 *      Author: gstu0908
 */

#ifndef CONFIGURATION3D_CPP_
#define CONFIGURATION3D_CPP_

#include "Configuration3D.h"
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <algorithm>
#include <limits>

#include <Matrix2D.h>
#include <MyException.h>

template<typename E>
void Configuration3D<E>::add(const vector<E> &elems)
{
	for (unsigned u = 0; u < elems.size(); u++)
		elements.push_back(elems[u]);
}

template<typename E>
double Configuration3D<E>::getTotalArea() const
{
	double ret = 0;
	for (unsigned i = 0; i < elements.size(); i++)
		ret += elements[i].getArea();
	return ret;
}

template<typename E>
Point3D Configuration3D<E>::getMaxPosition() const
{
	Point3D p = (elements.size() == 0 ? Point3D(0,0,0):elements[0].getCenter());
	for (unsigned i = 1; i < elements.size(); i++) {
		p.x = max(p.x, elements[i].getCenter().x);
		p.y = max(p.y, elements[i].getCenter().y);
		p.z = max(p.z, elements[i].getCenter().z);
	}
	return p;
}

template<typename E>
double Configuration3D<E>::getMinimalRadius() const
{
	double minradius = ((elements.size() == 0) ? 0 : elements[0].getRadius());
	for (unsigned i = 0; i < elements.size(); i++) {
		minradius = min(minradius, elements[i].getRadius());
	}
	return minradius;
}

template<typename E>
double Configuration3D<E>::getMaximalRadius() const
{
	double maxradius = ((elements.size() == 0) ? 0 : elements[0].getRadius());
	for (unsigned i = 0; i < elements.size(); i++) {
		maxradius = max(maxradius, elements[i].getRadius());
	}
	return maxradius;
}

template<typename E>
Point3D Configuration3D<E>::getMinPosition() const
{
	Point3D p = (elements.size() == 0 ? Point3D(0,0,0):elements[0].getCenter());
	for (unsigned i = 1; i < elements.size(); i++) {
		p.x = min(p.x, elements[i].getCenter().x);
		p.y = min(p.y, elements[i].getCenter().y);
		p.z = max(p.z, elements[i].getCenter().z);
	}
	return p;
}

template<typename E>
bool Configuration3D<E>::isInside(const Point3D &p) const
{
	for (unsigned i = 0; i < elements.size(); i++) {
		if(elements[i].isInside(p))return true;
	}
	return false;
}
template<typename E>
Cuboid Configuration3D<E>::getCircumscribedCuboid() const
{
	double maxRadius = getMaximalRadius();
	Point3D min 	 = getMinPosition();
	Point3D max 	 = getMaxPosition();
	Point3D shift(maxRadius,maxRadius,maxRadius);
	return Cuboid(min-shift, max+shift);
}

istream& operator >>(istream &is, Configuration3D<Sphere3D> &obj) {
	obj.elements.clear();
	Matrix2D m(1,1);
	is>>m;
	for (unsigned i = 0; i < m.lines; i++)
		obj.add(Sphere3D(Point3D(m(i,0),m(i,1),m(i,2)),m(i,3)));
	return is;
}

ostream& operator <<(ostream &os, const Configuration3D<Sphere3D> &obj) {
	Matrix2D m(obj.elements.size(),3);
	for (unsigned i = 0; i < obj.elements.size(); i++) {
		m(i,0) = obj.elements[i].getCenter().x;
		m(i,1) = obj.elements[i].getCenter().y;
		m(i,2) = obj.elements[i].getRadius();
	}
	os << m;
	return os;
}

#endif /* CONFIGURATION3D_CPP_*/
