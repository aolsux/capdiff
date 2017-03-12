/*
 * FieldPreCalculator.cpp
 *
 *  Created on: Aug 21, 2009
 *      Author: gstu0908
 */

#include <iostream>
#include <sstream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <algorithm>
#include <limits>
#include <boost/filesystem.hpp>

#include <Matrix2D.h>
#include <MyException.h>
#include <Rectangle.h>

#include "CapillaryConfiguration.h"

CapillaryConfiguration::CapillaryConfiguration()									{}
CapillaryConfiguration::CapillaryConfiguration(Capillary2D cap)						{addCapillary(cap);}
CapillaryConfiguration::CapillaryConfiguration(const vector<Capillary2D> &cap)		{addCapillaries(cap);}
CapillaryConfiguration::CapillaryConfiguration(const CapillaryConfiguration &ref)	{capillaries=ref.capillaries;}


CapillaryConfiguration& CapillaryConfiguration::operator=(const CapillaryConfiguration &ref)	{capillaries=ref.capillaries; return *this;}
void 					CapillaryConfiguration::addCapillary(const Capillary2D &cap)			{capillaries.push_back(cap);}
void 					CapillaryConfiguration::addCapillaries(const vector<Capillary2D> &cap)
{
	for (unsigned u = 0; u < cap.size(); u++)
		capillaries.push_back(cap[u]);
}

double CapillaryConfiguration::getCapillaryArea() const
{
	double ret = 0;
	for (unsigned i = 0; i < capillaries.size(); i++)
		ret += capillaries[i].getArea();
	return ret;
}

unsigned CapillaryConfiguration::getCapillaryCount() const
{
	return 	 capillaries.size();
}

Point CapillaryConfiguration::getMaxPosition() const
{
	Point p(0, 0);
	if (capillaries.size() == 0) return p;
	p = capillaries[0].location;
	for (unsigned i = 1; i < capillaries.size(); i++) {
		p.x = max(p.x, capillaries[i].location.x);
		p.y = max(p.y, capillaries[i].location.y);
	}
	return p;
}

Rectangle CapillaryConfiguration::getSystemSize() const{
	Point 	max	= getMaxPosition();//upper right corner
	Point 	min	= getMinPosition();//lower left corner
	double 	r	= getMaximalRadius();
	return Rectangle(Point(min.x-r,min.y-r), Point(max.x+r,max.y+r));
}

double CapillaryConfiguration::getMinimalRadius() const
{
	double minradius = ((capillaries.size() == 0) ? 0 : capillaries[0].getRadius());
	for (unsigned i = 0; i < capillaries.size(); i++) {
		minradius = min(minradius, capillaries[i].getRadius());
	}
	return minradius;
}

double CapillaryConfiguration::getMaximalRadius() const
{
	double maxradius = ((capillaries.size() == 0) ? 0 : capillaries[0].getRadius());
	for (unsigned i = 0; i < capillaries.size(); i++) {
		maxradius = max(maxradius, capillaries[i].getRadius());
	}
	return maxradius;
}

Point CapillaryConfiguration::getMinPosition() const
{
	Point p(0, 0);
	if (capillaries.size() == 0) return p;
	p = capillaries[0].location;
	for (unsigned i = 1; i < capillaries.size(); i++) {
		p.x = min(p.x, capillaries[i].location.x);
		p.y = min(p.y, capillaries[i].location.y);
	}
	return p;
}

bool CapillaryConfiguration::isInside(const Point &p) const
{
	for (unsigned i = 0; i < capillaries.size(); i++) {
		if(capillaries[i].isInside(p))return true;
	}
	return false;
}

CapillaryConfiguration::~CapillaryConfiguration()
{
}



Matrix2D CapillaryConfiguration::getSupplyAreaHistogram(Rectangle rect, double dx) const
{
	Matrix2D m(capillaries.size(), 4);
	for (unsigned c = 0; c < capillaries.size(); c++) {
		m(c, 0) = capillaries[c].getLocation().x;
		m(c, 1) = capillaries[c].getLocation().y;
		m(c, 2) = capillaries[c].getRadius();
	}
	for (double x = rect.getLowerLeftCorner().x; x < rect.getUpperRightCorner().x; x+=dx) {
		for (double y = rect.getLowerLeftCorner().y; y < rect.getUpperRightCorner().y; y+=dx) {
			Point p(x, y);
			double min_distance = numeric_limits<double>::max();
			unsigned min_index = 0;
			for (unsigned c = 0; c < capillaries.size(); c++) {
				double distance = (p - capillaries[c].getLocation()).abs() / capillaries[c].getRadius();
				if (distance < min_distance) {
					min_distance = distance;
					min_index = c;
				}
			}
			m(min_index, 3)++;
		}
	}
	return m;
}

Matrix2D CapillaryConfiguration::getSupplyAreaMatrix(Rectangle rect, double dx) const
{
	Matrix2D m(rect.getWidth()/dx, rect.getHeight()/dx);
	for (unsigned l = 0; l < m.lines; l++) {
		for (unsigned c = 0; c < m.columns; c++) {
			Point p = Point(l*dx,c*dx) + rect.getLowerLeftCorner();
			double min_distance = numeric_limits<double>::max();
			unsigned min_index = 0;
			for (unsigned cap = 0; cap < capillaries.size(); cap++) {
				double distance = (p - capillaries[cap].getLocation()).abs() / capillaries[cap].getRadius();
				if (distance < min_distance) {
					min_distance = distance;
					min_index = cap;
				}
			}
			m(l, c) = min_index;
		}
	}
	return m;
}


CapillaryConfiguration CapillaryConfiguration::readFile(string filename){
	if (!boost::filesystem::exists(filename)) throw MyException("The File " + filename + " does not exist.");
	cout << "Reading Capillary-File " << filename << " ...\n";
	Matrix2D m(1,1);
	CapillaryConfiguration c;
	ifstream stream(filename.c_str());
	stream >> m;
	std::cout << m.lines << "x" << m.columns << std::endl;
	std::cout << m << std::endl;
	for (unsigned i = 0; i < m.lines; i++) {
		Capillary2D cap(Point(m(i,0),m(i,1)),m(i,2),m(i,3));
		//cout << "\t#" << i << ": x = " << cap.getLocation().x << "\ty = " << cap.getLocation().y << "\tr = " << cap.getRadius() << "\ta=" << cap.getAngle() << endl;
		c.addCapillary(cap);
	}
	return c;
}

ostream& operator <<(ostream &os, const CapillaryConfiguration &obj) {
	Matrix2D m(obj.capillaries.size(),4);
	for (unsigned i = 0; i < obj.capillaries.size(); i++) {
		m(i,0) = obj.capillaries[i].getLocation().x;
		m(i,1) = obj.capillaries[i].getLocation().y;
		m(i,2) = obj.capillaries[i].getRadius();
		m(i,3) = obj.capillaries[i].getAngle();
	}
	os << m;
	return os;
}
