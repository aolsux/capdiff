//============================================================================
// Name        : voronoi.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
using namespace std;

int main() {
	cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!
	return 0;
}

template<unsigned D>
class Site {
	double location[D];

	BiSector<D> getBiSector(const Site<D> &otherSite) const;
};

template<unsigned D>
class Point {
	double location[D];
	/*
	 Let S be a set of n points in the plane, called sites. For p elem S, d_p: R2->R is the
	 (Euclidean) distance from a point in R2 to p, and d: R2->R is min_(p elem S)d_p. The
	 Voronoi circle at z elem R2 is the circle centered at z of radius d(z).
	 */
	VoronoiCircle<D> getVoronoiCircle() const;
};

template<unsigned D>
class VoronoiCircle {
	double center[D];
	double radius;
};

template<unsigned D>
class BiSector {
	double hangup;
	double vectors[D];
	HalfSpace<D> getUpperHalfSpace() const;
	HalfSpace<D> getLowerHalfSpace() const;
};
