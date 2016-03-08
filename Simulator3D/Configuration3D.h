/*
 * CapillaryConfiguration.h
 *
 *  Created on: Aug 21, 2009
 *      Author: gstu0908
 */

#ifndef CONFIGURATION3D_H_
#define CONFIGURATION3D_H_

#include <vector>
#include <iostream>
#include "Sphere3D.h"
#include <Rectangle.h>

using namespace std;

template<class Element_Type>
class Configuration3D
{
		template<class E>
		friend ostream& operator <<(ostream &os, const Configuration3D<E> &obj);
		template<class E>
		friend istream& operator >>(istream &os, const Configuration3D<E> &obj);

	public:
		Configuration3D()											{};
		Configuration3D(const Element_Type &elem)					{add(elem);};
		Configuration3D(const Configuration3D<Element_Type> &ref)	{elements=ref.elements;}
		Configuration3D(const vector<Element_Type> &elems)			{add(elems);}

		Configuration3D& operator=(const Configuration3D &ref)		{elements=ref.elements;return *this;}

		virtual ~Configuration3D();

		void add(const Element_Type &elem)							{elements.push_back(elem);}
		void add(const vector<Element_Type> &elems);

		/* Returns the total area of the contained elements.	 */
		double getTotalArea() const;

		/* Returns the number of elements contained. */
		unsigned getElementCount() const							{return elements.size();}

		/* Returns system size parameters. */
		Point3D getMinPosition() const;
		Point3D getMaxPosition() const;
		double getMinimalRadius() const;
		double getMaximalRadius() const;

		/*Simply calls isInside for all elements.*/
		bool isInside(const Point3D &p) const;

		Cuboid getCircumscribedCuboid() const;

	public:
		vector<Element_Type> elements;
};

#include "Configuration3D.cpp"

typedef Configuration3D<Sphere3D> SphereConfiguration;

#endif /* COLLISIONOBJECTCOLLECTION_H_ */
