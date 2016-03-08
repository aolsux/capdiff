/*
 * Compartement.h
 *
 *  Created on: 21.02.2011
 *      Author: Martin Rueckl
 */

#ifndef COMPARTEMENT3D_H_
#define COMPARTEMENT3D_H_

#include <vector>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include <Point.h>
#include <Rectangle.h>
#include "CollisionObject3D.h"

using namespace std;

class Compartement3D : public Cuboid
{
		friend ostream& operator <<(ostream &os, const Compartement3D &obj);

	private:
		enum Plane {XPlane,YPlane,ZPlane};

	public:
		/* recursive constructor*/
		Compartement3D(const Cuboid &cub, const vector<CollisionObject3D  const *> &objects);

		virtual ~Compartement3D();

		/* Returns NULL if no Collision was found. Otherwise the object which was collided is returned. */
		CollisionObject3D  const * getCollision(const Point3D &location) const;

	private:
		Plane plane;
		Compartement3D *child1, *child2;
		CollisionObject3D  const * collision;
};

#endif /* COMPARTEMENT3D_H_ */
