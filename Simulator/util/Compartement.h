/*
 * Compartement.h
 *
 *  Created on: 21.02.2011
 *      Author: Martin Rueckl
 */

#ifndef COMPARTEMENT_H_
#define COMPARTEMENT_H_

#include <vector>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include <Point.h>
#include <Rectangle.h>
#include "CollisionObject.h"

using namespace std;

class Compartement : public Rectangle
{
		friend istream& operator >>(istream &is, Compartement &obj);
		friend ostream& operator <<(ostream &os, const Compartement &obj);

	public:
		/* Recursive Constructor */
		Compartement(const Rectangle &_rect, const vector<CollisionObject const *> &objects);

		virtual ~Compartement();

		/* Returns NULL if no Collision was found. Otherwise the object which was collided is returned. */
		CollisionObject const * getCollision(const Point &location) const;

	private:
		Compartement *ul_sub, *ur_sub, *ll_sub, *lr_sub;
		CollisionObject const * collision;

		/* Checks for collisions of the Compartment with any of the objects in the list. */
		vector<CollisionObject const *> getCollisions(const Rectangle &_rect, const vector<CollisionObject const *> &objects) const;


};

#endif /* COMPARTEMENT_H_ */
