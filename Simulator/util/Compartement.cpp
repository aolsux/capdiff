/*
 * Compartement.cpp
 *
 *  Created on: 21.02.2011
 *      Author: Martin Rueckl
 */

#include "Compartement.h"

Compartement::Compartement(const Rectangle &_rect, const vector<CollisionObject const *> &objects):
	Rectangle	(_rect),
	ul_sub		(NULL),
	ur_sub		(NULL),
	ll_sub		(NULL),
	lr_sub		(NULL),
	collision	(NULL)
{
	if(objects.size() == 0){
		return; //no collision, this may happen when empty environment is initialized. this then is root of tree.
	} else if (objects.size() == 1) {//only colliding with one object -> abort recursion
		collision = objects[0];
	} else {
		vector<CollisionObject const *> ll_collisions = getCollisions(getLowerLeftQuadrant(), objects);
		vector<CollisionObject const *> ul_collisions = getCollisions(getUpperLeftQuadrant(), objects);
		vector<CollisionObject const *> lr_collisions = getCollisions(getLowerRightQuadrant(), objects);
		vector<CollisionObject const *> ur_collisions = getCollisions(getUpperRightQuadrant(), objects);

		if (ll_collisions.size()>0) 	ll_sub = new Compartement(getLowerLeftQuadrant(), ll_collisions);
		if (ul_collisions.size()>0) 	ul_sub = new Compartement(getUpperLeftQuadrant(), ul_collisions);
		if (lr_collisions.size()>0) 	lr_sub = new Compartement(getLowerRightQuadrant(), lr_collisions);
		if (ur_collisions.size()>0) 	ur_sub = new Compartement(getUpperRightQuadrant(), ur_collisions);
	}
}

Compartement::~Compartement(){

}

vector<CollisionObject const *> Compartement::getCollisions(const Rectangle &rect, const vector<CollisionObject const *> &objects) const
{
	vector<CollisionObject const *> collisions;
	for (unsigned u = 0; u < objects.size(); u++) {
		CollisionObject const &obj = *(objects[u]);
		bool left   = obj.isBoundaryCrossed(rect.getLowerLeftCorner(),rect.getUpperLeftCorner());
		bool top    = obj.isBoundaryCrossed(rect.getUpperLeftCorner(),rect.getUpperRightCorner());
		bool right  = obj.isBoundaryCrossed(rect.getUpperRightCorner(),rect.getLowerRightCorner());
		bool bottom = obj.isBoundaryCrossed(rect.getLowerRightCorner(),rect.getLowerLeftCorner());
		bool inside = rect.isInside(obj.getCenter());
		if (left || top || right || bottom || inside)	collisions.push_back(&obj);
	}
	return collisions;
}

CollisionObject const *Compartement::getCollision(const Point &location) const
{
	if (collision != NULL) return (collision->isInside(location) ? collision : NULL);
	Point diff = location - getLowerLeftCorner();
	bool left	= 			  0	<=	diff.x			&& 		diff.x <=	width_half;
	bool right	= 	 width_half <	diff.x			&& 		diff.x <	width;
	bool upper	= 	height_half	<	diff.y			&& 		diff.y <	height;
	bool lower	= 			  0	<=	diff.y			&& 		diff.y <=	height_half;
	if (left 	&& upper && ul_sub != NULL) return ul_sub->getCollision(location);
	if (left 	&& lower && ll_sub != NULL) return ll_sub->getCollision(location);
	if (right 	&& upper && ur_sub != NULL) return ur_sub->getCollision(location);
	if (right 	&& lower && lr_sub != NULL) return lr_sub->getCollision(location);
	return NULL;
}

ostream& operator <<(ostream &os, const Compartement &obj) {
	os << obj.ll.x << "," << obj.ll.y << "," << obj.ur.x << "," << obj.ur.y <<endl;
	if(obj.ul_sub != NULL)os << *obj.ul_sub;
	if(obj.ur_sub != NULL)os << *obj.ur_sub;
	if(obj.ll_sub != NULL)os << *obj.ll_sub;
	if(obj.lr_sub != NULL)os << *obj.lr_sub;
	return os;
}

