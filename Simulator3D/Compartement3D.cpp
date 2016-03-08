/*
 * Compartement.cpp
 *
 *  Created on: 21.02.2011
 *      Author: Martin Rueckl
 */

#include "Compartement3D.h"
#include "Rectangle.h"

Compartement3D::Compartement3D(const Cuboid &cub, const vector<CollisionObject3D const *> &objects):
	Cuboid		(cub),
	child1		(NULL),
	child2		(NULL),
	collision	(NULL)
{
	if(objects.size() == 0){
		return; //no collision, this might happen when empty environment is initialized. this then is root of tree.
	} else if (objects.size() == 1) {//only colliding with one object -> abort recursion
		collision = objects[0];
	} else {

		if		(width	>= height && width >= depth)	plane = XPlane;
		else if	(height >= depth)						plane = YPlane;
		else											plane = ZPlane;

		vector<CollisionObject3D const *> child1_collisions, child2_collisions;

		for (unsigned u = 0; u < objects.size(); u++) {
			CollisionObject3D const &obj = *(objects[u]);
			bool flg_child1		= false; //indicates whether object is located in half space 1
			bool collision		= false; //indicated whether object has collision with plane separating the half spaces
			switch (plane)
			{
				case XPlane:
					flg_child1 = obj.getCenter().x < llf.x+width_half;
					collision  = obj.isPlaneCrossed(llf+Point3D(width_half,0,0),Point3D(0,0,1),Point3D(0,1,0));
					break;
				case YPlane:
					flg_child1 = obj.getCenter().y < llf.y+height_half;
					collision  = obj.isPlaneCrossed(llf+Point3D(0,height_half,0),Point3D(1,0,0),Point3D(0,0,1));
					break;
				case ZPlane:
					flg_child1 = obj.getCenter().z < llf.z+depth_half;
					collision  = obj.isPlaneCrossed(llf+Point3D(0,0,depth_half),Point3D(1,0,0),Point3D(0,1,0));
					break;
			}
			if (collision){
				child1_collisions.push_back(&obj);
				child2_collisions.push_back(&obj);
			}else{
				if(flg_child1)	child1_collisions.push_back(&obj);
				else			child2_collisions.push_back(&obj);
			}
		}

		Cuboid c1_cub,c2_cub;
		switch (plane)
		{
			case XPlane:
					c1_cub = Cuboid(Point3D(llf.x				,llf.y				,llf.z				),
									Point3D(llf.x+width_half	,llf.y+height		,llf.z+depth		));

					c2_cub = Cuboid(Point3D(llf.x+width_half	,llf.y				,llf.z				),
									Point3D(llf.x+width			,llf.y+height		,llf.z+depth		));
				break;
			case YPlane:
					c1_cub = Cuboid(Point3D(llf.x				,llf.y				,llf.z				),
									Point3D(llf.x+width			,llf.y+height_half	,llf.z+depth		));

					c2_cub = Cuboid(Point3D(llf.x				,llf.y+height_half	,llf.z				),
									Point3D(llf.x+width			,llf.y+height	 	,llf.z+depth		));
				break;
			case ZPlane:
					c1_cub = Cuboid(Point3D(llf.x				,llf.y				,llf.z				),
									Point3D(llf.x+width			,llf.y+height		,llf.z+depth_half	));

					c2_cub = Cuboid(Point3D(llf.x				,llf.y				,llf.z+depth_half	),
									Point3D(llf.x+width			,llf.y+height		,llf.z+depth		));
				break;
		}

		if(child1_collisions.size()>0)	child1 = new Compartement3D(c1_cub, child1_collisions);
		if(child2_collisions.size()>0)	child2 = new Compartement3D(c2_cub, child2_collisions);

	}
}

Compartement3D::~Compartement3D(){

}

CollisionObject3D const *Compartement3D::getCollision(const Point3D &location) const
{
	if (collision != NULL) return (collision->isInside(location) ? collision : NULL);
	switch (plane)
	{
		case XPlane:
			if(location.x < llf.x+width_half) return (child1 == NULL ? NULL : child1->getCollision(location));
			else							  return (child2 == NULL ? NULL : child2->getCollision(location));
		case YPlane:
			if(location.y < llf.y+height_half)return (child1 == NULL ? NULL : child1->getCollision(location));
			else							  return (child2 == NULL ? NULL : child2->getCollision(location));
		case ZPlane:
			if(location.z < llf.z+depth_half) return (child1 == NULL ? NULL : child1->getCollision(location));
			else							  return (child2 == NULL ? NULL : child2->getCollision(location));

	}
	return NULL;
}

ostream& operator <<(ostream &os, const Compartement3D &obj) {
	os << obj.llf.x << "," << obj.llf.y << ","<< obj.llf.z << "," << obj.urb.x << "," << obj.urb.y <<"," << obj.urb.z <<endl;
	if(obj.child1 != NULL)os << *obj.child1;
	if(obj.child2 != NULL)os << *obj.child2;
	return os;
}

