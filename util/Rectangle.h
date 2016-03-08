/*
 * Rectangle.h
 *
 *  Created on: 21.02.2011
 *      Author: Martin Rueckl
 */

#ifndef RECTANGLE_H_
#define RECTANGLE_H_

#include "Point.h"

class Rectangle
{
	public:
		Rectangle():
			ll(0,0),
			ur(1,1),
			width		(1),
			height		(1),
			width_half	(0.5),
			height_half	(0.5)
		{}

		/*Copy Constructor*/
		Rectangle(const Rectangle &ref):
			ll(ref.ll),
			ur(ref.ur),
			width		(ref.width),
			height		(ref.height),
			width_half	(ref.width_half),
			height_half	(ref.height_half)
		{}

		/*Create a rectangle centered at the origin.*/
		Rectangle(double _width, double _height):
			ll(-_width/2,-_height/2),
			ur(_width/2,_height/2),
			width		(_width),
			height		(_height),
			width_half	(_width / 2),
			height_half	(_height / 2)
		{}

		/*Create a rectangle with lower left point _ll and upper right point _ur.*/
		Rectangle(const Point &_ll, const Point &_ur):
			ll			(_ll),
			ur			(_ur),
			width		(ur.x - ll.x),
			height		(ur.y - ll.y),
			width_half	(width / 2),
			height_half	(height / 2)
		{}

		virtual ~Rectangle(){}

		bool	isInside(const Point &p)	const {return ll.x <= p.x && ll.y <= p.y && p.x < ur.x && p.y < ur.y;}
		Point	getCenter() 				const {return Point(ll.x + width_half,	ll.y + height_half);}
		Point	getLowerLeftCorner()		const {return ll;}
		Point	getUpperLeftCorner()		const {return Point(ll.x,ur.y);}
		Point	getUpperRightCorner()		const {return ur;}
		Point	getLowerRightCorner()		const {return Point(ur.x,ll.y);}
		double	getWidth()					const {return width;}
		double	getHeight()					const {return height;}
		double	getWidthHalf()				const {return width_half;}
		double	getHeightHalf()				const {return height_half;}
		double  getArea()					const {return width*height;}

		Rectangle getUpperLeftQuadrant()  const	{return Rectangle(getLowerLeftQuadrant().getUpperLeftCorner(),getUpperRightQuadrant().getUpperLeftCorner());}
		Rectangle getLowerLeftQuadrant()  const	{return Rectangle(ll,getCenter());}
		Rectangle getUpperRightQuadrant() const	{return Rectangle(getCenter(),ur);}
		Rectangle getLowerRightQuadrant() const	{return Rectangle(getLowerLeftQuadrant().getLowerRightCorner(),getUpperRightQuadrant().getLowerRightCorner());}

	protected:
		Point ll, ur;
		double width, height, width_half, height_half;
};

class Cuboid
{
	public:

		Cuboid():
			llf			(0,0,0),
			urb			(1,1,1),
			width		(1),
			height		(1),
			depth		(1),
			width_half	(0.5),
			height_half	(0.5),
			depth_half	(0.5)
		{}

		/*Copy Constructor*/
		Cuboid(const Cuboid &ref):
			llf			(ref.llf),
			urb			(ref.urb),
			width		(ref.width),
			height		(ref.height),
			depth		(ref.depth),
			width_half	(ref.width_half),
			height_half	(ref.height_half),
			depth_half	(ref.depth_half)
		{}

		/*Create a Cuboid centered at the origin.*/
		Cuboid(double _width, double _height, double _depth):
			llf(-_width/2,-_height/2,-_depth/2),
			urb( _width/2, _height/2, _depth/2),
			width		(_width),
			height		(_height),
			depth		(_depth),
			width_half	(_width / 2),
			height_half	(_height / 2),
			depth_half	(_depth / 2)
		{}

		/*Create a Cuboid with lower left front point _llf and upper right back point _urb.*/
		Cuboid(const Point3D &_llf, const Point3D &_urb):
			llf			(_llf),
			urb			(_urb),
			width		(urb.x - llf.x),
			height		(urb.y - llf.y),
			depth		(urb.z - llf.z),
			width_half	(width / 2),
			height_half	(height / 2),
			depth_half	(depth / 2)
		{}

		virtual ~Cuboid(){}

		bool	isInside(const Point3D &p)	const {return (llf.x<=p.x) && (llf.y<=p.y) && (llf.z<=p.z) && (p.x<urb.x) && (p.y<urb.y) && (p.z<urb.z);}

		Point3D	getCenter() 				const {return Point3D(llf.x + width_half, llf.y + height_half, llf.z + depth_half);}

		Point3D	getLLFCorner()				const {return llf;}
		Point3D	getULFCorner()				const {return Point3D(llf.x 	   , llf.y + height, llf.z);}
		Point3D	getLRFCorner()				const {return Point3D(llf.x + width, llf.y		   , llf.z);}
		Point3D	getURFCorner()				const {return Point3D(llf.x + width, llf.y + height, llf.z);}

		Point3D	getLLBCorner()				const {return Point3D(llf.x 	   , llf.y 		   , llf.z + depth);}
		Point3D	getULBCorner()				const {return Point3D(llf.x 	   , llf.y + height, llf.z + depth);}
		Point3D	getLRBCorner()				const {return Point3D(llf.x + width, llf.y		   , llf.z + depth);}
		Point3D	getURBCorner()				const {return urb;}

		double	getWidth()					const {return width;}
		double	getHeight()					const {return height;}
		double	getDepth()					const {return depth;}
		double	getWidthHalf()				const {return width_half;}
		double	getHeightHalf()				const {return height_half;}
		double	getDepthHalf()				const {return depth_half;}
		double  getVolume()					const {return width*height*depth;}

		Cuboid getLLFQuadrant()		const	{return Cuboid(llf,getCenter());}
		Cuboid getULFQuadrant()		const	{return Cuboid(	Point3D(llf.x				, llf.y + height_half	, llf.z				),
															Point3D(llf.x + width_half	, llf.y + height		, llf.z + depth_half));}
		Cuboid getLRFQuadrant()		const	{return Cuboid(	Point3D(llf.x + width_half	, llf.y 				, llf.z				),
															Point3D(llf.x + width		, llf.y + height_half	, llf.z + depth_half));}
		Cuboid getURFQuadrant()		const	{return Cuboid(	Point3D(llf.x +width_half   , llf.y + height_half	, llf.z				),
															Point3D(llf.x + width		, llf.y + height		, llf.z + depth_half));}

		Cuboid getLLBQuadrant()		const	{return Cuboid(	Point3D(llf.x				, llf.y 				, llf.z	+ depth_half),
															Point3D(llf.x + width_half	, llf.y + height		, llf.z + depth		));}
		Cuboid getULBQuadrant()		const	{return Cuboid(	Point3D(llf.x				, llf.y + height_half	, llf.z	+ depth_half),
															Point3D(llf.x + width_half	, llf.y + height		, llf.z + depth		));}
		Cuboid getLRBQuadrant()		const	{return Cuboid(	Point3D(llf.x + width_half  , llf.y 				, llf.z	+ depth_half),
															Point3D(llf.x + width		, llf.y + height_half	, llf.z + depth		));}
		Cuboid getURBQuadrant()		const	{return Cuboid(getCenter(),urb);}

	protected:
		Point3D llf, urb;
		double width, height, depth, width_half, height_half, depth_half;
};


#endif /* RECTANGLE_H_ */
