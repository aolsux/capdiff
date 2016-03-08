/*
 * Field2D.h
 *
 *  Created on: 15.10.2010
 *      Author: Martin Rueckl
 */

#ifndef FIELD2D_H_
#define FIELD2D_H_

#include <Point.h>
#include "Describable.h"

/*
 * This is only a helper interface class to assure that all Field types implement getValue() methods.
 */
class Field2D: public Describable {
public:
	virtual double getValue(const double &x, const double &y) const=0;
	virtual double getValue(const Point &p) const =0;
};

#endif /* FIELD2D_H_ */
