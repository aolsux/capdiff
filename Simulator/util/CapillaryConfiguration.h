/*
 * CapillaryConfiguration.h
 *
 *  Created on: Aug 21, 2009
 *      Author: gstu0908
 */

#ifndef CAPILLARYCONFIGURATION_H_
#define CAPILLARYCONFIGURATION_H_

#include <vector>
#include <string.h>

#include <Matrix2D.h>
#include <Rectangle.h>
#include <MyException.h>

#include "Capillary2D.h"

using namespace std;

class CapillaryConfiguration
{
		friend ostream& operator <<(ostream &os, const CapillaryConfiguration &obj);

	public:
		CapillaryConfiguration();
		CapillaryConfiguration(Capillary2D capillary);
		CapillaryConfiguration(const CapillaryConfiguration &ref);
		CapillaryConfiguration(const vector<Capillary2D> &capillaries);

		CapillaryConfiguration& operator=(const CapillaryConfiguration &ref);

		virtual ~CapillaryConfiguration();

		/* Add capillaries to the system. */
		void addCapillary(const Capillary2D &cap);
		void addCapillaries(const vector<Capillary2D> &cap);

		/* Returns the total area of the contained capillaries.	 */
		double getCapillaryArea() const;

		/* Returns the number of capillaries contained. */
		unsigned getCapillaryCount() const;

		/* Returns system size parameters. */
		Point getMinPosition() const;
		Point getMaxPosition() const;
		Rectangle getSystemSize() const;
		double getMinimalRadius() const;
		double getMaximalRadius() const;

		/*Simply calls isInside for all capillaries.*/
		bool isInside(const Point &p) const;

		/*Returns a Matrix containing a histogram for the size of the supply areas.*/
		Matrix2D getSupplyAreaHistogram(Rectangle rect, double dx) const;
		Matrix2D getSupplyAreaMatrix(Rectangle rect, double dx) const;

		static CapillaryConfiguration readFile(string filename);

		vector<Capillary2D> capillaries;
};

#endif /* CAPILLARYCONFIGURATION_H_ */
