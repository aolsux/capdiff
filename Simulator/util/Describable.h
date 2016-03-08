/*
 * Describable.h
 *
 *  Created on: 18.10.2010
 *      Author: Martin Rueckl
 */

#ifndef DESCRIBABLE_H_
#define DESCRIBABLE_H_

#include <string>

using namespace std;

class Describable {
public:
	/*
	 * allows to set a description. descriptions will we added to information files.
	 */
	void setDescription(string _description) {
		description = _description;
	}

	string getDescription() const {
		return description;
	}

protected:
	string description;
};

#endif /* DESCRIBABLE_H_ */
