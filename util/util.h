/*
 * util.h
 *
 *  Created on: 17.11.2010
 *      Author: Martin Rueckl
 */

#ifndef UTIL_H_
#define UTIL_H_
#include <iomanip>
#include <sstream>

namespace util {
using namespace std;

template<typename T>
static string num2str(const T& num, unsigned prec = 8) {
	stringstream ss;
	ios_base::fmtflags ff = ss.flags();
	ff |= ios_base::floatfield;
	ff |= ios_base::fixed;
	ss.precision(prec);
	ss.flags(ff);
	ss << num;
	return ss.str();
}

}

#endif /* UTIL_H_ */
