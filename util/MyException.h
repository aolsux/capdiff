/*
 * MeanStepSizeException.h
 *
 *  Created on: Aug 18, 2009
 *      Author: gstu0908
 */

#ifndef MEANSTEPSIZEEXCEPTION_H_
#define MEANSTEPSIZEEXCEPTION_H_

#include <iostream>
#include <string>
using namespace std;

class MyException: public std::exception {
public:
	MyException(string _str) throw () :
		str(_str) {
	}

	virtual ~MyException() throw () {
	}

	virtual const char* what() const throw () {
		cout << str << endl;
		return str.c_str();
	}

private:
	string str;

};

#endif /* MEANSTEPSIZEEXCEPTION_H_ */

