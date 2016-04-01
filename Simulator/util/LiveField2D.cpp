//
// Created by mrueckl on 01/04/16.
//

#include "LiveField2D.h"

#include <iostream>

using namespace std;

ostream &operator<<(ostream &os, const LiveField2D &obj) {
    os << obj.capillaryConfiguration;
    return os;
}

