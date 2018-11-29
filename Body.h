//
// Created by Brian Dolle on 10/16/18.
//

#ifndef TWOBODYTEST_BODY_H
#define TWOBODYTEST_BODY_H


#include "Coords.h"
class Body {

public:
    double mass_;
    Coords bodycoords_;
    Body(double mass, double x, double y, double z, double vx, double vy, double vz);
    ~Body();
};



#endif //TWOBODYTEST_BODY_H
