//
// Created by Brian Dolle on 10/16/18.
//

#ifndef TWOBODYTEST_SOLARSYSTEM_H
#define TWOBODYTEST_SOLARSYSTEM_H

#include "Coords.h"

class SolarSystem {

public:
    SolarSystem(double M1, double M2); // constructor
    ~SolarSystem();

    double SunVxDot(double xs, double ys, double zs, double xp, double yp, double zp, double Ms, double Mp);
    double SunVyDot(double xs, double ys, double zs, double xp, double yp, double zp, double Ms, double Mp);
    double SunVzDot(double xs, double ys, double zs, double xp, double yp, double zp, double Ms, double Mp);

    double PlanetVxDot(double xs, double ys, double zs, double xp, double yp, double zp, double Ms, double Mp);
    double PlanetVyDot(double xs, double ys, double zs, double xp, double yp, double zp, double Ms, double Mp);
    double PlanetVzDot(double xs, double ys, double zs, double xp, double yp, double zp, double Ms, double Mp);



    void runsystem();

    Coords sun_;
    Coords planet_;

    double t_, h_;


};

#endif //TWOBODYTEST_SOLARSYSTEM_H
