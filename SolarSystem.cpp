//
// Created by Brian Dolle on 10/16/18.
//

#include "SolarSystem.h"
#include "Coords.h"
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

double G = 6.67e-11;

SolarSystem::SolarSystem(double M1, double M2) {
/*
    double x1 = 0;
    double y1 = 0;
    double z1 = 0;
    double vx1 = 0;
    double vy1 = 0;
    double vz1 = 0;

//    Body sun_(M1, x1, y1, z1, vx1, vy1, vz1);
//    sun_ = Sun;
    double x2 = 10;
    double y2 = 0;
    double z2 = 0;
    double vx2 = 0;
    double vy2 = 5;  // ultimately need to write code to make this a circular /near circular orbit  - function calls here
    double vz2 = 0;
//    Body planet_(M2, x2, y2, z2, vx2, vy2, vz2);
 //   planet_ = Planet;
 */


    sun_.mass_ = 2e+30;
    planet_.mass_ = 6e+24;
    sun_.x_ = 0.0;
    sun_.y_ = 0.0;
    sun_.z_ = 0.0;
    sun_.vx_ = 0.0;
    sun_.vy_ = 0.0;
    sun_.vz_ = 0.0;

    planet_.x_ = 0.0;
    planet_.y_ = 146e+9;
    planet_.z_ = 0.0;
    planet_.vx_ = 30300;
    planet_.vy_ = 0.0;
    planet_.vz_ = 0.0;
}

SolarSystem::~SolarSystem() {

}

void SolarSystem::runsystem() {
    t_ = 0.0;
    h_ = 3600.0; // depends on scale

    double k1Sx, k1Sy, k1Sz, k2Sx, k2Sy, k2Sz, k3Sx, k3Sy, k3Sz, k4Sx, k4Sy, k4Sz;
    double k1Px, k1Py, k1Pz, k2Px, k2Py, k2Pz, k3Px, k3Py, k3Pz, k4Px, k4Py, k4Pz;
    double k1Svx, k1Svy, k1Svz, k2Svx, k2Svy, k2Svz, k3Svx, k3Svy, k3Svz, k4Svx, k4Svy, k4Svz;
    double k1Pvx, k1Pvy, k1Pvz, k2Pvx, k2Pvy, k2Pvz, k3Pvx, k3Pvy, k3Pvz, k4Pvx, k4Pvy, k4Pvz;


    double sun_tempx, sun_tempy, sun_tempz, planet_tempx, planet_tempy, planet_tempz;


    for(int i = 0; i < 10000; i++){

        // Runge Kutta Code here for all 8 ODES

        k1Sx  = h_ * sun_.vx_;
        k1Sy  = h_ * sun_.vy_;
        k1Sz  = h_ * sun_.vz_;
        k1Px  = h_ * planet_.vx_;
        k1Py  = h_ * planet_.vy_;
        k1Pz  = h_ * planet_.vz_;



        k1Svx = h_ * SunVxDot(sun_.x_, sun_.y_, sun_.z_,
                              planet_.x_, planet_.y_, planet_.z_, sun_.mass_, planet_.mass_);
        k1Svy = h_ * SunVyDot(sun_.x_, sun_.y_, sun_.z_,
                              planet_.x_, planet_.y_, planet_.z_, sun_.mass_, planet_.mass_);
        k1Svz = h_ * SunVzDot(sun_.x_, sun_.y_, sun_.z_,
                              planet_.x_, planet_.y_, planet_.z_, sun_.mass_, planet_.mass_);
        k1Pvx = h_ * PlanetVxDot(sun_.x_, sun_.y_, sun_.z_,
                                 planet_.x_, planet_.y_, planet_.z_, sun_.mass_, planet_.mass_);
        k1Pvy = h_ * PlanetVyDot(sun_.x_, sun_.y_, sun_.z_,
                                 planet_.x_, planet_.y_, planet_.z_, sun_.mass_, planet_.mass_);
        k1Pvz = h_ * PlanetVzDot(sun_.x_, sun_.y_, sun_.z_,
                                 planet_.x_, planet_.y_, planet_.z_, sun_.mass_, planet_.mass_);



        sun_tempx = sun_.x_ + k1Sx/2.0;
        sun_tempy = sun_.y_ + k1Sy/2.0;
        sun_tempz = sun_.z_ + k1Sz/2.0;
        planet_tempx = planet_.x_ + k1Px/2.0;
        planet_tempy = planet_.y_ + k1Py/2.0;
        planet_tempz = planet_.z_ + k1Pz/2.0;

        k2Sx  = h_ * (sun_.vx_ + k1Svx/2.0);
        k2Sy  = h_ * (sun_.vy_ + k1Svy/2.0);
        k2Sz  = h_ * (sun_.vz_ + k1Svz/2.0);
        k2Px  = h_ * (planet_.vx_ + k1Pvx/2.0);
        k2Py  = h_ * (planet_.vy_ + k1Pvy/2.0);
        k2Pz  = h_ * (planet_.vz_ + k1Pvz/2.0);

        k2Svx = h_ * SunVxDot(sun_tempx, sun_tempy, sun_tempz, planet_tempx, planet_tempy, planet_tempz, sun_.mass_, planet_.mass_);
        k2Svy = h_ * SunVyDot(sun_tempx, sun_tempy, sun_tempz, planet_tempx, planet_tempy, planet_tempz, sun_.mass_, planet_.mass_);
        k2Svz = h_ * SunVzDot(sun_tempx, sun_tempy, sun_tempz, planet_tempx, planet_tempy, planet_tempz, sun_.mass_, planet_.mass_);
        k2Pvx = h_ * PlanetVxDot(sun_tempx, sun_tempy, sun_tempz, planet_tempx, planet_tempy, planet_tempz, sun_.mass_, planet_.mass_);
        k2Pvy = h_ * PlanetVyDot(sun_tempx, sun_tempy, sun_tempz, planet_tempx, planet_tempy, planet_tempz, sun_.mass_, planet_.mass_);
        k2Pvz = h_ * PlanetVzDot(sun_tempx, sun_tempy, sun_tempz, planet_tempx, planet_tempy, planet_tempz, sun_.mass_, planet_.mass_);

        sun_tempx = sun_.x_ + k2Sx/2.0;
        sun_tempy = sun_.y_ + k2Sy/2.0;
        sun_tempz = sun_.z_ + k2Sz/2.0;
        planet_tempx = planet_.x_ + k2Px/2.0;
        planet_tempy = planet_.y_ + k2Py/2.0;
        planet_tempz = planet_.z_ + k2Pz/2.0;

        k3Sx  = h_ * (sun_.vx_ + k2Svx/2.0);
        k3Sy  = h_ * (sun_.vy_ + k2Svy/2.0);
        k3Sz  = h_ * (sun_.vz_ + k2Svz/2.0);
        k3Px  = h_ * (planet_.vx_ + k2Pvx/2.0);
        k3Py  = h_ * (planet_.vy_ + k2Pvy/2.0);
        k3Pz  = h_ * (planet_.vz_ + k2Pvz/2.0);

        k3Svx = h_ * SunVxDot(sun_tempx, sun_tempy, sun_tempz, planet_tempx, planet_tempy, planet_tempz, sun_.mass_, planet_.mass_);
        k3Svy = h_ * SunVyDot(sun_tempx, sun_tempy, sun_tempz, planet_tempx, planet_tempy, planet_tempz, sun_.mass_, planet_.mass_);
        k3Svz = h_ * SunVzDot(sun_tempx, sun_tempy, sun_tempz, planet_tempx, planet_tempy, planet_tempz, sun_.mass_, planet_.mass_);
        k3Pvx = h_ * PlanetVxDot(sun_tempx, sun_tempy, sun_tempz, planet_tempx, planet_tempy, planet_tempz, sun_.mass_, planet_.mass_);
        k3Pvy = h_ * PlanetVyDot(sun_tempx, sun_tempy, sun_tempz, planet_tempx, planet_tempy, planet_tempz, sun_.mass_, planet_.mass_);
        k3Pvz = h_ * PlanetVzDot(sun_tempx, sun_tempy, sun_tempz, planet_tempx, planet_tempy, planet_tempz, sun_.mass_, planet_.mass_);

        sun_tempx = sun_.x_ + k3Sx;
        sun_tempy = sun_.y_ + k3Sy;
        sun_tempz = sun_.z_ + k3Sz;
        planet_tempx = planet_.x_ + k3Px;
        planet_tempy = planet_.y_ + k3Py;
        planet_tempz = planet_.z_ + k3Pz;

        k4Sx  = h_ * (sun_.vx_ + k3Svx);
        k4Sy  = h_ * (sun_.vy_ + k3Svy);
        k4Sz  = h_ * (sun_.vz_ + k3Svz);
        k4Px  = h_ * (planet_.vx_ + k3Pvx);
        k4Py  = h_ * (planet_.vy_ + k3Pvy);
        k4Pz  = h_ * (planet_.vz_ + k3Pvz);

        k4Svx = h_ * SunVxDot(sun_tempx, sun_tempy, sun_tempz, planet_tempx, planet_tempy, planet_tempz, sun_.mass_, planet_.mass_);
        k4Svy = h_ * SunVyDot(sun_tempx, sun_tempy, sun_tempz, planet_tempx, planet_tempy, planet_tempz, sun_.mass_, planet_.mass_);
        k4Svz = h_ * SunVzDot(sun_tempx, sun_tempy, sun_tempz, planet_tempx, planet_tempy, planet_tempz, sun_.mass_, planet_.mass_);
        k4Pvx = h_ * PlanetVxDot(sun_tempx, sun_tempy, sun_tempz, planet_tempx, planet_tempy, planet_tempz, sun_.mass_, planet_.mass_);
        k4Pvy = h_ * PlanetVyDot(sun_tempx, sun_tempy, sun_tempz, planet_tempx, planet_tempy, planet_tempz, sun_.mass_, planet_.mass_);
        k4Pvz = h_ * PlanetVzDot(sun_tempx, sun_tempy, sun_tempz, planet_tempx, planet_tempy, planet_tempz, sun_.mass_, planet_.mass_);

        sun_.x_ = sun_.x_ + (k1Sx + 2.0 * k2Sx + 2.0 * k3Sx + k4Sx)/6.0;
        sun_.y_ = sun_.y_ + (k1Sy + 2.0 * k2Sy + 2.0 * k3Sy + k4Sy)/6.0;
        sun_.z_ = sun_.z_ + (k1Sz + 2.0 * k2Sz + 2.0 * k3Sz + k4Sz)/6.0;
        planet_.x_ = planet_.x_ + (k1Px + 2.0 * k2Px + 2.0 * k3Px + k4Px)/6.0;
        planet_.y_ = planet_.y_ + (k1Py + 2.0 * k2Py + 2.0 * k3Py + k4Py)/6.0;
        planet_.z_ = planet_.z_ + (k1Pz + 2.0 * k2Pz + 2.0 * k3Pz + k4Pz)/6.0;

        sun_.vx_ = sun_.vx_ + (k1Svx + 2.0 * k2Svx + 2.0 * k3Svx + k4Svx)/6.0;
        sun_.vy_ = sun_.vy_ + (k1Svy + 2.0 * k2Svy + 2.0 * k3Svy + k4Svy)/6.0;
        sun_.vz_ = sun_.vz_ + (k1Svz + 2.0 * k2Svz + 2.0 * k3Svz + k4Svz)/6.0;
        planet_.vx_ = planet_.vx_ + (k1Pvx + 2.0 * k2Pvx + 2.0 * k3Pvx + k4Pvx)/6.0;
        planet_.vy_ = planet_.vy_ + (k1Pvy + 2.0 * k2Pvy + 2.0 * k3Pvy + k4Pvy)/6.0;
        planet_.vz_ = planet_.vz_ + (k1Pvz + 2.0 * k2Pvz + 2.0 * k3Pvz + k4Pvz)/6.0;

        if (i % 100 == 0){
          //  cout << i << "," << planet_.x_  <<","<< planet_.y_ << ","<< planet_.z_ << ","<< sun_.x_ <<","<< sun_.y_ <<
            //     ","<< sun_.z_ <<endl;

            ofstream fout;
            fout.open("Bodytest.csv", ios::out | ios::app);
            fout<< i << "," << planet_.x_  <<","<< planet_.y_ << ","<< planet_.z_ << ","<< sun_.x_ <<","<< sun_.y_ <<
                ","<< sun_.z_ <<endl;
        }
    }
}

double SolarSystem::SunVxDot(double xs, double ys, double zs, double xp, double yp, double zp, double Ms, double Mp) {
    double r = sqrt((xs-xp)*(xs-xp) + (ys-yp)*(ys-yp) + (zs-zp)*(zs-zp));
    return - G * Mp * (xs-xp)/pow(r,3);
}

double SolarSystem::PlanetVxDot(double xs, double ys, double zs, double xp, double yp, double zp, double Ms, double Mp) {
    double r = sqrt((xs-xp)*(xs-xp) + (ys-yp)*(ys-yp) + (zs-zp)*(zs-zp));
    return + G * Ms * (xs-xp)/r/r/r;
}

double SolarSystem::SunVyDot(double xs, double ys, double zs, double xp, double yp, double zp, double Ms, double Mp) {
    double r = sqrt((xs-xp)*(xs-xp) + (ys-yp)*(ys-yp) + (zs-zp)*(zs-zp));
    return - G * Mp * (ys-yp)/r/r/r;
}

double SolarSystem::PlanetVyDot(double xs, double ys, double zs, double xp, double yp, double zp, double Ms, double Mp) {
    double r = sqrt((xs-xp)*(xs-xp) + (ys-yp)*(ys-yp) + (zs-zp)*(zs-zp));
    return + G * Ms * (ys-yp)/r/r/r;
}

double SolarSystem::SunVzDot(double xs, double ys, double zs, double xp, double yp, double zp, double Ms, double Mp) {
    double r = sqrt((xs-xp)*(xs-xp) + (ys-yp)*(ys-yp) + (zs-zp)*(zs-zp));
    return - G * Mp * (zs-zp)/r/r/r;
}

double SolarSystem::PlanetVzDot(double xs, double ys, double zs, double xp, double yp, double zp, double Ms, double Mp) {
    double r = sqrt((xs-xp)*(xs-xp) + (ys-yp)*(ys-yp) + (zs-zp)*(zs-zp));
    return + G * Ms * (zs-zp)/r/r/r;
}