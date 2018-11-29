#include <iostream>
#include "SolarSystem.h"

int main() {
    std::cout << "Hello, World!" << std::endl;

    double Msun = 1;
    double Mplanet = 1e-5;
    SolarSystem twobody(Msun, Mplanet);
    twobody.runsystem();
    return 0;
}