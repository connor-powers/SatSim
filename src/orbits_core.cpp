#include <iostream>
#include "Satellite.h"

int main()
{
    Satellite test_sat("../input.json");
    std::array<double,3> position=test_sat.get_position();
    std::array<double,3> velocity=test_sat.get_velocity();

    std::cout << "Position: \n";
    for (double elem:position){
        std::cout << elem << "\n";
    }
    double orbital_radius=sqrt(pow(position.at(0),2) + pow(position.at(1),2) + pow(position.at(2),2));

    std::cout << "orbital radius: " << orbital_radius << "\n";
    std::cout << "Velocity: \n";

    for (double elem:velocity){
        std::cout << elem << "\n";
    }

    std::cout << "Speed:\n";
    double speed=sqrt(pow(velocity.at(0),2) + pow(velocity.at(1),2) + pow(velocity.at(2),2));
    std::cout << speed << "\n";

    return 0;
}