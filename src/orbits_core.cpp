#include <iostream>
#include "Satellite.h"
#include "utils.h"

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

    std::cout << "Orbital radius: " << orbital_radius << "\n";
    std::cout << "Velocity: \n";

    for (double elem:velocity){
        std::cout << elem << "\n";
    }

    std::cout << "Speed:\n";
    double speed=test_sat.get_speed();
    std::cout << speed << "\n";

    double test_timestep=1; //s
    int num_timesteps=1;
    test_sat.evolve_RK4(test_timestep,num_timesteps);

    std::array<double,3> evolved_position=test_sat.get_position();
    std::array<double,3> evolved_velocity=test_sat.get_velocity();

    double evolved_orbital_radius=sqrt(pow(evolved_position.at(0),2) + pow(evolved_position.at(1),2) + pow(evolved_position.at(2),2));
    std::cout << "Evolved orbital radius: " << evolved_orbital_radius << "\n";
    double evolved_speed=test_sat.get_speed();
    std::cout << "Evolved orbital speed: " << evolved_speed << "\n";


    return 0;
}