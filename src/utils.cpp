#include <iostream>
#include <cmath>
#include "Satellite.h"
//baselining cartesian coordinates in ECI frame

std::array<double,3> calculate_orbital_acceleration(const std::array<double,3> input_r_vec)
{
    //orbital acceleration = -G m_Earth/distance^3 * r_vec (just based on rearranging F=ma with a the acceleration due to gravitational attraction between satellite and Earth https://en.wikipedia.org/wiki/Newton%27s_law_of_universal_gravitation)
    //going to be assuming Earth's position doesn't change for now
    //also assuming Earth is spherical, can loosen this assumption in the future
    std::array<double,3> acceleration_vec=input_r_vec; //shouldn't need to explicitly call copy function because input_r_vec is passed by value not ref
    


    const double distance=sqrt(input_r_vec.at(0)*input_r_vec.at(0) + input_r_vec.at(1)*input_r_vec.at(1) + input_r_vec.at(2)*input_r_vec.at(2));
    const double overall_factor= -G*mass_Earth/pow(distance,3);



    for (size_t ind=0;ind<input_r_vec.size();ind++){
        acceleration_vec.at(ind)*=overall_factor;
    }

    return acceleration_vec;
}

std::array<double,6> RK4_deriv_function_orbit_position_and_velocity(std::array<double,6> input_position_and_velocity){
    std::array<double,6> derivative_of_input_y={};
    std::array<double,3> position_array={};

    for (size_t ind=0;ind<3;ind++){
        derivative_of_input_y.at(ind)=input_position_and_velocity.at(ind+3);
        position_array.at(ind)=input_position_and_velocity.at(ind);
    }
    
    std::array<double,3> calculated_orbital_acceleration=calculate_orbital_acceleration(position_array);

    for (size_t ind=3;ind<6;ind++){
        derivative_of_input_y.at(ind)=calculated_orbital_acceleration.at(ind-3);
    }

    return derivative_of_input_y;
}