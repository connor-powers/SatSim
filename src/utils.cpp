#include <iostream>
#include <cmath>
#include "utils.h"
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


template <int T> std::array<double, T> RK4_step(std::array<double, T> y_n, double input_step_size,std::function<std::array<double,T>(const std::array<double,T> input_y_vec)> input_derivative_function){
    //ref: https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#The_Runge%E2%80%93Kutta_method
    std::array<double,T> y_nplus1=y_n;

    //first, k=1;
    //going to be assuming derivative function does not have explicit time dependence
    std::array<double, T> k_1=input_derivative_function(y_n);

    //now k_2
    std::array<double, T> y_vec_for_k_2=y_n;
    for (size_t ind=0;ind<T;ind++){
        y_vec_for_k_2.at(ind)+=( (input_step_size/2) * k_1.at(ind));
    }
    
    std::array<double, T> k_2=input_derivative_function(y_vec_for_k_2);
    //now k_3
    std::array<double, T> y_vec_for_k_3=y_n;
    for (size_t ind=0;ind<T;ind++){
        y_vec_for_k_3.at(ind)+=( (input_step_size/2) * k_2.at(ind));
    }
    std::array<double, T> k_3=input_derivative_function(y_vec_for_k_3);

    //now k_4
    std::array<double, T> y_vec_for_k_4=y_n;
    for (size_t ind=0;ind<T;ind++){
        y_vec_for_k_4.at(ind)+=( input_step_size * k_3.at(ind));
    }
    std::array<double, T> k_4=input_derivative_function(y_vec_for_k_4);

    for (size_t ind=0;ind<T;ind++){
        y_nplus1.at(ind)+=( (input_step_size/6) * (k_1.at(ind) + 2*k_2.at(ind) + 2*k_3.at(ind) + k_4.at(ind)) );
    }
    
    return y_nplus1;
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

std::pair<std::array<double,3>,std::array<double,3>> RK4_step_orbital_position_and_velocity(std::array<double,3> input_position_array, std::array<double,3> input_velocity_array, double input_step_size){
    //format input position and velocity arrays into single array for RK4 step
    std::array<double,6> combined_initial_position_and_velocity_array={};
    std::pair<std::array<double,3>,std::array<double,3>> output_position_velocity_pair={};

    for (size_t ind=0;ind<3;ind++){
        combined_initial_position_and_velocity_array.at(ind) = input_position_array.at(ind);
    }
    for (size_t ind=3;ind<6;ind++){
        combined_initial_position_and_velocity_array.at(ind) = input_velocity_array.at(ind-3);
    }

    std::array<double,6> output_combined_initial_position_and_velocity_array= RK4_step<6>(combined_initial_position_and_velocity_array,input_step_size,RK4_deriv_function_orbit_position_and_velocity);
    for (size_t ind=0;ind<3;ind++){
        output_position_velocity_pair.first.at(ind) = output_combined_initial_position_and_velocity_array.at(ind);
        output_position_velocity_pair.second.at(ind) = output_combined_initial_position_and_velocity_array.at(ind+3);
    }

    return output_position_velocity_pair;
}
