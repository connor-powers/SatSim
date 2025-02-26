#include <iostream>
#include <functional>

std::array<double,3> calculate_orbital_acceleration(const std::array<double,3> input_r_vec);


std::array<double,6> RK4_deriv_function_orbit_position_and_velocity(std::array<double,6> input_position_and_velocity);

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