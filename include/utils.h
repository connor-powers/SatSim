#include <iostream>

std::array<double,3> calculate_orbital_acceleration(const std::array<double,3> input_r_vec);

template <int T> std::array<double, T> RK4_step(std::array<double, T> y_n, double input_step_size,std::function<std::array<double,T>(const std::array<double,T> input_y_vec)> input_derivative_function);

std::array<double,6> RK4_deriv_function_orbit_position_and_velocity(std::array<double,6> input_position_and_velocity);

std::pair<std::array<double,3>,std::array<double,3>> RK4_step_orbital_position_and_velocity(std::array<double,3> input_position_array, std::array<double,3> input_velocity_array, double input_step_size);
