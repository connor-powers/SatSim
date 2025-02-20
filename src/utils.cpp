#include <iostream>
#include <cmath>
//baselining cartesian coordinates in ECI frame

std::array<double,3> calculate_orbital_acceleration(const std::array<double,3> input_r_vec)
{
    //orbital acceleration = -G m_Earth/distance^3 * r_vec
    //going to be assuming Earth's position doesn't change for now
    std::array<double,3> acceleration_vec=input_r_vec; //shouldn't need to explicitly call copy function because input_r_vec is passed by value not ref
    
    const double G=6.674*pow(10,-11);
    const double mass_earth=5.9722*pow(10,24);


    const double distance=std::sqrt(input_r_vec.at(0)*input_r_vec.at(0) + input_r_vec.at(1)*input_r_vec.at(1) + input_r_vec.at(2)*input_r_vec.at(2));
    const double overall_factor= -G*mass_earth/pow(distance,3);



    for (size_t ind=0;ind<input_r_vec.size();ind++){
        acceleration_vec.at(ind)*=overall_factor;
    }

    return acceleration_vec;
}


std::array<std::array<double,3>,2> propagate_orbit_RK4(std::function<std::array<double,3>(const std::array<double,3> input_r_vec)> input_acceleration_calculation_function, const double input_initial_t,const std::array<double,3> input_initial_position_vec,const std::array<double,3> input_initial_velocity_vec,const double input_step_size)
{
    std::array<std::array<double,3>,2> position_and_vel_at_next_timestep={}; //format: {\vec{x}(t+\Delta t),\vec{v}(t+\Delta t)}

    double t_0=input_initial_t;
    double t_1=input_initial_t+input_step_size/2;
    double t_2=t_1;
    double t_3=input_initial_t+input_step_size;

    std::array<double,3> position_vec_0=input_initial_position_vec;
    std::array<double,3> velocity_vec_0=input_initial_velocity_vec;

    std::array<double,3> acceleration_vec_0=input_acceleration_calculation_function(position_vec_0);


    return position_and_vel_at_next_timestep;
}