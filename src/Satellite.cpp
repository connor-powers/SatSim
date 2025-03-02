
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <Eigen/Dense>
#include "Satellite.h"
#include "utils.h"
using Eigen::Matrix3d;


//Currently baselining use of cartesian coordinates in ECI frame (J2000 specifically, if that comes up later)

//implementations of methods for Satellite class
double Satellite::calculate_orbital_period(double input_semimajor_axis){
    //https://en.wikipedia.org/wiki/Orbital_period
    double T=2*M_PI*sqrt(pow(input_semimajor_axis,3)/(G*mass_Earth));
    return T;
}

std::pair<double,double> Satellite::calculate_eccentric_anomaly(const double input_eccentricity, const double input_true_anomaly,const double input_semimajor_axis){
    //https://en.wikipedia.org/wiki/Eccentric_anomaly
    double numerator=std::sqrt(1-pow(input_eccentricity,2))*sin(input_true_anomaly);
    double denominator=input_eccentricity+cos(input_true_anomaly);
    double eccentric_anomaly=atan2(numerator,denominator);

    std::pair<double,double> output_pair;
    output_pair.first=eccentric_anomaly;
    double computed_distance_from_Earth=input_semimajor_axis*(1-input_eccentricity*cos(eccentric_anomaly));
    output_pair.second=computed_distance_from_Earth;
    return output_pair;
}

std::array<double,3> Satellite::transform_orbital_plane_coords_to_3D_ECI_cartesian(double input_x, double input_y,double input_RAAN,double input_i,double input_arg_of_periapsis){
    //https://en.wikipedia.org/wiki/Orbital_elements#Euler_angle_transformations
    
    Matrix3d mat1 {
        {cos(input_RAAN),-sin(input_RAAN),0},
        {sin(input_RAAN),cos(input_RAAN),0},
        {0,0,1}
    };
    Matrix3d mat2 {
        {1,0,0},
        {0,cos(input_i),-sin(input_i)},
        {0,sin(input_i),cos(input_i)}
    };
    Matrix3d mat3 {
        {cos(input_arg_of_periapsis),-sin(input_arg_of_periapsis),0},
        {sin(input_arg_of_periapsis),cos(input_arg_of_periapsis),0},
        {0,0,1}
    };

    Matrix3d transformation_matrix;
    transformation_matrix=mat1*mat2*mat3;

    double transformed_x_coord=transformation_matrix(0,0)*input_x + transformation_matrix(0,1)*input_y;
    double transformed_y_coord=transformation_matrix(1,0)*input_x + transformation_matrix(1,1)*input_y;
    double transformed_z_coord=transformation_matrix(2,0)*input_x + transformation_matrix(2,1)*input_y;

    std::array<double,3> ECI_cartesian_coords={transformed_x_coord,transformed_y_coord,transformed_z_coord};
    return ECI_cartesian_coords;
}



std::pair<std::array<double,3>,std::array<double,3>> Satellite::calculate_position_and_velocity_from_orbit_params(const double input_semimajor_axis,const double input_eccentricity,const double input_true_anomaly,const double input_RAAN,const double input_i,const double input_arg_of_periapsis,const double input_orbital_period){
    //Using Kepler's equation to get position from orbital parameters
    //need mean anomaly to get eccentric anomaly, from which orbital plane positions can be determined
    //https://en.wikipedia.org/wiki/Kepler%27s_equation
    std::pair<double,double> output_pair=calculate_eccentric_anomaly(input_eccentricity,input_true_anomaly,input_semimajor_axis);
    double eccentric_anomaly=output_pair.first;
    double computed_distance_from_Earth=output_pair.second;
    double x=input_semimajor_axis*(cos(eccentric_anomaly)-input_eccentricity);

    double b=input_semimajor_axis*sqrt(1-pow(input_eccentricity,2));
    double y=b*sin(eccentric_anomaly);
    //but these are x,y coordinates in the orbital plane. Need to convert to 3D cartesian coordinates in ECI

    std::array<double,3> ECI_cartesian_coords=transform_orbital_plane_coords_to_3D_ECI_cartesian(x,y,input_RAAN,input_i,input_arg_of_periapsis);

    //https://en.wikipedia.org/wiki/Orbital_speed
    double distance=sqrt(pow(ECI_cartesian_coords.at(0),2)+pow(ECI_cartesian_coords.at(1),2)+pow(ECI_cartesian_coords.at(2),2));
    if (abs(distance-computed_distance_from_Earth)>=pow(10,-8)){
        std::cout << "distance diff was " << abs(distance-computed_distance_from_Earth) << "\n";
        //these should be the same distance
    }
    double initial_speed=sqrt(G*mass_Earth*(2/distance - 1/input_semimajor_axis));
    //But what's the orbital velocity direction?

    //http://mae-nas.eng.usu.edu/MAE_5540_Web/propulsion_systems/section2/section2.3.pdf
    //need to figure out how to calculate omega, the instantaneous rate of change of the true anomaly at that value of true anomaly
    //from Kepler's second law:
    double omega=2*pow(input_semimajor_axis,2)*M_PI*sqrt(1-pow(input_eccentricity,2))/(input_orbital_period*pow(distance,2));

    double v_r_orbital_plane=distance*omega*input_eccentricity*sin(input_true_anomaly)/(1+input_eccentricity*cos(input_true_anomaly));
    double v_theta_orbital_plane=distance*omega;
    double v_x_orbital_plane=v_r_orbital_plane*cos(input_true_anomaly) - v_theta_orbital_plane*sin(input_true_anomaly);
    double v_y_orbital_plane=v_theta_orbital_plane*cos(input_true_anomaly) + v_r_orbital_plane*sin(input_true_anomaly);

    std::array<double,3> ECI_cartesian_velocities=transform_orbital_plane_coords_to_3D_ECI_cartesian(v_x_orbital_plane,v_y_orbital_plane,input_RAAN,input_i,input_arg_of_periapsis);

    std::pair<std::array<double,3>,std::array<double,3>> output_array;
    output_array.first=ECI_cartesian_coords;
    output_array.second=ECI_cartesian_velocities;
    return output_array;
}


void Satellite::evolve_RK4(double input_step_size){
    //format input position and velocity arrays into single array for RK4 step
    std::array<double,6> combined_initial_position_and_velocity_array={};
    std::pair<std::array<double,3>,std::array<double,3>> output_position_velocity_pair={};

    for (size_t ind=0;ind<3;ind++){
        combined_initial_position_and_velocity_array.at(ind) = instantaneous_position_.at(ind);
    }
    for (size_t ind=3;ind<6;ind++){
        combined_initial_position_and_velocity_array.at(ind) = instantaneous_velocity_.at(ind-3);
    }

    std::array<double,6> output_combined_initial_position_and_velocity_array= RK4_step<6>(combined_initial_position_and_velocity_array,input_step_size,RK4_deriv_function_orbit_position_and_velocity);
    
    
    for (size_t ind=0;ind<3;ind++){
        instantaneous_position_.at(ind) = output_combined_initial_position_and_velocity_array.at(ind);
        instantaneous_velocity_.at(ind) = output_combined_initial_position_and_velocity_array.at(ind+3);
    }
    t_+=input_step_size;

    return;
}