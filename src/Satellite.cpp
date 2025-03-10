
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "Satellite.h"
#include "utils.h"
using Eigen::Matrix3d;
using Eigen::Vector3d;


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


std::array<double,3> Satellite::calculate_perifocal_position(){
    //Using approach from Fundamentals of Astrodynamics
    std::array<double,3> calculated_perifocal_position;
    double mu=G*mass_Earth;

    double p=a_*(1-eccentricity_)*(1+eccentricity_); //rearranging eq. 1-47

    double r=p/(1+eccentricity_*cos(true_anomaly_));

    double r_perifocal_P=r*cos(true_anomaly_);
    double r_perifocal_Q=r*sin(true_anomaly_);

    calculated_perifocal_position.at(0)=r_perifocal_P;
    calculated_perifocal_position.at(1)=r_perifocal_Q;
    calculated_perifocal_position.at(2)=0; //W component is 0 

    return calculated_perifocal_position;
}

std::array<double,3> Satellite::calculate_perifocal_velocity(){
    //Using approach from Fundamentals of Astrodynamics
    std::array<double,3> calculated_perifocal_velocity;
    double mu=G*mass_Earth;
    double p=a_*(1-eccentricity_)*(1+eccentricity_); //rearranging eq. 1-47

    double v_perifocal_P=sqrt(mu/p)*(-sin(true_anomaly_));
    double v_perifocal_Q=sqrt(mu/p)*(eccentricity_+cos(true_anomaly_));

    calculated_perifocal_velocity.at(0)=v_perifocal_P;
    calculated_perifocal_velocity.at(1)=v_perifocal_Q;
    calculated_perifocal_velocity.at(2)=0; //0 component in the W direction

    return calculated_perifocal_velocity;
}


std::array<double,3> Satellite::convert_perifocal_to_ECI(std::array<double,3> input_perifocal_vec){
    //method from Fundamentals of Astrodynamics
    //sounds like what they describe as their "Geocentric-Equatorial" coordinate system is ECI
    Matrix3d R_tilde_matrix;

    R_tilde_matrix(0,0)=cos(raan_)*cos(arg_of_periapsis_)-sin(raan_)*sin(arg_of_periapsis_)*cos(inclination_);

    R_tilde_matrix(0,1)=-cos(raan_)*sin(arg_of_periapsis_) - sin(raan_)*cos(arg_of_periapsis_)*cos(inclination_);

    R_tilde_matrix(0,2)=sin(raan_)*sin(inclination_);

    R_tilde_matrix(1,0)=sin(raan_)*cos(arg_of_periapsis_) + cos(raan_)*sin(arg_of_periapsis_)*cos(inclination_);

    R_tilde_matrix(1,1)=-sin(raan_)*sin(arg_of_periapsis_) + cos(raan_)*cos(arg_of_periapsis_)*cos(inclination_);

    R_tilde_matrix(1,2) = -cos(raan_)*sin(inclination_);

    R_tilde_matrix(2,0)=sin(arg_of_periapsis_)*sin(inclination_);

    R_tilde_matrix(2,1)= cos(arg_of_periapsis_)*sin(inclination_);

    R_tilde_matrix(2,2) = cos(inclination_);

    //convert to Eigen vectors just to make sure matrix-vector multiplication is done correctly
    Vector3d vector_ijk(0,0,0);
    Vector3d vector_pqw(0,0,0);

    vector_pqw(0)=input_perifocal_vec.at(0);
    vector_pqw(1)=input_perifocal_vec.at(1);
    vector_pqw(2)=input_perifocal_vec.at(2);


    vector_ijk=R_tilde_matrix*vector_pqw;
    std::array<double,3> output_vector_ijk={vector_ijk(0),vector_ijk(1),vector_ijk(2)};

    return output_vector_ijk;
}




std::array<double,3> Satellite::convert_ECI_to_perifocal(std::array<double,3> input_ECI_vec){
    //method from Fundamentals of Astrodynamics, and using https://en.wikipedia.org/wiki/Perifocal_coordinate_system
    //sounds like what they describe as their "Geocentric-Equatorial" coordinate system is ECI
    Matrix3d R_tilde_matrix;

    R_tilde_matrix(0,0)=cos(raan_)*cos(arg_of_periapsis_)-sin(raan_)*sin(arg_of_periapsis_)*cos(inclination_);

    R_tilde_matrix(0,1)=-cos(raan_)*sin(arg_of_periapsis_) - sin(raan_)*cos(arg_of_periapsis_)*cos(inclination_);

    R_tilde_matrix(0,2)=sin(raan_)*sin(inclination_);

    R_tilde_matrix(1,0)=sin(raan_)*cos(arg_of_periapsis_) + cos(raan_)*sin(arg_of_periapsis_)*cos(inclination_);

    R_tilde_matrix(1,1)=-sin(raan_)*sin(arg_of_periapsis_) + cos(raan_)*cos(arg_of_periapsis_)*cos(inclination_);

    R_tilde_matrix(1,2) = -cos(raan_)*sin(inclination_);

    R_tilde_matrix(2,0)=sin(arg_of_periapsis_)*sin(inclination_);

    R_tilde_matrix(2,1)= cos(arg_of_periapsis_)*sin(inclination_);

    R_tilde_matrix(2,2) = cos(inclination_);

    //convert to Eigen vectors just to make sure matrix-vector multiplication is done correctly
    Vector3d vector_ijk(0,0,0);
    Vector3d vector_pqw(0,0,0);
    vector_ijk(0)=input_ECI_vec.at(0);
    vector_ijk(1)=input_ECI_vec.at(1);
    vector_ijk(2)=input_ECI_vec.at(2);

    vector_pqw=(R_tilde_matrix.inverse())*vector_ijk;
    std::array<double,3> output_vector_pqw={vector_pqw(0),vector_pqw(1),vector_pqw(2)};

    return output_vector_pqw;
}












void Satellite::evolve_RK4(double input_step_size){
    //format input position and velocity arrays into single array for RK4 step
    std::array<double,6> combined_initial_position_and_velocity_array={};
    std::pair<std::array<double,3>,std::array<double,3>> output_position_velocity_pair={};

    for (size_t ind=0;ind<3;ind++){
        combined_initial_position_and_velocity_array.at(ind) = ECI_position_.at(ind);
    }
    for (size_t ind=3;ind<6;ind++){
        combined_initial_position_and_velocity_array.at(ind) = ECI_velocity_.at(ind-3);
    }

    //populate list of thrust forces at half a timestep past, for RK4 calculations
    std::vector<std::array<double,3>> list_of_LVLH_forces_at_half_timestep_past={};
    std::vector<std::array<double,3>> list_of_ECI_forces_at_half_timestep_past={};

    for (ThrustProfileLVLH thrust_profile : thrust_profile_list_){

        if (((t_+(input_step_size/2))>=thrust_profile.t_start_)&&((t_+(input_step_size/2))<=thrust_profile.t_end_)){
            list_of_LVLH_forces_at_half_timestep_past.push_back(thrust_profile.LVLH_force_vec_);
            std::array<double,3> ECI_thrust_vector=convert_LVLH_to_ECI_manual(thrust_profile.LVLH_force_vec_);
            list_of_ECI_forces_at_half_timestep_past.push_back(ECI_thrust_vector);
        }
    }

    //populate list of thrust forces at a timestep past, for RK4 calculations
    std::vector<std::array<double,3>> list_of_LVLH_forces_at_one_timestep_past={};
    std::vector<std::array<double,3>> list_of_ECI_forces_at_one_timestep_past={};

    for (ThrustProfileLVLH thrust_profile : thrust_profile_list_){

        if (((t_+input_step_size)>=thrust_profile.t_start_)&&((t_+input_step_size)<=thrust_profile.t_end_)){
            list_of_LVLH_forces_at_one_timestep_past.push_back(thrust_profile.LVLH_force_vec_);
            std::array<double,3> ECI_thrust_vector=convert_LVLH_to_ECI_manual(thrust_profile.LVLH_force_vec_);
            list_of_ECI_forces_at_one_timestep_past.push_back(ECI_thrust_vector);
        }
    }
    

    std::array<double,6> output_combined_position_and_velocity_array= RK4_step<6>(combined_initial_position_and_velocity_array,input_step_size,RK4_deriv_function_orbit_position_and_velocity,m_,list_of_ECI_forces_at_this_time_,list_of_ECI_forces_at_half_timestep_past,list_of_ECI_forces_at_one_timestep_past);
    // std::array<double,6> output_combined_angular_array= RK4_step<6>(combined_initial_angular_array,input_step_size,RK4_deriv_function_angular,I_,list_of_body_frame_torques_at_this_time_,list_of_body_frame_torques_at_half_timestep_past,list_of_body_frame_torques_at_one_timestep_past);

    
    for (size_t ind=0;ind<3;ind++){
        ECI_position_.at(ind) = output_combined_position_and_velocity_array.at(ind);
        ECI_velocity_.at(ind) = output_combined_position_and_velocity_array.at(ind+3);
        //also update the perifocal versions
        perifocal_position_=convert_ECI_to_perifocal(ECI_position_);
        perifocal_velocity_=convert_ECI_to_perifocal(ECI_velocity_);
    }
    t_+=input_step_size;

    list_of_LVLH_forces_at_this_time_=list_of_LVLH_forces_at_one_timestep_past;
    list_of_ECI_forces_at_this_time_=list_of_ECI_forces_at_one_timestep_past;

    return;
}


// std::array<double,3> Satellite::convert_LVLH_to_ECI(std::array<double,3> input_LVLH_vec){
//     Vector3d input_LVLH_vec_eigen;
//     input_LVLH_vec_eigen << input_LVLH_vec.at(0),input_LVLH_vec.at(1),input_LVLH_vec.at(2);
//     //ref: https://ntrs.nasa.gov/api/citations/20205003902/downloads/Introduction%20to%20Orbital%20Mechanics%20and%20Spacecraft%20Attitudes%20for%20Thermal%20Engineers%20CHARTS%20PDF.pdf
//     Matrix3d ref_mat;
//     ref_mat << 0,0,-1,
//                1,0,0,
//                0,-1,0;
    
//     Matrix3d omega_mat=z_rot_matrix(raan_);
//     Matrix3d inclination_mat=y_rot_matrix(inclination_);
//     Matrix3d arg_periapsis_mat=z_rot_matrix(arg_of_periapsis_);
//     Matrix3d true_anomaly_mat=z_rot_matrix(true_anomaly_);

//     Matrix3d pitch_angle_mat=y_rot_matrix(pitch_angle_);
//     Matrix3d yaw_angle_mat=z_rot_matrix(yaw_angle_);
//     Matrix3d roll_angle_mat=x_rot_matrix(roll_angle_);


//     Vector3d ECI_vec_eigen=omega_mat*inclination_mat*arg_periapsis_mat*true_anomaly_mat*ref_mat*pitch_angle_mat*yaw_angle_mat*roll_angle_mat*input_LVLH_vec_eigen;

//     std::array<double,3> ECI_vec_output={ECI_vec_eigen(0),ECI_vec_eigen(1),ECI_vec_eigen(2)};
//     return ECI_vec_output;
// }


void Satellite::add_LVLH_thrust_profile(std::array<double,3> input_LVLH_thrust_vector,double input_thrust_start_time, double input_thrust_end_time){
    ThrustProfileLVLH new_thrust_profile(input_thrust_start_time, input_thrust_end_time, input_LVLH_thrust_vector);
    thrust_profile_list_.push_back(new_thrust_profile);
    if (input_thrust_start_time==0){
        list_of_LVLH_forces_at_this_time_.push_back(input_LVLH_thrust_vector);
        std::array<double,3> ECI_thrust_vector=convert_LVLH_to_ECI_manual(input_LVLH_thrust_vector);
        list_of_ECI_forces_at_this_time_.push_back(ECI_thrust_vector);
    }
}

void Satellite::add_LVLH_thrust_profile(std::array<double,3> input_LVLH_normalized_thrust_direction,double input_LVLH_thrust_magnitude,double input_thrust_start_time, double input_thrust_end_time){
    ThrustProfileLVLH new_thrust_profile(input_thrust_start_time, input_thrust_end_time, input_LVLH_normalized_thrust_direction,input_LVLH_thrust_magnitude);
    thrust_profile_list_.push_back(new_thrust_profile);

    std::array<double,3> LVLH_thrust_vec={0,0,0};

    for (size_t ind=0;ind<3;ind++){
        LVLH_thrust_vec.at(ind)=input_LVLH_thrust_magnitude*input_LVLH_normalized_thrust_direction.at(ind);
    }

    if (input_thrust_start_time==0){
        list_of_LVLH_forces_at_this_time_.push_back(LVLH_thrust_vec);
        std::array<double,3> ECI_thrust_vector=convert_LVLH_to_ECI_manual(LVLH_thrust_vec);
        list_of_ECI_forces_at_this_time_.push_back(ECI_thrust_vector);
    }
}









//"manual" version, via dot products and cross products with position and velocity vectors
std::array<double,3> Satellite::convert_LVLH_to_ECI_manual(std::array<double,3> input_LVLH_vec){
    //LVLH x-axis is defined as in the direction of motion
    //LVLH z-axis is defined as pointing back towards Earth, so along the reversed direction of the position vector from the center of the Earth
    //y-axis determined from a cross product

    Vector3d ECI_unit_vec_x={1,0,0};
    Vector3d ECI_unit_vec_y={0,1,0};
    Vector3d ECI_unit_vec_z={0,0,1};

    std::array<double,3> current_ECI_position_array=get_ECI_position();
    Vector3d current_ECI_position_unit_vec;
    current_ECI_position_unit_vec << current_ECI_position_array.at(0),current_ECI_position_array.at(1),current_ECI_position_array.at(2);
    current_ECI_position_unit_vec.normalize(); //to make it actually a unit vector

    std::array<double,3> current_ECI_velocity_array=get_ECI_velocity();
    Vector3d current_ECI_velocity_unit_vec;
    current_ECI_velocity_unit_vec << current_ECI_velocity_array.at(0),current_ECI_velocity_array.at(1),current_ECI_velocity_array.at(2);
    current_ECI_velocity_unit_vec.normalize();

    Vector3d LVLH_x_unit_vec=current_ECI_velocity_unit_vec;
    Vector3d LVLH_z_unit_vec=(-1)*current_ECI_position_unit_vec;

    Vector3d LVLH_y_unit_vec=LVLH_z_unit_vec.cross(LVLH_x_unit_vec);
    //Should already be normalized, just in case though
    LVLH_y_unit_vec.normalize();

    

    double v_x_ECI= input_LVLH_vec.at(0)*ECI_unit_vec_x.dot(LVLH_x_unit_vec) + input_LVLH_vec.at(1)*ECI_unit_vec_x.dot(LVLH_y_unit_vec) + input_LVLH_vec.at(2)*ECI_unit_vec_x.dot(LVLH_z_unit_vec);

    double v_y_ECI= input_LVLH_vec.at(0)*ECI_unit_vec_y.dot(LVLH_x_unit_vec) + input_LVLH_vec.at(1)*ECI_unit_vec_y.dot(LVLH_y_unit_vec) + input_LVLH_vec.at(2)*ECI_unit_vec_y.dot(LVLH_z_unit_vec);

    double v_z_ECI= input_LVLH_vec.at(0)*ECI_unit_vec_z.dot(LVLH_x_unit_vec) + input_LVLH_vec.at(1)*ECI_unit_vec_z.dot(LVLH_y_unit_vec) + input_LVLH_vec.at(2)*ECI_unit_vec_z.dot(LVLH_z_unit_vec);


    std::array<double,3> output_ECI_arr={v_x_ECI,v_y_ECI,v_z_ECI};
    return output_ECI_arr;
}

std::array<double,3> Satellite::convert_ECI_to_LVLH_manual(std::array<double,3> input_ECI_vec){
    //LVLH x-axis is defined as in the direction of motion
    //LVLH z-axis is defined as pointing back towards Earth, so along the reversed direction of the position vector from the center of the Earth
    //y-axis determined from a cross product

    Vector3d ECI_unit_vec_x={1,0,0};
    Vector3d ECI_unit_vec_y={0,1,0};
    Vector3d ECI_unit_vec_z={0,0,1};

    std::array<double,3> current_ECI_position_array=get_ECI_position();
    Vector3d current_ECI_position_unit_vec;
    current_ECI_position_unit_vec << current_ECI_position_array.at(0),current_ECI_position_array.at(1),current_ECI_position_array.at(2);
    current_ECI_position_unit_vec.normalize(); //to make it actually a unit vector

    std::array<double,3> current_ECI_velocity_array=get_ECI_velocity();
    Vector3d current_ECI_velocity_unit_vec;
    current_ECI_velocity_unit_vec << current_ECI_velocity_array.at(0),current_ECI_velocity_array.at(1),current_ECI_velocity_array.at(2);
    current_ECI_velocity_unit_vec.normalize();

    Vector3d LVLH_x_unit_vec=current_ECI_velocity_unit_vec;
    Vector3d LVLH_z_unit_vec=(-1)*current_ECI_position_unit_vec;

    Vector3d LVLH_y_unit_vec=LVLH_z_unit_vec.cross(LVLH_x_unit_vec);
    //Should already be normalized, just in case though
    LVLH_y_unit_vec.normalize();


    

    double v_x_LVLH= input_ECI_vec.at(0)*LVLH_x_unit_vec.dot(ECI_unit_vec_x) + input_ECI_vec.at(1)*LVLH_x_unit_vec.dot(ECI_unit_vec_y) + input_ECI_vec.at(2)*LVLH_x_unit_vec.dot(ECI_unit_vec_z);

    double v_y_LVLH= input_ECI_vec.at(0)*LVLH_y_unit_vec.dot(ECI_unit_vec_x) + input_ECI_vec.at(1)*LVLH_y_unit_vec.dot(ECI_unit_vec_y) + input_ECI_vec.at(2)*LVLH_y_unit_vec.dot(ECI_unit_vec_z);

    double v_z_LVLH= input_ECI_vec.at(0)*LVLH_z_unit_vec.dot(ECI_unit_vec_x) + input_ECI_vec.at(1)*LVLH_z_unit_vec.dot(ECI_unit_vec_y) + input_ECI_vec.at(2)*LVLH_z_unit_vec.dot(ECI_unit_vec_z);


    std::array<double,3> output_LVLH_arr={v_x_LVLH,v_y_LVLH,v_z_LVLH};
    return output_LVLH_arr;
}

// std::array<double,3> Satellite::convert_body_frame_to_LVLH(std::array<double,3> input_body_frame_vec){
//     //taking LVLH to be the "unrotated" body frame

//     std::array<double,3> LVLH_vec=convert_rotated_body_frame_to_unrotated_body_frame(input_body_frame_vec, theta_, phi_, psi_);
//     return LVLH_vec;
// }

// std::array<double,3> Satellite::convert_body_frame_to_ECI(std::array<double,3> input_body_frame_vec){
//     //Current strategy: first convert rotated body frame to unrotated body frame, which I'm taking to be the LVLH frame
//     //Then convert the unrotated LVLH frame into ECI
//     std::array<double,3> LVLH_vec=convert_body_frame_to_LVLH(input_body_frame_vec);
//     std::array<double,3> ECI_vec=convert_LVLH_to_ECI_manual(LVLH_vec);
// }
