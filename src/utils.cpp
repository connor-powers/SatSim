#include <iostream>
#include <cmath>
#include "Satellite.h"
#include <Eigen/Dense>

using Eigen::Matrix3d;
using Eigen::Vector3d;
using Eigen::Vector4d;
using Eigen::Matrix4d;



//"manual" version, via dot products and cross products with position and velocity vectors
std::array<double,3> convert_LVLH_to_ECI_manual(std::array<double,3> input_LVLH_vec,std::array<double,3> input_position_vec,std::array<double,3> input_velocity_vec){
    //LVLH x-axis is defined as in the direction of motion
    //LVLH z-axis is defined as pointing back towards Earth, so along the reversed direction of the position vector from the center of the Earth
    //y-axis determined from a cross product

    Vector3d ECI_unit_vec_x={1.0,0.0,0.0};
    Vector3d ECI_unit_vec_y={0.0,1.0,0.0};
    Vector3d ECI_unit_vec_z={0.0,0.0,1.0};

    std::array<double,3> current_ECI_position_array=input_position_vec;
    Vector3d current_ECI_position_unit_vec;
    current_ECI_position_unit_vec << current_ECI_position_array.at(0),current_ECI_position_array.at(1),current_ECI_position_array.at(2);
    current_ECI_position_unit_vec.normalize(); //to make it actually a unit vector

    std::array<double,3> current_ECI_velocity_array=input_velocity_vec;
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

std::array<double,3> convert_ECI_to_LVLH_manual(std::array<double,3> input_ECI_vec,std::array<double,3> input_position_vec,std::array<double,3> input_velocity_vec){
    //LVLH x-axis is defined as in the direction of motion
    //LVLH z-axis is defined as pointing back towards Earth, so along the reversed direction of the position vector from the center of the Earth
    //y-axis determined from a cross product

    Vector3d ECI_unit_vec_x={1.0,0.0,0.0};
    Vector3d ECI_unit_vec_y={0.0,1.0,0.0};
    Vector3d ECI_unit_vec_z={0.0,0.0,1.0};

    std::array<double,3> current_ECI_position_array=input_position_vec;
    Vector3d current_ECI_position_unit_vec;
    current_ECI_position_unit_vec << current_ECI_position_array.at(0),current_ECI_position_array.at(1),current_ECI_position_array.at(2);
    current_ECI_position_unit_vec.normalize(); //to make it actually a unit vector

    std::array<double,3> current_ECI_velocity_array=input_velocity_vec;
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



//baselining cartesian coordinates in ECI frame

std::array<double,3> calculate_orbital_acceleration(const std::array<double,3> input_r_vec,const double input_spacecraft_mass,std::vector<std::array<double,3>> input_vec_of_force_vectors_in_ECI)
{
    //orbital acceleration = -G m_Earth/distance^3 * r_vec (just based on rearranging F=ma with a the acceleration due to gravitational attraction between satellite and Earth https://en.wikipedia.org/wiki/Newton%27s_law_of_universal_gravitation)
    //going to be assuming Earth's position doesn't change for now
    //also assuming Earth is spherical, can loosen this assumption in the future
    //note: this is in ECI frame

    std::array<double,3> acceleration_vec_due_to_gravity=input_r_vec; //shouldn't need to explicitly call copy function because input_r_vec is passed by value not ref
    
    //F=ma
    //a=F/m = (F_grav + F_ext)/m = (F_grav/m) + (F_ext/m) = -G*M_Earth/distance^3 + ...

    const double distance=sqrt(input_r_vec.at(0)*input_r_vec.at(0) + input_r_vec.at(1)*input_r_vec.at(1) + input_r_vec.at(2)*input_r_vec.at(2));
    const double overall_factor= -G*mass_Earth/pow(distance,3);



    for (size_t ind=0;ind<input_r_vec.size();ind++){
        acceleration_vec_due_to_gravity.at(ind)*=overall_factor;
    }

    //now add effects from externally-applied forces, e.g., thrusters, if any
    std::array<double,3> acceleration_vec=acceleration_vec_due_to_gravity;

    for (std::array<double,3> external_force_vec_in_ECI : input_vec_of_force_vectors_in_ECI){
        acceleration_vec.at(0)+=(external_force_vec_in_ECI.at(0)/input_spacecraft_mass);
        acceleration_vec.at(1)+=(external_force_vec_in_ECI.at(1)/input_spacecraft_mass);
        acceleration_vec.at(2)+=(external_force_vec_in_ECI.at(2)/input_spacecraft_mass);
    }
    return acceleration_vec;
}



std::array<double,3> calculate_orbital_acceleration(const std::array<double,3> input_r_vec,const double input_spacecraft_mass,std::vector<ThrustProfileLVLH> input_list_of_thrust_profiles_LVLH,double input_evaluation_time,const std::array<double,3> input_velocity_vec)
{
    //orbital acceleration = -G m_Earth/distance^3 * r_vec (just based on rearranging F=ma with a the acceleration due to gravitational attraction between satellite and Earth https://en.wikipedia.org/wiki/Newton%27s_law_of_universal_gravitation)
    //going to be assuming Earth's position doesn't change for now
    //also assuming Earth is spherical, can loosen this assumption in the future
    //note: this is in ECI frame (r_vec and velocity vec should also be in ECI frame)

    std::array<double,3> acceleration_vec_due_to_gravity=input_r_vec; //shouldn't need to explicitly call copy function because input_r_vec is passed by value not ref
    
    //F=ma
    //a=F/m = (F_grav + F_ext)/m = (F_grav/m) + (F_ext/m) = -G*M_Earth/distance^3 + ...

    const double distance=sqrt(input_r_vec.at(0)*input_r_vec.at(0) + input_r_vec.at(1)*input_r_vec.at(1) + input_r_vec.at(2)*input_r_vec.at(2));
    const double overall_factor= -G*mass_Earth/pow(distance,3);



    for (size_t ind=0;ind<input_r_vec.size();ind++){
        acceleration_vec_due_to_gravity.at(ind)*=overall_factor;
    }

    std::array<double,3> acceleration_vec=acceleration_vec_due_to_gravity;

    //now add effects from externally-applied forces, e.g., thrusters, if any

    std::vector<std::array<double,3>> list_of_LVLH_forces_at_evaluation_time={};
    std::vector<std::array<double,3>> list_of_ECI_forces_at_evaluation_time={};

    for (ThrustProfileLVLH thrust_profile : input_list_of_thrust_profiles_LVLH){

        if ((input_evaluation_time>=thrust_profile.t_start_)&&(input_evaluation_time<=thrust_profile.t_end_)){
            list_of_LVLH_forces_at_evaluation_time.push_back(thrust_profile.LVLH_force_vec_);
            std::array<double,3> ECI_thrust_vector=convert_LVLH_to_ECI_manual(thrust_profile.LVLH_force_vec_,input_r_vec,input_velocity_vec);
            list_of_ECI_forces_at_evaluation_time.push_back(ECI_thrust_vector);
        }
    }

    for (std::array<double,3> external_force_vec_in_ECI : list_of_ECI_forces_at_evaluation_time){
        acceleration_vec.at(0)+=(external_force_vec_in_ECI.at(0)/input_spacecraft_mass);
        acceleration_vec.at(1)+=(external_force_vec_in_ECI.at(1)/input_spacecraft_mass);
        acceleration_vec.at(2)+=(external_force_vec_in_ECI.at(2)/input_spacecraft_mass);
    }
    return acceleration_vec;
}









std::array<double,6> RK4_deriv_function_orbit_position_and_velocity(std::array<double,6> input_position_and_velocity,const double input_spacecraft_mass,std::vector<std::array<double,3>> input_vec_of_force_vectors_in_ECI){
    std::array<double,6> derivative_of_input_y={};
    std::array<double,3> position_array={};

    for (size_t ind=0;ind<3;ind++){
        derivative_of_input_y.at(ind)=input_position_and_velocity.at(ind+3);
        position_array.at(ind)=input_position_and_velocity.at(ind);
    }
    
    std::array<double,3> calculated_orbital_acceleration=calculate_orbital_acceleration(position_array,input_spacecraft_mass,input_vec_of_force_vectors_in_ECI);

    for (size_t ind=3;ind<6;ind++){
        derivative_of_input_y.at(ind)=calculated_orbital_acceleration.at(ind-3);
    }

    return derivative_of_input_y;
}


std::array<double,6> RK45_deriv_function_orbit_position_and_velocity(std::array<double,6> input_position_and_velocity,const double input_spacecraft_mass,std::vector<ThrustProfileLVLH> input_list_of_thrust_profiles_LVLH,double input_evaluation_time){
    std::array<double,6> derivative_of_input_y={};
    std::array<double,3> position_array={};
    std::array<double,3> velocity_array={};

    for (size_t ind=0;ind<3;ind++){
        derivative_of_input_y.at(ind)=input_position_and_velocity.at(ind+3);
        velocity_array.at(ind)=input_position_and_velocity.at(ind+3);
        position_array.at(ind)=input_position_and_velocity.at(ind);
    }
    
    std::array<double,3> calculated_orbital_acceleration=calculate_orbital_acceleration(position_array,input_spacecraft_mass,input_list_of_thrust_profiles_LVLH,input_evaluation_time,velocity_array);

    for (size_t ind=3;ind<6;ind++){
        derivative_of_input_y.at(ind)=calculated_orbital_acceleration.at(ind-3);
    }

    return derivative_of_input_y;
}


// std::array<double,6> RK4_deriv_function_angular(std::array<double,6> input_angular_vec,const double input_spacecraft_MOI,std::vector<std::array<double,3>> input_vec_of_body_frame_torque_vectors){
//     std::array<double,6> derivative_of_input_y={};

//     for (size_t ind=0;ind<3;ind++){
//         derivative_of_input_y.at(ind)=input_angular_vec.at(ind+3);
//     }
    
//     std::array<double,3> calculated_angular_acceleration=calculate_body_frame_angular_acceleration(input_spacecraft_MOI,input_vec_of_body_frame_torque_vectors);

//     for (size_t ind=3;ind<6;ind++){
//         derivative_of_input_y.at(ind)=calculated_angular_acceleration.at(ind-3);
//     }

//     return derivative_of_input_y;
// }






void sim_and_draw_orbit_gnuplot(std::vector<Satellite> input_satellite_vector,double input_timestep, double input_total_sim_time, double input_epsilon){
    if (input_satellite_vector.size()<1){
        std::cout << "No input Satellite objects\n";
        return;
    }

    //first, open "pipe" to gnuplot
    FILE *gnuplot_pipe = popen("gnuplot -persist", "w");
    //if it exists
    if (gnuplot_pipe){


        //formatting
        fprintf(gnuplot_pipe,"set xlabel 'x'\n");
        fprintf(gnuplot_pipe,"set ylabel 'y'\n");
        fprintf(gnuplot_pipe,"set zlabel 'z'\n");
        fprintf(gnuplot_pipe,"set title 'Simulated orbits up to time %.2f s'\n",input_total_sim_time);
        // fprintf(gnuplot_pipe,"set view 70,1,1,1\n");        
        fprintf(gnuplot_pipe,"set view equal xyz\n");        

        fprintf(gnuplot_pipe,"unset colorbox\n");    
        fprintf(gnuplot_pipe,"set style fill transparent solid 1.0\n");    

        fprintf(gnuplot_pipe,"set key\n");   
        fprintf(gnuplot_pipe,"set hidden3d front\n");   


        //plotting
        //first let's set the stage for plotting the Earth
        fprintf(gnuplot_pipe,"R_Earth=%f\n",radius_Earth);
        fprintf(gnuplot_pipe,"set isosamples 50,50\n");
        fprintf(gnuplot_pipe,"set parametric\n");
        fprintf(gnuplot_pipe,"set urange [-pi/2:pi/2]\n");
        fprintf(gnuplot_pipe,"set vrange [0:2*pi]\n");


        //first satellite
        Satellite current_satellite=input_satellite_vector.at(0);
        if (input_satellite_vector.size()==1){
            if (current_satellite.plotting_color_.size()>0){
                fprintf(gnuplot_pipe,"splot '-' with lines lw 1 lc rgb '%s' title '%s' \\\n",current_satellite.plotting_color_.c_str(),current_satellite.get_name().c_str());
            }
            else {
                fprintf(gnuplot_pipe,"splot '-' with lines lw 1 title '%s' \\\n",current_satellite.get_name().c_str());
            }
            
        }

        else {
            if (current_satellite.plotting_color_.size()>0){
                fprintf(gnuplot_pipe,"splot '-' with lines lw 1 lc rgb '%s' title '%s'\\\n",current_satellite.plotting_color_.c_str(),current_satellite.get_name().c_str());
            }
            else {
                fprintf(gnuplot_pipe,"splot '-' with lines lw 1 title '%s'\\\n",current_satellite.get_name().c_str());
            }
        }

        for (size_t satellite_index=1;satellite_index<input_satellite_vector.size();satellite_index++){
            current_satellite=input_satellite_vector.at(satellite_index);
            if (satellite_index<input_satellite_vector.size()-1){
                if (current_satellite.plotting_color_.size()>0){
                    fprintf(gnuplot_pipe,",'-' with lines lw 1 lc rgb '%s' title '%s' \\\n",current_satellite.plotting_color_.c_str(),current_satellite.get_name().c_str());
                }
                else {
                    fprintf(gnuplot_pipe,",'-' with lines lw 1 title '%s' \\\n",current_satellite.get_name().c_str());
                }
                
            }

            else {
                if (current_satellite.plotting_color_.size()>0){
                    fprintf(gnuplot_pipe,",'-' with lines lw 1 lc rgb '%s' title '%s'\\\n",current_satellite.plotting_color_.c_str(),current_satellite.get_name().c_str());
                }
                else {
                    fprintf(gnuplot_pipe,",'-' with lines lw 1 title '%s'\\\n",current_satellite.get_name().c_str());
                }
            }
        }
        fprintf(gnuplot_pipe,",R_Earth*cos(u)*cos(v),R_Earth*cos(u)*sin(v),R_Earth*sin(u) notitle with pm3d fillcolor rgbcolor 'navy'\n");

        //now the orbit data, inline, one satellite at a time
        for (size_t satellite_index=0;satellite_index<input_satellite_vector.size();satellite_index++){
            Satellite current_satellite=input_satellite_vector.at(satellite_index);
            std::array<double,3> initial_position=current_satellite.get_ECI_position();
            fprintf(gnuplot_pipe,"%f %f %f\n",initial_position.at(0),initial_position.at(1),initial_position.at(2));
    
    
            std::array<double,3> evolved_position={};
    
            double timestep_to_use=input_timestep;
            double current_satellite_time=current_satellite.get_instantaneous_time();
            while (current_satellite_time<input_total_sim_time){
                double new_timestep=current_satellite.evolve_RK45(input_epsilon,timestep_to_use);
                timestep_to_use=new_timestep;
                evolved_position=current_satellite.get_ECI_position();
                current_satellite_time=current_satellite.get_instantaneous_time();
                fprintf(gnuplot_pipe,"%f %f %f\n",evolved_position.at(0),evolved_position.at(1),evolved_position.at(2));
            }
            fprintf(gnuplot_pipe,"e\n");

        }
        fprintf(gnuplot_pipe,"pause mouse keypress\n");

        fprintf(gnuplot_pipe,"exit \n");
        pclose(gnuplot_pipe);


    }
    else {
        std::cout << "gnuplot not found";
    }

    return;
}




// Matrix3d z_rot_matrix(double input_angle){
//     Matrix3d z_rotation_matrix;
//     z_rotation_matrix << cos(input_angle), -sin(input_angle), 0,
//                  sin(input_angle), cos(input_angle), 0,
//                  0,0,1;

//     return z_rotation_matrix;
// }

// Matrix3d y_rot_matrix(double input_angle){
//     Matrix3d y_rotation_matrix;
//     y_rotation_matrix << cos(input_angle), 0, sin(input_angle),
//                         0,1,0,
//                         -sin(input_angle), 0, cos(input_angle);

//     return y_rotation_matrix;
// }


// Matrix3d x_rot_matrix(double input_angle){
//     Matrix3d x_rotation_matrix;
//     x_rotation_matrix << 1,0,0,
//                          0, cos(input_angle), -sin(input_angle),
//                          0, sin(input_angle), cos(input_angle);

//     return x_rotation_matrix;
// }


// std::array<double,3> convert_rotated_body_frame_to_unrotated_body_frame(std::array<double,3> input_rotated_body_frame_array,double theta, double phi, double psi){
//     Vector3d rotated_body_frame_vec;
//     rotated_body_frame_vec << input_rotated_body_frame_array.at(0), input_rotated_body_frame_array.at(1), input_rotated_body_frame_array.at(2);
//     //build A matrix
//     Matrix3d A_mat;
//     A_mat << cos(phi)*cos(psi) - cos(theta)*sin(phi)*sin(psi), -cos(phi)*sin(psi) - cos(theta)*sin(phi)*cos(psi), sin(theta)*sin(phi),
//              sin(phi)*cos(psi) + cos(theta)*cos(phi)*sin(psi), -sin(phi)*sin(psi) + cos(theta)*cos(phi)*cos(psi), -sin(theta)*cos(phi),
//              sin(theta)*sin(psi), sin(theta)*cos(psi), cos(theta);

//     Vector3d unrotated_body_frame_vec=A_mat*rotated_body_frame_vec;
//     std::array<double,3> unrotated_body_frame_array;
//     for (size_t ind=0;ind<3;ind++){
//         unrotated_body_frame_array.at(ind)=unrotated_body_frame_vec(ind);
//     }
//     return unrotated_body_frame_array;
// }




// std::array<double,3> calculate_body_frame_angular_acceleration(const double input_spacecraft_moi,std::vector<std::array<double,3>> input_vec_of_torque_vectors_in_body_frame)
// {
//     //Body-frame angular acceleration calculated from alpha_i=tau_i/I for torque tau, MOI I

//     std::array<double,3> body_frame_angular_acceleration={0,0,0};

//     for (std::array<double,3> body_frame_torque_vec : input_vec_of_torque_vectors_in_body_frame){
//         for (size_t coord_ind=0;coord_ind<3;coord_ind++){
//             acceleration_vec.at(coord_ind)+=(body_frame_torque_vec.at(coord_ind)/input_spacecraft_moi);
//         }

//     }
//     return body_frame_angular_acceleration;
// }




std::array<double,4> bodyframe_quaternion_deriv(std::array<double,4> input_bodyframe_quaternion,double input_w_1, double input_w_2, double input_w_3){

    //Assumes quaternions are input as [q_0, q_1, q_2, q_3] and omega is the spacecraft "body rate", represented in the body frame, w.r.t. a chosen reference frame
    //Ref: https://ntrs.nasa.gov/api/citations/20240009554/downloads/Space%20Attitude%20Development%20Control.pdf

    //Let's do this with Eigen so we can use matrix multiplication easily
    Vector4d input_body_quaternion;
    input_body_quaternion << input_bodyframe_quaternion.at(0),input_bodyframe_quaternion.at(1),input_bodyframe_quaternion.at(2),input_bodyframe_quaternion.at(3);

    Matrix4d coeff_mat;
    coeff_mat << 0,-input_w_1,-input_w_2,-input_w_3,
                 input_w_1, 0, input_w_3, -input_w_2,
                 input_w_2, -input_w_3, 0, input_w_1,
                 input_w_3, input_w_2, -input_w_1, 0;

    Vector4d quaternion_deriv_vec= 0.5*coeff_mat*input_body_quaternion;
    std::array<double,4> quaternion_deriv_array;
    for (size_t ind=0;ind<quaternion_deriv_array.size();ind++){
        quaternion_deriv_array.at(ind)=quaternion_deriv_vec(ind);
    }
    return quaternion_deriv_array;
}   
