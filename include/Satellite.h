#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

using json=nlohmann::json;

//Define constants
const double G=6.674*pow(10,-11);
const double mass_Earth=5.9722*pow(10,24);

class Satellite
{
    private:
        double inclination_;
        double raan_;
        double arg_of_periapsis_;
        double eccentricity_;
        double a_;
        double true_anomaly_;
        double orbital_period_;

        std::array<double,3> instantaneous_position_;
        std::array<double,3> instantaneous_velocity_;

        std::pair<std::array<double,3>,std::array<double,3>> calculate_position_and_velocity_from_orbit_params(const double input_semimajor_axis,const double input_eccentricity,const double input_true_anomaly,const double input_RAAN,const double input_i,const double input_arg_of_periapsis,const double input_orbital_period);
        std::pair<double,double> calculate_eccentric_anomaly(const double input_eccentricity, const double input_true_anomaly,const double input_semimajor_axis);
        double calculate_orbital_period(double input_semimajor_axis);
        std::array<double,3> transform_orbital_plane_coords_to_3D_ECI_cartesian(double input_x, double input_y,double input_RAAN,double input_i,double input_arg_of_periapsis);

    public:
        Satellite(std::string input_file_name){
            //baselining JSON input file format specifying initial orbital parameters of satellite
            //semimajor axis is read in units of km
            std::ifstream input_filestream(input_file_name);
            json input_data=json::parse(input_filestream);

            inclination_=input_data["Inclination"];
            raan_=input_data["RAAN"];
            arg_of_periapsis_=input_data["Argument of Periapsis"];
            eccentricity_=input_data["Eccentricity"];
            a_=input_data["Semimajor Axis"];
            a_*=1000; //converting from km to m
            true_anomaly_=input_data["True Anomaly"];
            orbital_period_=calculate_orbital_period(a_);

            std::pair<std::array<double,3>,std::array<double,3>> initial_position_and_vel=calculate_position_and_velocity_from_orbit_params(a_,eccentricity_,true_anomaly_,raan_,inclination_,arg_of_periapsis_,orbital_period_);
            std::array<double,3> initial_cartesian_ECI_position=initial_position_and_vel.first;
            std::array<double,3> initial_cartesian_ECI_velocity=initial_position_and_vel.second;
            instantaneous_position_=initial_cartesian_ECI_position;
            instantaneous_velocity_=initial_cartesian_ECI_velocity;

        }

        std::array<double,3> get_position(){
            return instantaneous_position_;
        }
        std::array<double,3> get_velocity(){
            return instantaneous_velocity_;
        }

};