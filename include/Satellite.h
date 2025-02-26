#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

using json=nlohmann::json;

//Define constants
const double G=6.674*pow(10,-11); //https://en.wikipedia.org/wiki/Gravitational_constant
const double mass_Earth=5.9722*pow(10,24); //https://en.wikipedia.org/wiki/Earth_mass

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
        double m_;
        double t_;

        std::array<double,3> instantaneous_position_;
        std::array<double,3> instantaneous_velocity_;

        double instantaneous_orbital_radius_;
        double instantaneous_orbital_speed_;

        std::pair<std::array<double,3>,std::array<double,3>> calculate_position_and_velocity_from_orbit_params(const double input_semimajor_axis,const double input_eccentricity,const double input_true_anomaly,const double input_RAAN,const double input_i,const double input_arg_of_periapsis,const double input_orbital_period);
        std::pair<double,double> calculate_eccentric_anomaly(const double input_eccentricity, const double input_true_anomaly,const double input_semimajor_axis);
        double calculate_orbital_period(double input_semimajor_axis);
        std::array<double,3> transform_orbital_plane_coords_to_3D_ECI_cartesian(double input_x, double input_y,double input_RAAN,double input_i,double input_arg_of_periapsis);
    public:
        Satellite(std::string input_file_name){
            //baselining JSON input file format specifying initial orbital parameters of satellite
            //semimajor axis is read in units of km
            //angles are read in units of degrees, then internally translated to radians
            std::ifstream input_filestream(input_file_name);
            json input_data=json::parse(input_filestream);

            inclination_=input_data["Inclination"];
            //convert to radians
            inclination_*=(M_PI/180);

            raan_=input_data["RAAN"];
            //convert to radians
            raan_*=(M_PI/180);

            arg_of_periapsis_=input_data["Argument of Periapsis"];
            //convert to radians
            arg_of_periapsis_*=(M_PI/180);

            eccentricity_=input_data["Eccentricity"];

            a_=input_data["Semimajor Axis"];
            a_*=1000; //converting from km to m

            true_anomaly_=input_data["True Anomaly"];
            //convert to radians
            true_anomaly_*=(M_PI/180);

            m_=input_data["Mass"];
            t_=0; //for now, assuming satellites are initialized at time t=0;

            orbital_period_=calculate_orbital_period(a_);

            std::pair<std::array<double,3>,std::array<double,3>> initial_position_and_vel=calculate_position_and_velocity_from_orbit_params(a_,eccentricity_,true_anomaly_,raan_,inclination_,arg_of_periapsis_,orbital_period_);
            std::array<double,3> initial_cartesian_ECI_position=initial_position_and_vel.first;
            std::array<double,3> initial_cartesian_ECI_velocity=initial_position_and_vel.second;
            instantaneous_position_=initial_cartesian_ECI_position;
            instantaneous_velocity_=initial_cartesian_ECI_velocity;

            instantaneous_orbital_radius_=sqrt(pow(instantaneous_position_.at(0),2)+pow(instantaneous_position_.at(1),2)+pow(instantaneous_position_.at(2),2));
            instantaneous_orbital_speed_=sqrt(pow(instantaneous_velocity_.at(0),2)+pow(instantaneous_velocity_.at(1),2)+pow(instantaneous_velocity_.at(2),2));

        }

        std::array<double,3> get_position(){
            return instantaneous_position_;
        }
        std::array<double,3> get_velocity(){
            return instantaneous_velocity_;
        }
        double get_speed(){
            return instantaneous_orbital_speed_;
        }
        double get_radius(){
            return instantaneous_orbital_radius_;
        }
        double get_total_energy(){
            double orbital_radius=get_radius();
            double gravitational_potential_energy=-G*mass_Earth*m_/orbital_radius;

            double orbital_speed=get_speed();
            double kinetic_energy=(1/2)*m_*(orbital_speed*orbital_speed);

            return (gravitational_potential_energy+kinetic_energy);
        }

        void evolve_RK4(double input_timestep, int input_num_timesteps);


};