#ifndef SATELLITE_HEADER
#define SATELLITE_HEADER

#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

using json=nlohmann::json;

//Define constants
const double G=6.674*pow(10,-11); //https://en.wikipedia.org/wiki/Gravitational_constant
const double mass_Earth=5.9722*pow(10,24); //https://en.wikipedia.org/wiki/Earth_mass
const double radius_Earth=6378137; //https://en.wikipedia.org/wiki/Earth_radius

class ThrustProfileLVLH
{
    //Note: for now, thrust forces are assumed to act through center of mass of satellite.
    public:
        double t_start_={0};
        double t_end_={0};
        std::array<double,3> LVLH_force_vec_={0,0,0};
        ThrustProfileLVLH(double t_start, double t_end, std::array<double,3> LVLH_force_vec){
            t_start_=t_start;
            t_end_=t_end;
            LVLH_force_vec_=LVLH_force_vec;
        }
        ThrustProfileLVLH(double t_start, double t_end, std::array<double,3> LVLH_normalized_force_direction_vec,double input_force_magnitude){
            t_start_=t_start;
            t_end_=t_end;
            for (size_t ind=0;ind<3;ind++){
                LVLH_force_vec_.at(ind)=input_force_magnitude*LVLH_normalized_force_direction_vec.at(ind);
            }
        }
};

class Satellite
{
    private:
        double inclination_={0};
        double raan_={0}; //Assuming RAAN can be used interchangeably with longitude of ascending node for the Earth-orbiting satellites simulated here
        double arg_of_periapsis_={0};
        double eccentricity_={0};
        double a_={0};
        double true_anomaly_={0};
        double orbital_period_={0};
        double m_={1}; //default value to prevent infinities in acceleration calculations from a=F/m
        // double I_={1}; //moment of inertia, taken to be same for all 3 principal axes, set to default value for same reasons as mass
        double t_={0};


        //Now body-frame attributes
        // double pitch_angle_={0};
        // double roll_angle_={0};
        // double yaw_angle_={0};

        // double theta_={0};
        // double phi_={0};
        // double psi_={0};

        std::string name_="";

        std::array<double,3> perifocal_position_={0,0,0};
        std::array<double,3> perifocal_velocity_={0,0,0};

        std::array<double,3> ECI_position_={0,0,0};
        std::array<double,3> ECI_velocity_={0,0,0};

        std::vector<ThrustProfileLVLH> thrust_profile_list_={};

        std::vector<std::array<double,3>> list_of_LVLH_forces_at_this_time_={};
        std::vector<std::array<double,3>> list_of_ECI_forces_at_this_time_={};
        // std::vector<std::array<double,3>> list_of_body_frame_torques_at_this_time_={};

        std::pair<double,double> calculate_eccentric_anomaly(const double input_eccentricity, const double input_true_anomaly,const double input_semimajor_axis);
        double calculate_orbital_period(double input_semimajor_axis);
    public:
        std::string plotting_color_="";
        Satellite(std::string input_file_name){
            //baselining JSON input file format specifying initial orbital parameters of satellite
            //semimajor axis is read in units of km
            //angles are read in units of degrees, then internally translated to radians
            std::ifstream input_filestream(input_file_name);
            json input_data=json::parse(input_filestream);

            inclination_=input_data.at("Inclination");
            //convert to radians
            inclination_*=(M_PI/180);

            raan_=input_data.at("RAAN");
            //convert to radians
            raan_*=(M_PI/180);

            arg_of_periapsis_=input_data.at("Argument of Periapsis");
            //convert to radians
            arg_of_periapsis_*=(M_PI/180);


            eccentricity_=input_data.at("Eccentricity");
            //If circular orbit, arg of periapsis is undefined, using convention of setting it to 0 in this case
            if (eccentricity_==0){
                arg_of_periapsis_=0;
            }


            a_=input_data.at("Semimajor Axis");
            a_*=1000; //converting from km to m

            true_anomaly_=input_data.at("True Anomaly");
            //convert to radians
            true_anomaly_*=(M_PI/180);


            // pitch_angle_=input_data.at("Pitch Angle");
            // //convert to radians
            // pitch_angle_*=(M_PI/180);

            // roll_angle_=input_data.at("Roll Angle");
            // //convert to radians
            // roll_angle_*=(M_PI/180);


            // yaw_angle_=input_data.at("Yaw Angle");
            // //convert to radians
            // yaw_angle_*=(M_PI/180);


            m_=input_data.at("Mass");
            name_=input_data.at("Name");

            //making plotting color an optional parameter
            if (input_data.find("Plotting Color")!=input_data.end()){
                plotting_color_=input_data.at("Plotting Color");
            }
            
            t_=0; //for now, assuming satellites are initialized at time t=0;

            orbital_period_=calculate_orbital_period(a_);

            //updated workflow
            perifocal_position_=calculate_perifocal_position();
            perifocal_velocity_=calculate_perifocal_velocity();

            ECI_position_=convert_perifocal_to_ECI(perifocal_position_);
            ECI_velocity_=convert_perifocal_to_ECI(perifocal_velocity_);


        }

        std::array<double,3> get_ECI_position(){
            return ECI_position_;
        }
        std::array<double,3> get_ECI_velocity(){
            return ECI_velocity_;
        }
        double get_speed(){
            //shouldn't matter which frame I use, might as well use perifocal coords since it's fewer operations (no W-direction component so can omit that term, whereas there's x,y,z components in ECI)
            return sqrt(pow(perifocal_velocity_.at(0),2)+pow(perifocal_velocity_.at(1),2));
        }
        double get_radius(){
            //shouldn't matter which frame I use, might as well use perifocal coords since it's fewer operations (no W-direction component so can omit that term, whereas there's x,y,z components in ECI)
            return sqrt(pow(perifocal_position_.at(0),2)+pow(perifocal_position_.at(1),2));
        }
        double get_total_energy(){
            double orbital_radius=get_radius();
            double gravitational_potential_energy=-G*mass_Earth*m_/orbital_radius;

            double orbital_speed=get_speed();
            double kinetic_energy=(1/2)*m_*(orbital_speed*orbital_speed);

            return (gravitational_potential_energy+kinetic_energy);
        }

        double get_instantaneous_time(){
            return t_;
        }
        std::string get_name(){
            return name_;
        }
        void evolve_RK4(double input_timestep);

        std::array<double,3> body_frame_to_ECI(std::array<double,3> input_vector);

        std::array<double,3> ECI_to_body_frame(std::array<double,3> input_vector);

        std::array<double,3> calculate_perifocal_position();

        std::array<double,3> calculate_perifocal_velocity();

        std::array<double,3> convert_perifocal_to_ECI(std::array<double,3> input_perifocal_vec);
        std::array<double,3> convert_ECI_to_perifocal(std::array<double,3> input_ECI_vec);

        // std::array<double,3> convert_LVLH_to_ECI(std::array<double,3> input_LVLH_vec);

        void add_LVLH_thrust_profile(std::array<double,3> input_LVLH_normalized_thrust_direction,double input_LVLH_thrust_magnitude,double input_thrust_start_time, double input_thrust_end_time);
        void add_LVLH_thrust_profile(std::array<double,3> input_LVLH_thrust_vector,double input_thrust_start_time, double input_thrust_end_time);
        


        // std::array<double,3> convert_body_frame_to_LVLH(std::array<double,3> input_body_frame_vec);

        // std::array<double,3> convert_body_frame_to_LVLH(std::array<double,3> input_body_frame_vec);

        // std::array<double,3> convert_body_frame_to_ECI(std::array<double,3> input_body_frame_vec);

        int update_orbital_elements_from_position_and_velocity();
        std::array<double,6> get_orbital_elements();

        std::pair<double,int> evolve_RK45(double input_epsilon,double input_initial_timestep,bool perturbation=true);

        double get_orbital_element(std::string orbital_element_name);
};





#endif