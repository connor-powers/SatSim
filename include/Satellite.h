#include <iostream>

class Satellite
{
    private:
        double inclination_;
        double raan_;
        double arg_of_periapsis_;
        double eccentricity_;
        double a_;
        double true_anomaly_;

        void propagate_orbit(double input_propagation_time);

    public:
        Satellite(std::string input_file_name){
            //baselining JSON input file format specifying initial orbital parameters of satellite

        }

        void set_orbital_altitude(double input_orbital_altitude);
};