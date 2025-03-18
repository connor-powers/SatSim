#include <iostream>
#include "Satellite.h"
#include "utils.h"

int main()
{
    Satellite test_sat_1("../input.json");

    std::array<double,3> LVLH_thrust_direction={1,0,0};
    double thrust_magnitude=100; //N
    double t_thrust_start=5000;
    double t_thrust_end=6000;
    test_sat_1.add_LVLH_thrust_profile(LVLH_thrust_direction,thrust_magnitude,t_thrust_start,t_thrust_end);

    std::array<double,3> LVLH_thrust_vec_2={0.51,20,-5};
    double t_thrust_2_start=5500;
    double t_thrust_2_end=6500;
    test_sat_1.add_LVLH_thrust_profile(LVLH_thrust_vec_2,t_thrust_2_start,t_thrust_2_end);
    

    Satellite test_sat_2("../input_2.json");

    Satellite test_sat_3("../input_3.json");
    
    std::vector<Satellite> satellite_vector_1={test_sat_1, test_sat_2,test_sat_3};

    double timestep=2;
    double total_sim_time=25000;
    double epsilon=pow(10,-11);
    sim_and_draw_orbit_gnuplot(satellite_vector_1,timestep,total_sim_time,epsilon);

    //Now some demonstrations of plotting orbital parameters
    Satellite test_sat_4("../input_4.json");
    Satellite test_sat_5("../input_5.json");
    Satellite test_sat_6("../input_6.json");


    std::vector<Satellite> satellite_vector_2={test_sat_4,test_sat_5,test_sat_6};
    total_sim_time=50000;
    sim_and_plot_orbital_elem_gnuplot(satellite_vector_2,timestep,total_sim_time,epsilon,"Argument of Periapsis");
    return 0;
}