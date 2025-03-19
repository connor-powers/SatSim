#include <iostream>
#include "Satellite.h"
#include "utils.h"

int main()
{
    Satellite test_sat_1("../input.json");

    // double epsilon=pow(10,-9);
    // double test_timestep=1;
    // for (size_t ind=0;ind<10;ind++){
    //     std::cout << "Timestep index: " << ind << "\n";
    //     std::pair<double,int> output_pair=test_sat_1.evolve_RK45(epsilon,test_timestep);
    //     double next_timestep=output_pair.first;
    //     std::cout << "next timestep: " << next_timestep << "\n";
    //     int error_code=output_pair.second;
    //     test_timestep=next_timestep;
    //     std::cout << "Roll: " << test_sat_1.get_attitude_val("Roll") << "\n";
    //     std::cout << "Pitch: " << test_sat_1.get_attitude_val("Pitch")<< "\n";
    //     std::cout << "Yaw: " << test_sat_1.get_attitude_val("Yaw")<< "\n";
    //     std::cout << "omega_x: " << test_sat_1.get_attitude_val("omega_x")<< "\n";
    //     std::cout << "omega_y: " << test_sat_1.get_attitude_val("omega_y")<< "\n";
    //     std::cout << "omega_z: " << test_sat_1.get_attitude_val("omega_z")<< "\n";
    //     std::cout << "q_0: " << test_sat_1.get_attitude_val("q_0")<< "\n";
    //     std::cout << "q_1: " << test_sat_1.get_attitude_val("q_1")<< "\n";
    //     std::cout << "q_2: " << test_sat_1.get_attitude_val("q_2")<< "\n";
    //     std::cout << "q_3: " << test_sat_1.get_attitude_val("q_3")<< "\n";
    // }
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
    double epsilon=pow(10,-12);
    // sim_and_draw_orbit_gnuplot(satellite_vector_1,timestep,total_sim_time,epsilon);

    //Now some demonstrations of plotting orbital parameters
    Satellite test_sat_4("../input_4.json");
    Satellite test_sat_5("../input_5.json");
    Satellite test_sat_6("../input_6.json");


    std::vector<Satellite> satellite_vector_2={test_sat_4,test_sat_5,test_sat_6};
    total_sim_time=9952;
    // sim_and_plot_orbital_elem_gnuplot(satellite_vector_2,timestep,total_sim_time,epsilon,"Argument of Periapsis");
    Satellite test_sat_7("../input_7.json");

    std::vector<Satellite> satellite_vector_3={test_sat_7};
        sim_and_draw_orbit_gnuplot(satellite_vector_3,timestep,total_sim_time,epsilon);

    // sim_and_plot_orbital_elem_gnuplot(satellite_vector_3,timestep,total_sim_time,epsilon,"True Anomaly");
    sim_and_plot_orbital_elem_gnuplot(satellite_vector_3,timestep,total_sim_time,epsilon,"Eccentricity");

    // sim_and_plot_orbital_elem_gnuplot(satellite_vector_3,timestep,total_sim_time,epsilon,"Orbital Rate");
    // sim_and_plot_orbital_elem_gnuplot(satellite_vector_3,timestep,total_sim_time,epsilon,"Orbital Angular Acceleration");

    // sim_and_plot_attitude_evolution_gnuplot(satellite_vector_3,timestep,total_sim_time,epsilon,"omega_y",false);

    // sim_and_plot_attitude_evolution_gnuplot(satellite_vector_3,timestep,total_sim_time,epsilon,"q_0",false);
    // sim_and_plot_attitude_evolution_gnuplot(satellite_vector_3,timestep,total_sim_time,epsilon,"q_1",false);
    // sim_and_plot_attitude_evolution_gnuplot(satellite_vector_3,timestep,total_sim_time,epsilon,"q_2",false);
    // sim_and_plot_attitude_evolution_gnuplot(satellite_vector_3,timestep,total_sim_time,epsilon,"q_3",false);
    // sim_and_plot_attitude_evolution_gnuplot(satellite_vector_3,timestep,total_sim_time,epsilon,"Pitch",false);
    // sim_and_plot_attitude_evolution_gnuplot(satellite_vector_3,timestep,total_sim_time,epsilon,"Yaw",false);
    // sim_and_plot_attitude_evolution_gnuplot(satellite_vector_3,timestep,total_sim_time,epsilon,"Roll",false);

    return 0;
}