#include <gtest/gtest.h>

#include <iostream>

#include "Satellite.h"
#include "utils.h"

double epsilon = pow(10, -7);

TEST(UtilsTests, OrbitalElementPlottingTests) {
    // Going to have one satellite with two thrust profiles and two torque profiles (added in different ways),
    // one satellite with just two thrust profiles, one with just two torque profiles
    // and one with none
    Satellite test_sat_both("../tests/elliptical_orbit_test_1.json");
    // Define parameters for an LVLH frame thrust profile
    std::array<double, 3> LVLH_thrust_direction = {1, 0, 0};
    double thrust_magnitude = 100;  // N
    double t_thrust_start = 0;
    double t_thrust_end = 100;
    // Add the thrust profile to the satellite object
    test_sat_both.add_LVLH_thrust_profile(LVLH_thrust_direction, thrust_magnitude,
                                       t_thrust_start, t_thrust_end);
    
    std::array<double, 3> LVLH_thrust_vec_2 = {0.51, 20, -5};
    double t_thrust_2_start = 50;
    double t_thrust_2_end = 150;
    test_sat_both.add_LVLH_thrust_profile(LVLH_thrust_vec_2, t_thrust_2_start,
                                        t_thrust_2_end);


    std::array<double, 3> torque_direction = {0, -1, 0};
    double torque_magnitude = 0.0001;  // N
    double t_torque_start = 101;
    double t_torque_end = 103;
    test_sat_both.add_bodyframe_torque_profile(torque_direction, torque_magnitude,
                                            t_torque_start, t_torque_end);

    std::array<double,3> another_torque_vec = {-0.00001,0.00002,0.00003};
    test_sat_both.add_bodyframe_torque_profile(another_torque_vec,t_torque_start,t_torque_end);
                                        
    Satellite test_sat_thrust("../tests/elliptical_orbit_test_1.json");
    test_sat_thrust.add_LVLH_thrust_profile(LVLH_thrust_direction, thrust_magnitude,
        t_thrust_start, t_thrust_end);
    test_sat_thrust.add_LVLH_thrust_profile(LVLH_thrust_vec_2, t_thrust_2_start,
        t_thrust_2_end);

    Satellite test_sat_torque("../tests/elliptical_orbit_test_4.json");
    test_sat_torque.add_bodyframe_torque_profile(torque_direction, torque_magnitude,
        t_torque_start, t_torque_end);
    test_sat_torque.add_bodyframe_torque_profile(another_torque_vec,t_torque_start,t_torque_end);

    Satellite test_sat_none_1("../tests/elliptical_orbit_test_1.json");
    Satellite test_sat_none_2("../tests/elliptical_orbit_test_4.json");

    std::vector<Satellite> satellite_vector_1 = {test_sat_none_1,test_sat_none_2,test_sat_thrust,test_sat_torque,test_sat_both};
    std::vector<Satellite> satellite_vector_2 = {test_sat_none_1,test_sat_none_2,test_sat_thrust,test_sat_torque,test_sat_both,test_sat_none_2};
    std::vector<Satellite> satellite_vector_3 = {test_sat_none_2,test_sat_none_1,test_sat_thrust,test_sat_torque,test_sat_both};
    std::vector<Satellite> satellite_vector_4 = {test_sat_none_2,test_sat_none_1,test_sat_thrust,test_sat_torque,test_sat_both,test_sat_none_2};

    std::vector<std::vector<Satellite>> satellite_vector_vector = {satellite_vector_1,satellite_vector_2,satellite_vector_3,satellite_vector_4};
    std::vector<Satellite> single_satellite_vec_1 = {test_sat_both};
    std::vector<Satellite> single_satellite_vec_2 = {test_sat_torque};
    satellite_vector_vector.push_back(single_satellite_vec_1);
    satellite_vector_vector.push_back(single_satellite_vec_2);
    double timestep = 2;
    double total_sim_time = 300;

    // Drag parameters
    double F_10 = 100;  // Solar radio ten centimeter flux
    double A_p = 120;   // Geomagnetic A_p index

    // Collect drag parameters into a tuple with F_10 first and A_p second
    std::pair<double, double> drag_elements = {F_10, A_p};
    std::string file_name = "test_plot";
    std::vector<std::string> orbital_elements = {"Semimajor Axis",
        "Eccentricity","Inclination","RAAN",
        "Argument of Periapsis","True Anomaly",
        "Orbital Rate","Orbital Angular Acceleration",
        "Total Energy"};
    for (std::vector<Satellite> satellite_vector : satellite_vector_vector){
        for (std::string element_name : orbital_elements){
            sim_and_plot_orbital_elem_gnuplot(satellite_vector, timestep,
                total_sim_time, epsilon, element_name, file_name,
                false, false);
            sim_and_plot_orbital_elem_gnuplot(satellite_vector, timestep,
                total_sim_time, epsilon, element_name, file_name,
                true, false);
            sim_and_plot_orbital_elem_gnuplot(satellite_vector, timestep,
                total_sim_time, epsilon, element_name, file_name,
                true, true, drag_elements);
        }
    }

  }


// Now a similar test for attitude plotting
TEST(UtilsTests, AttitudeElementPlottingTests) {
    // Going to have one satellite with two thrust profiles and two torque profiles (added in different ways),
    // one satellite with just two thrust profiles, one with just two torque profiles
    // and one with none
    Satellite test_sat_both("../tests/elliptical_orbit_test_1.json");
    // Define parameters for an LVLH frame thrust profile
    std::array<double, 3> LVLH_thrust_direction = {1, 0, 0};
    double thrust_magnitude = 100;  // N
    double t_thrust_start = 50;
    double t_thrust_end = 150;
    // Add the thrust profile to the satellite object
    test_sat_both.add_LVLH_thrust_profile(LVLH_thrust_direction, thrust_magnitude,
                                       t_thrust_start, t_thrust_end);
    
    std::array<double, 3> LVLH_thrust_vec_2 = {0.51, 20, -5};
    double t_thrust_2_start = 0;
    double t_thrust_2_end = 100;
    test_sat_both.add_LVLH_thrust_profile(LVLH_thrust_vec_2, t_thrust_2_start,
                                        t_thrust_2_end);


    std::array<double, 3> torque_direction = {0, -1, 0};
    double torque_magnitude = 0.0001;  // N
    double t_torque_start = 101;
    double t_torque_end = 103;
    test_sat_both.add_bodyframe_torque_profile(torque_direction, torque_magnitude,
                                            t_torque_start, t_torque_end);

    std::array<double,3> another_torque_vec = {-0.00001,0.00002,0.00003};
    test_sat_both.add_bodyframe_torque_profile(another_torque_vec,0,2);
                                        
    Satellite test_sat_thrust("../tests/elliptical_orbit_test_1.json");
    test_sat_thrust.add_LVLH_thrust_profile(LVLH_thrust_direction, thrust_magnitude,
        t_thrust_start, t_thrust_end);
    test_sat_thrust.add_LVLH_thrust_profile(LVLH_thrust_vec_2, t_thrust_2_start,
        t_thrust_2_end);

    Satellite test_sat_torque("../tests/elliptical_orbit_test_4.json");
    test_sat_torque.add_bodyframe_torque_profile(torque_direction, torque_magnitude,
        t_torque_start, t_torque_end);
    test_sat_torque.add_bodyframe_torque_profile(another_torque_vec,t_torque_start,t_torque_end);

    Satellite test_sat_none_1("../tests/elliptical_orbit_test_1.json");
    Satellite test_sat_none_2("../tests/elliptical_orbit_test_4.json");

    std::vector<Satellite> satellite_vector_1 = {test_sat_none_1,test_sat_none_2,test_sat_thrust,test_sat_torque,test_sat_both};
    std::vector<Satellite> satellite_vector_2 = {test_sat_none_1,test_sat_none_2,test_sat_thrust,test_sat_torque,test_sat_both,test_sat_none_2};
    std::vector<Satellite> satellite_vector_3 = {test_sat_none_2,test_sat_none_1,test_sat_thrust,test_sat_torque,test_sat_both};
    std::vector<Satellite> satellite_vector_4 = {test_sat_none_2,test_sat_none_1,test_sat_thrust,test_sat_torque,test_sat_both,test_sat_none_2};

    std::vector<std::vector<Satellite>> satellite_vector_vector = {satellite_vector_1,satellite_vector_2,satellite_vector_3,satellite_vector_4};
    std::vector<Satellite> single_satellite_vec_1 = {test_sat_both};
    std::vector<Satellite> single_satellite_vec_2 = {test_sat_torque};
    satellite_vector_vector.push_back(single_satellite_vec_1);
    satellite_vector_vector.push_back(single_satellite_vec_2);

    double timestep = 2;
    double total_sim_time = 300;

    // Drag parameters
    double F_10 = 100;  // Solar radio ten centimeter flux
    double A_p = 120;   // Geomagnetic A_p index

    // Collect drag parameters into a tuple with F_10 first and A_p second
    std::pair<double, double> drag_elements = {F_10, A_p};
    std::string file_name = "test_plot";
    std::vector<std::string> orbital_elements = {"Roll",
        "Pitch","Yaw","omega_x",
        "omega_y","omega_z",
        "q_0","q_1","q_2","q_3"};
    for (std::vector<Satellite> satellite_vector : satellite_vector_vector){
        for (std::string element_name : orbital_elements){
            sim_and_plot_attitude_evolution_gnuplot(satellite_vector, timestep,
                total_sim_time, epsilon, element_name, file_name,
                false, false);
            sim_and_plot_attitude_evolution_gnuplot(satellite_vector, timestep,
                total_sim_time, epsilon, element_name, file_name,
                true, false);
            sim_and_plot_attitude_evolution_gnuplot(satellite_vector, timestep,
                total_sim_time, epsilon, element_name, file_name,
                true, true, drag_elements);
        }
    }
  }

  // Now one for the 3D orbital plot
  TEST(UtilsTests, ThreeDimensionalOrbitPlotTest) {
    // Going to have one satellite with two thrust profiles and two torque profiles (added in different ways),
    // one satellite with just two thrust profiles, one with just two torque profiles
    // and one with none
    Satellite test_sat_both("../tests/elliptical_orbit_test_1.json");
    // Define parameters for an LVLH frame thrust profile
    std::array<double, 3> LVLH_thrust_direction = {1, 0, 0};
    double thrust_magnitude = 100;  // N
    double t_thrust_start = 0;
    double t_thrust_end = 100;
    // Add the thrust profile to the satellite object
    test_sat_both.add_LVLH_thrust_profile(LVLH_thrust_direction, thrust_magnitude,
                                       t_thrust_start, t_thrust_end);
    
    std::array<double, 3> LVLH_thrust_vec_2 = {0.51, 20, -5};
    double t_thrust_2_start = 50;
    double t_thrust_2_end = 150;
    test_sat_both.add_LVLH_thrust_profile(LVLH_thrust_vec_2, t_thrust_2_start,
                                        t_thrust_2_end);


    std::array<double, 3> torque_direction = {0, -1, 0};
    double torque_magnitude = 0.0001;  // N
    double t_torque_start = 101;
    double t_torque_end = 103;
    test_sat_both.add_bodyframe_torque_profile(torque_direction, torque_magnitude,
                                            t_torque_start, t_torque_end);

    std::array<double,3> another_torque_vec = {-0.00001,0.00002,0.00003};
    test_sat_both.add_bodyframe_torque_profile(another_torque_vec,t_torque_start,t_torque_end);
                                        
    Satellite test_sat_thrust("../tests/elliptical_orbit_test_1.json");
    test_sat_thrust.add_LVLH_thrust_profile(LVLH_thrust_direction, thrust_magnitude,
        t_thrust_start, t_thrust_end);
    test_sat_thrust.add_LVLH_thrust_profile(LVLH_thrust_vec_2, t_thrust_2_start,
        t_thrust_2_end);

    Satellite test_sat_torque("../tests/elliptical_orbit_test_4.json");
    test_sat_torque.add_bodyframe_torque_profile(torque_direction, torque_magnitude,
        t_torque_start, t_torque_end);
    test_sat_torque.add_bodyframe_torque_profile(another_torque_vec,t_torque_start,t_torque_end);

    Satellite test_sat_none_1("../tests/elliptical_orbit_test_1.json");
    Satellite test_sat_none_2("../tests/elliptical_orbit_test_4.json");

    std::vector<Satellite> satellite_vector_1 = {test_sat_none_1,test_sat_none_2,test_sat_thrust,test_sat_torque,test_sat_both};
    std::vector<Satellite> satellite_vector_2 = {test_sat_none_1,test_sat_none_2,test_sat_thrust,test_sat_torque,test_sat_both,test_sat_none_2};
    std::vector<Satellite> satellite_vector_3 = {test_sat_none_2,test_sat_none_1,test_sat_thrust,test_sat_torque,test_sat_both};
    std::vector<Satellite> satellite_vector_4 = {test_sat_none_2,test_sat_none_1,test_sat_thrust,test_sat_torque,test_sat_both,test_sat_none_2};

    std::vector<std::vector<Satellite>> satellite_vector_vector = {satellite_vector_1,satellite_vector_2,satellite_vector_3,satellite_vector_4};
    std::vector<Satellite> single_satellite_vec_1 = {test_sat_both};
    std::vector<Satellite> single_satellite_vec_2 = {test_sat_torque};
    satellite_vector_vector.push_back(single_satellite_vec_1);
    satellite_vector_vector.push_back(single_satellite_vec_2);

    double timestep = 2;
    double total_sim_time = 300;

    // Drag parameters
    double F_10 = 100;  // Solar radio ten centimeter flux
    double A_p = 120;   // Geomagnetic A_p index

    const std::string plotting_term = "png"; // qt terminal opens window, probably not suited to running on a remote Github actions runner
    const std::string output_file_name = "test_plot";
    // Collect drag parameters into a tuple with F_10 first and A_p second
    std::pair<double, double> drag_elements = {F_10, A_p};
    for (std::vector<Satellite> satellite_vector : satellite_vector_vector){
        sim_and_draw_orbit_gnuplot(satellite_vector, timestep, total_sim_time,
            epsilon,false,false,drag_elements,plotting_term,output_file_name);
        sim_and_draw_orbit_gnuplot(satellite_vector, timestep, total_sim_time,
                    epsilon,true,false,drag_elements,plotting_term,output_file_name);
        sim_and_draw_orbit_gnuplot(satellite_vector, timestep, total_sim_time,
                        epsilon,true,true,drag_elements,plotting_term,output_file_name);
    }

  }