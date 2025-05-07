#include <iostream>

#include "Satellite.h"
#include "utils.h"
#include "PhasedArrayGroundStation.h"

int main() {
  // This file demonstrates a few different ways you can run and
  // visualize data from satellite simulations.
  // First let's initialize the struct containing parameters for simulation and plotting
  SimParameters sim_parameters("../sim_parameters_example.json");
  // Initialize satellite object from an input JSON file
  Satellite test_sat_1("../example_satellite_input_files/input.json");
  // Define parameters for an LVLH frame thrust profile
  std::array<double, 3> LVLH_thrust_direction = {1, 0, 0};
  double thrust_magnitude = 100;  // N
  double t_thrust_start = 5000;
  double t_thrust_end = 6000;
  // Add the thrust profile to the satellite object
  test_sat_1.add_LVLH_thrust_profile(LVLH_thrust_direction, thrust_magnitude,
                                     t_thrust_start, t_thrust_end);

  std::array<double, 3> LVLH_thrust_vec_2 = {0.51, 20, -5};
  double t_thrust_2_start = 5500;
  double t_thrust_2_end = 6500;
  test_sat_1.add_LVLH_thrust_profile(LVLH_thrust_vec_2, t_thrust_2_start,
                                     t_thrust_2_end);

  Satellite test_sat_2("../example_satellite_input_files/input_2.json");

  Satellite test_sat_3("../example_satellite_input_files/input_3.json");

  std::vector<Satellite> satellite_vector_1 = {test_sat_1, test_sat_2,
                                               test_sat_3};

  sim_parameters.initial_timestep_guess = 2;
  sim_parameters.total_sim_time = 25000;
  sim_parameters.epsilon = pow(10, -12);
  sim_and_draw_orbit_gnuplot(satellite_vector_1, sim_parameters);

  // Now some demonstrations of plotting orbital parameters
  Satellite test_sat_4("../example_satellite_input_files/input_4.json");
  Satellite test_sat_5("../example_satellite_input_files/input_5.json");
  Satellite test_sat_6("../example_satellite_input_files/input_6.json");

  std::vector<Satellite> satellite_vector_2 = {test_sat_4, test_sat_5,
                                               test_sat_6};
  sim_parameters.total_sim_time = 9952;
  std::string file_name = "Arg of Periapsis Plot";
  sim_parameters.perturbation_bool = false;
  sim_and_plot_orbital_elem_gnuplot(satellite_vector_2, sim_parameters,
                                    "Argument of Periapsis", file_name);
  Satellite test_sat_7("../example_satellite_input_files/input_7.json");
  std::array<double, 3> torque_direction = {0, -1, 0};
  double torque_magnitude = 0.0005;  // N
  double t_torque_start = 3000;
  double t_torque_end = 3002;
  test_sat_7.add_bodyframe_torque_profile(torque_direction, torque_magnitude,
                                          t_torque_start, t_torque_end);
  std::vector<Satellite> satellite_vector_3 = {test_sat_7};
  file_name = "Pitch Plot";
  sim_and_plot_attitude_evolution_gnuplot(
      satellite_vector_3,sim_parameters, "Pitch", file_name);

  // Now let's demonstrate effect of atmospheric drag approximation
  Satellite test_sat_8("../example_satellite_input_files/input_8.json");
  Satellite test_sat_9("../example_satellite_input_files/input_9.json");
  std::vector<Satellite> satellite_vector_4 = {test_sat_8, test_sat_9};

  // Drag parameters
  double F_10 = 100;  // Solar radio ten centimeter flux
  double A_p = 120;   // Geomagnetic A_p index

  // Collect drag parameters into a tuple with F_10 first and A_p second
  std::pair<double, double> drag_elements = {F_10, A_p};
  sim_parameters.total_sim_time = 10000;
  sim_parameters.epsilon = pow(10,-14);
  sim_parameters.drag_bool = true;
  file_name = "Eccentricity Plot";
  sim_and_plot_orbital_elem_gnuplot(satellite_vector_4,sim_parameters,
                                 "Eccentricity", file_name);
  file_name = "Semimajor Axis";
  sim_and_plot_orbital_elem_gnuplot(satellite_vector_4, sim_parameters, 
                                  "Semimajor Axis",file_name);


  // Now demonstrating phased array ground station coverage
  Satellite test_sat_10("../example_satellite_input_files/input_10.json");
  Satellite test_sat_11("../example_satellite_input_files/input_11.json");
  Satellite test_sat_12("../example_satellite_input_files/input_12.json");
  Satellite test_sat_13("../example_satellite_input_files/input_13.json");
  double gs_latitude = 2;
  double gs_longitude = 5;
  double gs_altitude = 10;
  double max_beam_angle_from_normal = 65; //deg
  int num_beams = 5;
  sim_parameters.total_sim_time = 40000;
  sim_parameters.perturbation_bool = true;
  PhasedArrayGroundStation example_ground_station_1(gs_latitude,gs_longitude,
    gs_altitude,max_beam_angle_from_normal,num_beams);
  std::vector<Satellite> long_sat_vec = {test_sat_1,test_sat_2,test_sat_3,test_sat_10,test_sat_11,test_sat_12,test_sat_13};
  std::string coverage_file_name = "Ground station connectivity with " + std::to_string(num_beams) + " beams";
  sim_and_plot_gs_connectivity_gnuplot(example_ground_station_1,
    long_sat_vec, sim_parameters,
    coverage_file_name);
  
  num_beams = 2;
  PhasedArrayGroundStation example_ground_station_2(gs_latitude,gs_longitude,
    gs_altitude,max_beam_angle_from_normal,num_beams);
  coverage_file_name = "Ground station connectivity with " + std::to_string(num_beams) + " beams";
  sim_and_plot_gs_connectivity_gnuplot(example_ground_station_2,
    long_sat_vec, sim_parameters,
    coverage_file_name);

  // Now demonstrating low-thrust orbital transfer
  Satellite orbit_raising_demo_sat("../example_satellite_input_files/circular_orbit_raising_test_input.json");
  Satellite orbit_lowering_demo_sat("../example_satellite_input_files/circular_orbit_lowering_test_input.json");
  double raising_final_orbit_semimajor_axis = 20000; // km
  thrust_magnitude = 1; // N
  double transfer_initiation_time = 0; // s
  int error_code_raising = add_lowthrust_orbit_transfer(orbit_raising_demo_sat, raising_final_orbit_semimajor_axis,
      thrust_magnitude, transfer_initiation_time);

  double lowering_final_orbit_semimajor_axis = 10000; // km
  int error_code_lowering = add_lowthrust_orbit_transfer(orbit_lowering_demo_sat, lowering_final_orbit_semimajor_axis,
    thrust_magnitude, transfer_initiation_time);
  std::vector<Satellite> orbit_transfer_demo_vec = {orbit_raising_demo_sat,orbit_lowering_demo_sat};
  sim_parameters.initial_timestep_guess = 2;
  sim_parameters.total_sim_time = 400000;
  sim_parameters.epsilon = pow(10, -12);
  sim_parameters.perturbation_bool = false;
  sim_parameters.drag_bool = false;
  // Axis tick increments to make axes less cluttered
  sim_parameters.x_increment = pow(10,7);
  sim_parameters.y_increment = pow(10,7);
  sim_parameters.z_increment = 5*pow(10,6);
  // file_name = "Semimajor axis transfer";
  // sim_and_plot_orbital_elem_gnuplot(orbit_transfer_demo_vec, sim_parameters, "Semimajor Axis", file_name);
  sim_and_draw_orbit_gnuplot(orbit_transfer_demo_vec,sim_parameters);


  // Now let's demonstrate changing the argument of periapsis
  // Calibration strategy: there are some inherent oscillations in, e.g., arg of periapsis, even in absence of external forces or perturbations
  // This oscillation magnitude can be different in initial and final orbits, so current strategy is to find the mean val of these oscillations
  // at initial and final orbits, use those to get an offset, which can be applied to the final argument of periapsis to try to make the
  // oscillations at the final orbit be close to the target value

  Satellite arg_periapsis_change_sat("../example_satellite_input_files/arg_of_periapsis_test_input.json");
  sim_parameters.epsilon = pow(10,-14);
  sim_parameters.total_sim_time = 400000; 
  sim_parameters.perturbation_bool = false;
  sim_parameters.drag_bool = false;

  t_thrust_start = 50000;
  double final_arg_of_periapsis_deg = 25;

  thrust_magnitude = 0.1; // N
  // std::cout << "orbital period: " << arg_periapsis_change_sat.calculate_orbital_period() << "\n";
  arg_periapsis_change_sat.add_maneuver("Argument of Periapsis Change",t_thrust_start,final_arg_of_periapsis_deg,thrust_magnitude);
  file_name = "Arg of periapsis transfer";
  std::vector<Satellite> arg_of_periapsis_transfer_vec = {arg_periapsis_change_sat};
  sim_and_plot_orbital_elem_gnuplot(arg_of_periapsis_transfer_vec, sim_parameters, "Argument of Periapsis", file_name);

  return 0;
}