#include <iostream>

#include "Satellite.h"
#include "utils.h"

int main() {
  // Initialize satellite object from an input JSON file
  Satellite test_sat_1("../input.json");
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

  Satellite test_sat_2("../input_2.json");

  Satellite test_sat_3("../input_3.json");

  std::vector<Satellite> satellite_vector_1 = {test_sat_1, test_sat_2,
                                               test_sat_3};

  double timestep = 2;
  double total_sim_time = 25000;
  double epsilon = pow(10, -12);
  // sim_and_draw_orbit_gnuplot(satellite_vector_1, timestep, total_sim_time,
  //                            epsilon);

  // Now some demonstrations of plotting orbital parameters
  Satellite test_sat_4("../input_4.json");
  Satellite test_sat_5("../input_5.json");
  Satellite test_sat_6("../input_6.json");

  std::vector<Satellite> satellite_vector_2 = {test_sat_4, test_sat_5,
                                               test_sat_6};
  total_sim_time = 9952;
  std::string file_name = "Arg of Periapsis Plot";
  sim_and_plot_orbital_elem_gnuplot(satellite_vector_2, timestep,
                                    total_sim_time, epsilon,
                                    "Argument of Periapsis", file_name, true);
  Satellite test_sat_7("../input_7.json");
  std::array<double, 3> torque_direction = {0, -1, 0};
  double torque_magnitude = 0.0005;  // N
  double t_torque_start = 3000;
  double t_torque_end = 3002;
  test_sat_7.add_bodyframe_torque_profile(torque_direction, torque_magnitude,
                                          t_torque_start, t_torque_end);
  std::vector<Satellite> satellite_vector_3 = {test_sat_7};
  file_name = "Pitch Plot";
  sim_and_plot_attitude_evolution_gnuplot(
      satellite_vector_3, timestep, total_sim_time, epsilon, "Pitch", file_name, false);

  // Now let's demonstrate effect of atmospheric drag approximation
  Satellite test_sat_8("../input_8.json");
  Satellite test_sat_9("../input_9.json");
  std::vector<Satellite> satellite_vector_4 = {test_sat_8, test_sat_9};

  // Drag parameters
  double F_10 = 100;  // Solar radio ten centimeter flux
  double A_p = 120;   // Geomagnetic A_p index

  // Collect drag parameters into a tuple with F_10 first and A_p second
  std::pair<double, double> drag_elements = {F_10, A_p};
  total_sim_time = 10000;
  epsilon = pow(10,-14);
  file_name = "Eccentricity Plot";
  sim_and_plot_orbital_elem_gnuplot(satellite_vector_4, timestep,
                                    total_sim_time, epsilon, "Eccentricity",
                                    file_name, false, true, drag_elements);
  file_name = "Semimajor Axis";
  sim_and_plot_orbital_elem_gnuplot(satellite_vector_4, timestep,
                                    total_sim_time, epsilon, "Semimajor Axis",
                                    file_name, false, true, drag_elements);
  return 0;
}