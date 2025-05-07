#include <gtest/gtest.h>

#include <iostream>

#include "Satellite.h"
#include "utils.h"

TEST(UtilsTests, OrbitalElementPlottingTests) {
  SimParameters sim_parameters("../tests/sim_parameters_baseline.json");
  // Going to have one satellite with two thrust profiles and two torque
  // profiles (added in different ways), one satellite with just two thrust
  // profiles, one with just two torque profiles and one with none
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

  std::array<double, 3> another_torque_vec = {-0.00001, 0.00002, 0.00003};
  test_sat_both.add_bodyframe_torque_profile(another_torque_vec, t_torque_start,
                                             t_torque_end);

  Satellite test_sat_thrust("../tests/elliptical_orbit_test_1.json");
  test_sat_thrust.add_LVLH_thrust_profile(
      LVLH_thrust_direction, thrust_magnitude, t_thrust_start, t_thrust_end);
  test_sat_thrust.add_LVLH_thrust_profile(LVLH_thrust_vec_2, t_thrust_2_start,
                                          t_thrust_2_end);

  Satellite test_sat_torque("../tests/elliptical_orbit_test_4.json");
  test_sat_torque.add_bodyframe_torque_profile(
      torque_direction, torque_magnitude, t_torque_start, t_torque_end);
  test_sat_torque.add_bodyframe_torque_profile(another_torque_vec,
                                               t_torque_start, t_torque_end);

  Satellite test_sat_none_1("../tests/elliptical_orbit_test_1.json");
  Satellite test_sat_none_2("../tests/elliptical_orbit_test_4.json");

  std::vector<Satellite> satellite_vector_1 = {test_sat_none_1, test_sat_none_2,
                                               test_sat_thrust, test_sat_torque,
                                               test_sat_both};
  std::vector<Satellite> satellite_vector_2 = {
      test_sat_none_1, test_sat_none_2, test_sat_thrust,
      test_sat_torque, test_sat_both,   test_sat_none_2};
  std::vector<Satellite> satellite_vector_3 = {test_sat_none_2, test_sat_none_1,
                                               test_sat_thrust, test_sat_torque,
                                               test_sat_both};
  std::vector<Satellite> satellite_vector_4 = {
      test_sat_none_2, test_sat_none_1, test_sat_thrust,
      test_sat_torque, test_sat_both,   test_sat_none_2};

  std::vector<std::vector<Satellite>> satellite_vector_vector = {
      satellite_vector_1, satellite_vector_2, satellite_vector_3,
      satellite_vector_4};
  std::vector<Satellite> single_satellite_vec_1 = {test_sat_both};
  std::vector<Satellite> single_satellite_vec_2 = {test_sat_torque};
  satellite_vector_vector.push_back(single_satellite_vec_1);
  satellite_vector_vector.push_back(single_satellite_vec_2);
  sim_parameters.initial_timestep_guess = 2;
  sim_parameters.total_sim_time = 300;

  // Drag parameters
  sim_parameters.F_10 = 100;  // Solar radio ten centimeter flux
  sim_parameters.A_p = 120;   // Geomagnetic A_p index

  // Collect drag parameters into a tuple with F_10 first and A_p second
  std::string file_name = "test_plot";
  std::vector<std::string> orbital_elements = {
      "Semimajor Axis",        "Eccentricity",
      "Inclination",           "RAAN",
      "Argument of Periapsis", "True Anomaly",
      "Orbital Rate",          "Orbital Angular Acceleration",
      "Total Energy"};
  for (std::vector<Satellite> satellite_vector : satellite_vector_vector) {
    for (std::string element_name : orbital_elements) {
      sim_parameters.perturbation_bool = false;
      sim_parameters.drag_bool = false;
      sim_and_plot_orbital_elem_gnuplot(satellite_vector, sim_parameters,
                                        element_name, file_name);
      sim_parameters.perturbation_bool = true;
      sim_parameters.drag_bool = false;
      sim_and_plot_orbital_elem_gnuplot(satellite_vector, sim_parameters,
                                        element_name, file_name);
      sim_parameters.perturbation_bool = true;
      sim_parameters.drag_bool = true;
      sim_and_plot_orbital_elem_gnuplot(satellite_vector, sim_parameters,
                                        element_name, file_name);
    }
  }
}

// Now a similar test for attitude plotting
TEST(UtilsTests, AttitudeElementPlottingTests) {
  // Going to have one satellite with two thrust profiles and two torque
  // profiles (added in different ways), one satellite with just two thrust
  // profiles, one with just two torque profiles and one with none
  SimParameters sim_parameters("../tests/sim_parameters_baseline.json");
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

  std::array<double, 3> another_torque_vec = {-0.00001, 0.00002, 0.00003};
  test_sat_both.add_bodyframe_torque_profile(another_torque_vec, 0, 2);

  Satellite test_sat_thrust("../tests/elliptical_orbit_test_1.json");
  test_sat_thrust.add_LVLH_thrust_profile(
      LVLH_thrust_direction, thrust_magnitude, t_thrust_start, t_thrust_end);
  test_sat_thrust.add_LVLH_thrust_profile(LVLH_thrust_vec_2, t_thrust_2_start,
                                          t_thrust_2_end);

  Satellite test_sat_torque("../tests/elliptical_orbit_test_4.json");
  test_sat_torque.add_bodyframe_torque_profile(
      torque_direction, torque_magnitude, t_torque_start, t_torque_end);
  test_sat_torque.add_bodyframe_torque_profile(another_torque_vec,
                                               t_torque_start, t_torque_end);

  Satellite test_sat_none_1("../tests/elliptical_orbit_test_1.json");
  Satellite test_sat_none_2("../tests/elliptical_orbit_test_4.json");

  std::vector<Satellite> satellite_vector_1 = {test_sat_none_1, test_sat_none_2,
                                               test_sat_thrust, test_sat_torque,
                                               test_sat_both};
  std::vector<Satellite> satellite_vector_2 = {
      test_sat_none_1, test_sat_none_2, test_sat_thrust,
      test_sat_torque, test_sat_both,   test_sat_none_2};
  std::vector<Satellite> satellite_vector_3 = {test_sat_none_2, test_sat_none_1,
                                               test_sat_thrust, test_sat_torque,
                                               test_sat_both};
  std::vector<Satellite> satellite_vector_4 = {
      test_sat_none_2, test_sat_none_1, test_sat_thrust,
      test_sat_torque, test_sat_both,   test_sat_none_2};

  std::vector<std::vector<Satellite>> satellite_vector_vector = {
      satellite_vector_1, satellite_vector_2, satellite_vector_3,
      satellite_vector_4};
  std::vector<Satellite> single_satellite_vec_1 = {test_sat_both};
  std::vector<Satellite> single_satellite_vec_2 = {test_sat_torque};
  satellite_vector_vector.push_back(single_satellite_vec_1);
  satellite_vector_vector.push_back(single_satellite_vec_2);

  sim_parameters.initial_timestep_guess = 2;
  sim_parameters.total_sim_time = 300;

  // Drag parameters
  sim_parameters.F_10 = 100;  // Solar radio ten centimeter flux
  sim_parameters.A_p = 120;   // Geomagnetic A_p index

  std::string file_name = "test_plot";
  std::vector<std::string> attitude_elements = {
      "Roll",    "Pitch", "Yaw", "omega_x", "omega_y",
      "omega_z", "q_0",   "q_1", "q_2",     "q_3"};
  for (std::vector<Satellite> satellite_vector : satellite_vector_vector) {
    for (std::string element_name : attitude_elements) {
      sim_parameters.perturbation_bool = false;
      sim_parameters.drag_bool = false;
      sim_and_plot_attitude_evolution_gnuplot(satellite_vector, sim_parameters,
                                              element_name, file_name);
      sim_parameters.perturbation_bool = true;
      sim_parameters.drag_bool = false;
      sim_and_plot_attitude_evolution_gnuplot(satellite_vector, sim_parameters,
                                              element_name, file_name);
      sim_parameters.perturbation_bool = true;
      sim_parameters.drag_bool = true;
      sim_and_plot_attitude_evolution_gnuplot(satellite_vector, sim_parameters,
                                              element_name, file_name);
    }
  }
}

// Now one for the 3D orbital plot
TEST(UtilsTests, ThreeDimensionalOrbitPlotTest) {
  // Going to have one satellite with two thrust profiles and two torque
  // profiles (added in different ways), one satellite with just two thrust
  // profiles, one with just two torque profiles and one with none
  SimParameters sim_parameters("../tests/sim_parameters_baseline.json");
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

  std::array<double, 3> another_torque_vec = {-0.00001, 0.00002, 0.00003};
  test_sat_both.add_bodyframe_torque_profile(another_torque_vec, t_torque_start,
                                             t_torque_end);

  Satellite test_sat_thrust("../tests/elliptical_orbit_test_1.json");
  test_sat_thrust.add_LVLH_thrust_profile(
      LVLH_thrust_direction, thrust_magnitude, t_thrust_start, t_thrust_end);
  test_sat_thrust.add_LVLH_thrust_profile(LVLH_thrust_vec_2, t_thrust_2_start,
                                          t_thrust_2_end);

  Satellite test_sat_torque("../tests/elliptical_orbit_test_4.json");
  test_sat_torque.add_bodyframe_torque_profile(
      torque_direction, torque_magnitude, t_torque_start, t_torque_end);
  test_sat_torque.add_bodyframe_torque_profile(another_torque_vec,
                                               t_torque_start, t_torque_end);

  Satellite test_sat_none_1("../tests/elliptical_orbit_test_1.json");
  Satellite test_sat_none_2("../tests/elliptical_orbit_test_4.json");

  std::vector<Satellite> satellite_vector_1 = {test_sat_none_1, test_sat_none_2,
                                               test_sat_thrust, test_sat_torque,
                                               test_sat_both};
  std::vector<Satellite> satellite_vector_2 = {
      test_sat_none_1, test_sat_none_2, test_sat_thrust,
      test_sat_torque, test_sat_both,   test_sat_none_2};
  std::vector<Satellite> satellite_vector_3 = {test_sat_none_2, test_sat_none_1,
                                               test_sat_thrust, test_sat_torque,
                                               test_sat_both};
  std::vector<Satellite> satellite_vector_4 = {
      test_sat_none_2, test_sat_none_1, test_sat_thrust,
      test_sat_torque, test_sat_both,   test_sat_none_2};

  std::vector<std::vector<Satellite>> satellite_vector_vector = {
      satellite_vector_1, satellite_vector_2, satellite_vector_3,
      satellite_vector_4};
  std::vector<Satellite> single_satellite_vec_1 = {test_sat_both};
  std::vector<Satellite> single_satellite_vec_2 = {test_sat_torque};
  satellite_vector_vector.push_back(single_satellite_vec_1);
  satellite_vector_vector.push_back(single_satellite_vec_2);

  sim_parameters.initial_timestep_guess = 2;
  sim_parameters.total_sim_time = 300;

  // Drag parameters
  sim_parameters.F_10 = 100;  // Solar radio ten centimeter flux
  sim_parameters.A_p = 120;   // Geomagnetic A_p index

  sim_parameters.terminal_name_3D =
      "png";  // qt terminal opens window, probably not suited to running on a
              // remote Github actions runner
  const std::string output_file_name = "test_plot";

  for (std::vector<Satellite> satellite_vector : satellite_vector_vector) {
    sim_parameters.perturbation_bool = false;
    sim_parameters.drag_bool = false;
    sim_and_draw_orbit_gnuplot(satellite_vector, sim_parameters,
                               output_file_name);
    sim_parameters.perturbation_bool = true;
    sim_parameters.drag_bool = false;
    sim_and_draw_orbit_gnuplot(satellite_vector, sim_parameters,
                               output_file_name);
    sim_parameters.perturbation_bool = true;
    sim_parameters.drag_bool = true;
    sim_and_draw_orbit_gnuplot(satellite_vector, sim_parameters,
                               output_file_name);
  }
}

TEST(UtilsTests, GroundStationConnectivityDistancePlotTests) {
  // Going to have one satellite with two thrust profiles and two torque
  // profiles (added in different ways), one satellite with just two thrust
  // profiles, one with just two torque profiles and one with none
  SimParameters sim_parameters("../tests/sim_parameters_baseline.json");
  Satellite test_sat_both("../tests/input_2.json");
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

  std::array<double, 3> another_torque_vec = {-0.00001, 0.00002, 0.00003};
  test_sat_both.add_bodyframe_torque_profile(another_torque_vec, t_torque_start,
                                             t_torque_end);

  Satellite test_sat_thrust("../tests/input_2.json");
  test_sat_thrust.add_LVLH_thrust_profile(
      LVLH_thrust_direction, thrust_magnitude, t_thrust_start, t_thrust_end);
  test_sat_thrust.add_LVLH_thrust_profile(LVLH_thrust_vec_2, t_thrust_2_start,
                                          t_thrust_2_end);

  Satellite test_sat_torque("../tests/input_3.json");
  test_sat_torque.add_bodyframe_torque_profile(
      torque_direction, torque_magnitude, t_torque_start, t_torque_end);
  test_sat_torque.add_bodyframe_torque_profile(another_torque_vec,
                                               t_torque_start, t_torque_end);

  Satellite test_sat_none_1("../tests/input_2.json");
  Satellite test_sat_none_2("../tests/input_3.json");

  std::vector<Satellite> satellite_vector_1 = {test_sat_none_1, test_sat_none_2,
                                               test_sat_thrust, test_sat_torque,
                                               test_sat_both};
  std::vector<Satellite> satellite_vector_2 = {
      test_sat_none_1, test_sat_none_2, test_sat_thrust,
      test_sat_torque, test_sat_both,   test_sat_none_2};
  std::vector<Satellite> satellite_vector_3 = {test_sat_none_2, test_sat_none_1,
                                               test_sat_thrust, test_sat_torque,
                                               test_sat_both};
  std::vector<Satellite> satellite_vector_4 = {
      test_sat_none_2, test_sat_none_1, test_sat_thrust,
      test_sat_torque, test_sat_both,   test_sat_none_2};

  std::vector<std::vector<Satellite>> satellite_vector_vector = {
      satellite_vector_1, satellite_vector_2, satellite_vector_3,
      satellite_vector_4};
  std::vector<Satellite> single_satellite_vec_1 = {test_sat_both};
  std::vector<Satellite> single_satellite_vec_2 = {test_sat_torque};
  satellite_vector_vector.push_back(single_satellite_vec_1);
  satellite_vector_vector.push_back(single_satellite_vec_2);
  sim_parameters.initial_timestep_guess = 2;
  sim_parameters.total_sim_time = 300;

  // Drag parameters
  sim_parameters.F_10 = 100;  // Solar radio ten centimeter flux
  sim_parameters.A_p = 120;   // Geomagnetic A_p index

  double gs_latitude = 2;
  double gs_longitude = 5;
  double gs_altitude = 10;
  double max_beam_angle_from_normal = 65;  // deg
  int num_beams = 2;
  sim_parameters.total_sim_time = 25000;
  PhasedArrayGroundStation example_ground_station_1(
      gs_latitude, gs_longitude, gs_altitude, max_beam_angle_from_normal,
      num_beams);

  std::string file_name = "test_plot";
  std::vector<std::string> orbital_elements = {
      "Semimajor Axis",        "Eccentricity",
      "Inclination",           "RAAN",
      "Argument of Periapsis", "True Anomaly",
      "Orbital Rate",          "Orbital Angular Acceleration",
      "Total Energy"};
  for (std::vector<Satellite> satellite_vector : satellite_vector_vector) {
    for (std::string element_name : orbital_elements) {
      sim_parameters.perturbation_bool = false;
      sim_parameters.drag_bool = false;
      sim_and_plot_gs_connectivity_distance_gnuplot(example_ground_station_1,
                                                    satellite_vector,
                                                    sim_parameters, file_name);
      sim_parameters.perturbation_bool = true;
      sim_parameters.drag_bool = false;
      sim_and_plot_gs_connectivity_distance_gnuplot(example_ground_station_1,
                                                    satellite_vector,
                                                    sim_parameters, file_name);
      sim_parameters.perturbation_bool = true;
      sim_parameters.drag_bool = true;
      sim_and_plot_gs_connectivity_distance_gnuplot(example_ground_station_1,
                                                    satellite_vector,
                                                    sim_parameters, file_name);
    }
  }
}

TEST(UtilsTests, GroundStationConnectivityPlotTests) {
  // Going to have one satellite with two thrust profiles and two torque
  // profiles (added in different ways), one satellite with just two thrust
  // profiles, one with just two torque profiles and one with none
  SimParameters sim_parameters("../tests/sim_parameters_baseline.json");
  Satellite test_sat_both("../tests/input_2.json");
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

  std::array<double, 3> another_torque_vec = {-0.00001, 0.00002, 0.00003};
  test_sat_both.add_bodyframe_torque_profile(another_torque_vec, t_torque_start,
                                             t_torque_end);

  Satellite test_sat_thrust("../tests/input_2.json");
  test_sat_thrust.add_LVLH_thrust_profile(
      LVLH_thrust_direction, thrust_magnitude, t_thrust_start, t_thrust_end);
  test_sat_thrust.add_LVLH_thrust_profile(LVLH_thrust_vec_2, t_thrust_2_start,
                                          t_thrust_2_end);

  Satellite test_sat_torque("../tests/input_3.json");
  test_sat_torque.add_bodyframe_torque_profile(
      torque_direction, torque_magnitude, t_torque_start, t_torque_end);
  test_sat_torque.add_bodyframe_torque_profile(another_torque_vec,
                                               t_torque_start, t_torque_end);

  Satellite test_sat_none_1("../tests/input_2.json");
  Satellite test_sat_none_2("../tests/input_3.json");

  std::vector<Satellite> satellite_vector_1 = {test_sat_none_1, test_sat_none_2,
                                               test_sat_thrust, test_sat_torque,
                                               test_sat_both};
  std::vector<Satellite> satellite_vector_2 = {
      test_sat_none_1, test_sat_none_2, test_sat_thrust,
      test_sat_torque, test_sat_both,   test_sat_none_2};
  std::vector<Satellite> satellite_vector_3 = {test_sat_none_2, test_sat_none_1,
                                               test_sat_thrust, test_sat_torque,
                                               test_sat_both};
  std::vector<Satellite> satellite_vector_4 = {
      test_sat_none_2, test_sat_none_1, test_sat_thrust,
      test_sat_torque, test_sat_both,   test_sat_none_2};

  std::vector<std::vector<Satellite>> satellite_vector_vector = {
      satellite_vector_1, satellite_vector_2, satellite_vector_3,
      satellite_vector_4};
  std::vector<Satellite> single_satellite_vec_1 = {test_sat_both};
  std::vector<Satellite> single_satellite_vec_2 = {test_sat_torque};
  satellite_vector_vector.push_back(single_satellite_vec_1);
  satellite_vector_vector.push_back(single_satellite_vec_2);
  sim_parameters.initial_timestep_guess = 2;
  sim_parameters.total_sim_time = 300;

  // Drag parameters
  sim_parameters.F_10 = 100;  // Solar radio ten centimeter flux
  sim_parameters.A_p = 120;   // Geomagnetic A_p index

  double gs_latitude = 2;
  double gs_longitude = 5;
  double gs_altitude = 10;
  double max_beam_angle_from_normal = 65;  // deg
  int num_beams = 2;
  sim_parameters.total_sim_time = 25000;
  PhasedArrayGroundStation example_ground_station_1(
      gs_latitude, gs_longitude, gs_altitude, max_beam_angle_from_normal,
      num_beams);

  std::string file_name = "test_plot";
  std::vector<std::string> orbital_elements = {
      "Semimajor Axis",        "Eccentricity",
      "Inclination",           "RAAN",
      "Argument of Periapsis", "True Anomaly",
      "Orbital Rate",          "Orbital Angular Acceleration",
      "Total Energy"};
  for (std::vector<Satellite> satellite_vector : satellite_vector_vector) {
    for (std::string element_name : orbital_elements) {
      sim_parameters.perturbation_bool = false;
      sim_parameters.drag_bool = false;
      sim_and_plot_gs_connectivity_gnuplot(example_ground_station_1,
                                           satellite_vector, sim_parameters,
                                           file_name);
      sim_parameters.perturbation_bool = true;
      sim_parameters.drag_bool = false;
      sim_and_plot_gs_connectivity_gnuplot(example_ground_station_1,
                                           satellite_vector, sim_parameters,
                                           file_name);
      sim_parameters.perturbation_bool = true;
      sim_parameters.drag_bool = true;
      sim_and_plot_gs_connectivity_gnuplot(example_ground_station_1,
                                           satellite_vector, sim_parameters,
                                           file_name);
    }
  }
}

TEST(UtilsTests, LowThrustTransferTest1) {
  Satellite test_circular_sat("../tests/circular_orbit_test_2_input.json");
  double final_orbit_semimajor_axis = 12500;  // km
  double thrust_magnitude = 0.1;              // N
  double transfer_initiation_time = 10;       // s
  int error_code = add_lowthrust_orbit_transfer(
      test_circular_sat, final_orbit_semimajor_axis, thrust_magnitude,
      transfer_initiation_time);

  EXPECT_TRUE(error_code == 0)
      << "Non-circular orbit flag got thrown incorrectly\n";
}

TEST(UtilsTests, LowThrustTransferTest2) {
  Satellite test_circular_sat("../tests/elliptical_orbit_test_1.json");
  double final_orbit_semimajor_axis = 50000;  // km
  double thrust_magnitude = 0.1;              // N
  double transfer_initiation_time = 10;       // s
  int error_code = add_lowthrust_orbit_transfer(
      test_circular_sat, final_orbit_semimajor_axis, thrust_magnitude,
      transfer_initiation_time);

  EXPECT_TRUE(error_code == 1)
      << "Non-circular orbit should have been thrown here\n";
}