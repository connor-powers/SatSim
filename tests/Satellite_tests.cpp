#include <gtest/gtest.h>

#include <iostream>

#include "Satellite.h"
#include "utils.h"

const double tolerance = pow(10.0, -7);
// Setting a different tolerance for semimajor axis and orbital radius than the other orbital
// parameters since there appears to be a minimum error associated with
// converting position and velocity to semimajor axis, best guess is this has to
// do with the scale of distances and/or velocities being dealt with here
const double length_tolerance = pow(10.0, -6);
const double epsilon = pow(10.0, -12);
const double energy_cons_relative_tolerance = pow(10.0, -5);

// Elliptical orbit tests

TEST(EllipticalOrbitTests, EvolvedOrbitalSpeed1) {
  // Starting at true anomaly=0 means it's starting at perigee, which is where
  // its orbital speed should be maximum
  Satellite test_satellite("../tests/elliptical_orbit_test_1.json");
  double calculated_initial_speed = test_satellite.get_speed();
  double test_timestep = 0.1;  // s
  double sim_time = 1;         // s
  double current_time = test_satellite.get_instantaneous_time();
  double next_timestep = 0;
  while (current_time < sim_time) {
    std::pair<double, int> new_timestep_and_error_code =
        test_satellite.evolve_RK45(epsilon, test_timestep);
    next_timestep = new_timestep_and_error_code.first;
    int error_code = new_timestep_and_error_code.second;
    test_timestep = next_timestep;
    current_time = test_satellite.get_instantaneous_time();
  }
  double calculated_evolved_speed = test_satellite.get_speed();

  EXPECT_TRUE(calculated_initial_speed > calculated_evolved_speed)
      << "Perigee speed not larger than calculated evolved speed. Difference: "
      << calculated_initial_speed - calculated_evolved_speed << "\n";
}

TEST(EllipticalOrbitTests, EvolvedOrbitalSpeed2) {
  // Starting at true anomaly=180 means it's starting at apogee, which is where
  // its orbital speed should be minimum
  Satellite test_satellite("../tests/elliptical_orbit_test_2.json");
  double calculated_initial_speed = test_satellite.get_speed();
  double test_timestep = 1;  // s
  std::pair<double, int> new_timestep_and_error_code =
      test_satellite.evolve_RK45(epsilon, test_timestep);
  double next_timestep = new_timestep_and_error_code.first;
  int error_code = new_timestep_and_error_code.second;
  double calculated_evolved_speed = test_satellite.get_speed();

  EXPECT_TRUE(calculated_initial_speed < calculated_evolved_speed)
      << "Apogee speed not smaller than calculated evolved speed. Difference: "
      << calculated_initial_speed - calculated_evolved_speed << "\n";
}

TEST(EllipticalOrbitTests, ConstantEvolvedOrbitalElementsTest) {
  // The idea behind this test is that after evolving a timestep, orbital
  // elements besides true anomaly should be constant (when J2 perturbation is
  // not taken into account)

  Satellite test_satellite("../tests/elliptical_orbit_test_2.json");
  std::array<double, 6> initial_orbit_elements =
      test_satellite.get_orbital_elements();

  double test_timestep = 1;  // s
  bool perturbation_bool = false;
  std::pair<double, int> new_timestep_and_error_code =
      test_satellite.evolve_RK45(epsilon, test_timestep, perturbation_bool);
  double next_timestep = new_timestep_and_error_code.first;
  int error_code = new_timestep_and_error_code.second;
  std::array<double, 6> evolved_orbit_elements =
      test_satellite.get_orbital_elements();

  std::array<std::string, 6> orbital_element_name_array;
  orbital_element_name_array.at(0) = "Semimajor Axis";
  orbital_element_name_array.at(1) = "Eccentricity";
  orbital_element_name_array.at(2) = "Inclination";
  orbital_element_name_array.at(3) = "RAAN";
  orbital_element_name_array.at(4) = "Argument of Periapsis";
  orbital_element_name_array.at(5) = "True Anomaly";

  for (size_t orbital_elem_index = 0; orbital_elem_index < 5;
       orbital_elem_index++) {
    // True anomaly shouldn't be constant over evolution
    if (orbital_elem_index == 0) {
      EXPECT_TRUE(abs(initial_orbit_elements.at(orbital_elem_index) -
                      evolved_orbit_elements.at(orbital_elem_index)) <
                      length_tolerance)
          << orbital_element_name_array.at(orbital_elem_index)
          << " was not constant within tolerance. Diff:"
          << initial_orbit_elements.at(orbital_elem_index) -
                 evolved_orbit_elements.at(orbital_elem_index)
          << "\n";
    } else {
      EXPECT_TRUE(abs(initial_orbit_elements.at(orbital_elem_index) -
                      evolved_orbit_elements.at(orbital_elem_index)) <
                  tolerance)
          << orbital_element_name_array.at(orbital_elem_index)
          << " was not constant within tolerance. Diff:"
          << initial_orbit_elements.at(orbital_elem_index) -
                 evolved_orbit_elements.at(orbital_elem_index)
          << "\n";
    }
  }
}

TEST(EllipticalOrbitTests, BasicOrbitalElementsTest) {
  Satellite test_satellite("../tests/elliptical_orbit_test_2.json");
  std::array<double, 6> initial_orbit_elements =
      test_satellite.get_orbital_elements();

  test_satellite.update_orbital_elements_from_position_and_velocity();
  std::array<double, 6> recalculated_orbit_elements =
      test_satellite.get_orbital_elements();

  std::array<std::string, 6> orbital_element_name_array;
  orbital_element_name_array.at(0) = "Semimajor Axis";
  orbital_element_name_array.at(1) = "Eccentricity";
  orbital_element_name_array.at(2) = "Inclination";
  orbital_element_name_array.at(3) = "RAAN";
  orbital_element_name_array.at(4) = "Argument of Periapsis";
  orbital_element_name_array.at(5) = "True Anomaly";

  for (size_t orbital_elem_index = 0; orbital_elem_index < 6;
       orbital_elem_index++) {
    if (orbital_elem_index == 0) {
      EXPECT_TRUE(abs(initial_orbit_elements.at(orbital_elem_index) -
                      recalculated_orbit_elements.at(orbital_elem_index)) <
                      length_tolerance)
          << orbital_element_name_array.at(orbital_elem_index)
          << " was not constant within tolerance. Diff:"
          << initial_orbit_elements.at(orbital_elem_index) -
                 recalculated_orbit_elements.at(orbital_elem_index)
          << "\n";
    } else {
      EXPECT_TRUE(abs(initial_orbit_elements.at(orbital_elem_index) -
                      recalculated_orbit_elements.at(orbital_elem_index)) <
                  tolerance)
          << orbital_element_name_array.at(orbital_elem_index)
          << " was not constant within tolerance. Diff:"
          << initial_orbit_elements.at(orbital_elem_index) -
                 recalculated_orbit_elements.at(orbital_elem_index)
          << "\n";
    }
  }
}

TEST(EllipticalOrbitTests,OrbitalRadiusCalcs1) {
  Satellite test_satellite("../tests/elliptical_orbit_test_2.json");
  double orbital_radius_perifocal=test_satellite.get_radius();
  double orbital_radius_ECI=test_satellite.get_radius_ECI();
  EXPECT_TRUE(abs(orbital_radius_perifocal - orbital_radius_ECI) < length_tolerance)
      << "Difference between orbital radii calculated with "
      " perifocal and ECI coordinates: "
      << orbital_radius_perifocal - orbital_radius_ECI <<"\n";
}

TEST(EllipticalOrbitTests,OrbitalRadiusCalcs2) {
  Satellite test_satellite("../tests/elliptical_orbit_test_2.json");
  double test_timestep = 0.1;  // s
  bool perturbation_bool = false;
  std::pair<double, int> new_timestep_and_error_code =
      test_satellite.evolve_RK45(epsilon, test_timestep, perturbation_bool);
  double next_timestep = new_timestep_and_error_code.first;
  int error_code = new_timestep_and_error_code.second;
  double orbital_radius_perifocal=test_satellite.get_radius();
  double orbital_radius_ECI=test_satellite.get_radius_ECI();

  EXPECT_TRUE(abs(orbital_radius_perifocal - orbital_radius_ECI) < length_tolerance)
      << "Difference between evolved orbital radii calculated with "
      " perifocal and ECI coordinates: "
      << orbital_radius_perifocal - orbital_radius_ECI <<"\n";
}

TEST(EllipticalOrbitTests,OrbitalSpeedCalcs1) {
  Satellite test_satellite("../tests/elliptical_orbit_test_2.json");
  double orbital_speed_perifocal=test_satellite.get_speed();
  double orbital_speed_ECI=test_satellite.get_speed_ECI();
  EXPECT_TRUE(abs(orbital_speed_ECI - orbital_speed_perifocal) < tolerance)
      << "Difference between orbital speeds calculated with "
      " perifocal and ECI coordinates: "
      << orbital_speed_ECI - orbital_speed_perifocal <<"\n";
}

TEST(EllipticalOrbitTests,OrbitalSpeedCalcs2) {
  Satellite test_satellite("../tests/elliptical_orbit_test_2.json");
  double test_timestep = 0.1;  // s
  bool perturbation_bool = false;
  std::pair<double, int> new_timestep_and_error_code =
      test_satellite.evolve_RK45(epsilon, test_timestep, perturbation_bool);
  double next_timestep = new_timestep_and_error_code.first;
  int error_code = new_timestep_and_error_code.second;
  double orbital_speed_perifocal=test_satellite.get_speed();
  double orbital_speed_ECI=test_satellite.get_speed_ECI();

  EXPECT_TRUE(abs(orbital_speed_ECI - orbital_speed_perifocal) < tolerance)
      << "Difference between evolved orbital speeds calculated with "
      " perifocal and ECI coordinates: "
      << orbital_speed_ECI - orbital_speed_perifocal <<"\n";
}

TEST(EllipticalOrbitTests, TotalEnergyTimestep1) {
  Satellite test_satellite("../tests/elliptical_orbit_test_1.json");
  double initial_energy = test_satellite.get_total_energy();
  double test_timestep = 0.1;  // s
  bool perturbation_bool = true;
  std::pair<double, int> new_timestep_and_error_code =
      test_satellite.evolve_RK45(epsilon, test_timestep, perturbation_bool);
  double next_timestep = new_timestep_and_error_code.first;
  int error_code = new_timestep_and_error_code.second;
  double evolved_energy = test_satellite.get_total_energy();
  EXPECT_TRUE(abs(initial_energy - evolved_energy)/initial_energy < energy_cons_relative_tolerance)
      << "Total energy not preserved within relative tolerance. Relative difference: "
      << abs(initial_energy - evolved_energy)/initial_energy << "\n";
}

TEST(EllipticalOrbitTests, TotalEnergyTimestep2) {
  Satellite test_satellite("../tests/elliptical_orbit_test_2.json");
  double initial_energy = test_satellite.get_total_energy();
  double test_timestep = 0.1;  // s
  bool perturbation_bool = true;
  std::pair<double, int> new_timestep_and_error_code =
      test_satellite.evolve_RK45(epsilon, test_timestep, perturbation_bool);
  double next_timestep = new_timestep_and_error_code.first;
  int error_code = new_timestep_and_error_code.second;
  double evolved_energy = test_satellite.get_total_energy();
  EXPECT_TRUE(abs(initial_energy - evolved_energy)/initial_energy < energy_cons_relative_tolerance)
      << "Total energy not preserved within relative tolerance. Relative difference: "
      << abs(initial_energy - evolved_energy)/initial_energy << "\n";
}

TEST(EllipticalOrbitTests, TotalEnergyTimestep3) {
  Satellite test_satellite("../tests/elliptical_orbit_test_3.json");
  double initial_energy = test_satellite.get_total_energy();
  double test_timestep = 0.1;  // s
  bool perturbation_bool = true;
  std::pair<double, int> new_timestep_and_error_code =
      test_satellite.evolve_RK45(epsilon, test_timestep, perturbation_bool);
  double next_timestep = new_timestep_and_error_code.first;
  int error_code = new_timestep_and_error_code.second;
  double evolved_energy = test_satellite.get_total_energy();
  EXPECT_TRUE(abs(initial_energy - evolved_energy)/initial_energy < energy_cons_relative_tolerance)
      << "Total energy not preserved within relative tolerance. Relative difference: "
      << abs(initial_energy - evolved_energy)/initial_energy << "\n";
}

TEST(EllipticalOrbitTests, DragTest1) {
  Satellite test_satellite_withdrag("../tests/elliptical_orbit_test_4.json");
  Satellite test_satellite_nodrag("../tests/elliptical_orbit_test_4.json");
  // Drag parameters
  double F_10 = 100;  // Solar radio ten centimeter flux
  double A_p = 120;   // Geomagnetic A_p index
  double temp_epsilon = pow(10,-14);
  // Collect drag parameters into a pair with F_10 first and A_p second
  std::pair<double, double> drag_elements = {F_10, A_p};
  double test_timestep = 0.01;  // s
  bool perturbation_bool = true;
  double total_sim_time = 10; // s
  double current_time = test_satellite_nodrag.get_instantaneous_time();
  double orbital_radius = 0;
  while (current_time < total_sim_time) {
  std::pair<double, int> new_timestep_and_error_code =
      test_satellite_nodrag.evolve_RK45(temp_epsilon, test_timestep, perturbation_bool,
        false);
      orbital_radius = test_satellite_nodrag.get_radius();
      double altitude = (orbital_radius - radius_Earth)/1000.0;
      double next_timestep = new_timestep_and_error_code.first;
      test_timestep = next_timestep;
      int error_code = new_timestep_and_error_code.second;
      current_time = test_satellite_nodrag.get_instantaneous_time();
  }
  double no_drag_semimajor_axis = test_satellite_nodrag.get_orbital_element("Semimajor Axis");

  current_time = test_satellite_withdrag.get_instantaneous_time();

  while (current_time < total_sim_time) {
  std::pair<double, int> new_timestep_and_error_code =
  test_satellite_withdrag.evolve_RK45(temp_epsilon, test_timestep, perturbation_bool,
        true,drag_elements);
      double next_timestep = new_timestep_and_error_code.first;
      test_timestep = next_timestep;
      int error_code = new_timestep_and_error_code.second;
      current_time = test_satellite_withdrag.get_instantaneous_time();
  }
  double with_drag_semimajor_axis = test_satellite_withdrag.get_orbital_element("Semimajor Axis");

  EXPECT_TRUE(no_drag_semimajor_axis > with_drag_semimajor_axis)
      << "Semimajor axis after evolution wasn't lower when drag was introduced. "
      "This isn't expected behavior. Difference: " << no_drag_semimajor_axis - with_drag_semimajor_axis << "\n";
}

// Now trying to test the 140 < altitude < 180 km altitude section
TEST(EllipticalOrbitTests, DragTest2) {
  Satellite test_satellite_withdrag("../tests/circular_orbit_test_1_input.json");
  Satellite test_satellite_nodrag("../tests/circular_orbit_test_1_input.json");
  // Drag parameters
  double F_10 = 100;  // Solar radio ten centimeter flux
  double A_p = 120;   // Geomagnetic A_p index
  double temp_epsilon = pow(10,-14);
  // Collect drag parameters into a pair with F_10 first and A_p second
  std::pair<double, double> drag_elements = {F_10, A_p};
  double test_timestep = 0.01;  // s
  bool perturbation_bool = true;
  double total_sim_time = 10; // s
  double current_time = test_satellite_nodrag.get_instantaneous_time();
  double orbital_radius = 0;
  while (current_time < total_sim_time) {
  std::pair<double, int> new_timestep_and_error_code =
      test_satellite_nodrag.evolve_RK45(temp_epsilon, test_timestep, perturbation_bool,
        false);
      orbital_radius = test_satellite_nodrag.get_radius();
      double altitude = (orbital_radius - radius_Earth)/1000.0;
      double next_timestep = new_timestep_and_error_code.first;
      test_timestep = next_timestep;
      int error_code = new_timestep_and_error_code.second;
      current_time = test_satellite_nodrag.get_instantaneous_time();
  }
  double no_drag_semimajor_axis = test_satellite_nodrag.get_orbital_element("Semimajor Axis");

  current_time = test_satellite_withdrag.get_instantaneous_time();

  while (current_time < total_sim_time) {
  std::pair<double, int> new_timestep_and_error_code =
  test_satellite_withdrag.evolve_RK45(temp_epsilon, test_timestep, perturbation_bool,
        true,drag_elements);
      double next_timestep = new_timestep_and_error_code.first;
      test_timestep = next_timestep;
      int error_code = new_timestep_and_error_code.second;
      current_time = test_satellite_withdrag.get_instantaneous_time();
  }
  double with_drag_semimajor_axis = test_satellite_withdrag.get_orbital_element("Semimajor Axis");

  EXPECT_TRUE(no_drag_semimajor_axis > with_drag_semimajor_axis)
      << "Semimajor axis after evolution wasn't lower when drag was introduced. "
      "This isn't expected behavior. Difference: " << no_drag_semimajor_axis - with_drag_semimajor_axis << "\n";
}


// Circular orbit tests

// Testing calculated orbital speed based on input orbital parameters at e=0
// against known formula for circular orbital speed
// https://en.wikipedia.org/wiki/Circular_orbit#Velocity

TEST(CircularOrbitTests, OrbitalSpeed1) {
  Satellite test_satellite("../tests/circular_orbit_test_1_input.json");
  double calculated_radius = test_satellite.get_radius();
  double calculated_speed = test_satellite.get_speed();
  double expected_circular_orbital_speed =
      sqrt(G * mass_Earth / calculated_radius);
  EXPECT_TRUE(abs(calculated_speed - expected_circular_orbital_speed) <
              tolerance)
      << "Calculated orbital speed did not match expected value within "
         "tolerance. Difference: "
      << calculated_speed - expected_circular_orbital_speed << "\n";
}

TEST(CircularOrbitTests, OrbitalSpeed2) {
  Satellite test_satellite("../tests/circular_orbit_test_2_input.json");
  double calculated_radius = test_satellite.get_radius();
  double calculated_speed = test_satellite.get_speed();
  double expected_circular_orbital_speed =
      sqrt(G * mass_Earth / calculated_radius);
  EXPECT_TRUE(abs(calculated_speed - expected_circular_orbital_speed) <
              tolerance)
      << "Calculated orbital speed did not match expected value within "
         "tolerance. Difference: "
      << calculated_speed - expected_circular_orbital_speed << "\n";
}

TEST(CircularOrbitTests, TotalEnergyTimestep1) {
  Satellite test_satellite("../tests/circular_orbit_test_2_input.json");
  double initial_energy = test_satellite.get_total_energy();
  double test_timestep = 1;  // s
  bool perturbation_bool = true;
  std::pair<double, int> new_timestep_and_error_code =
      test_satellite.evolve_RK45(epsilon, test_timestep, perturbation_bool);
  double next_timestep = new_timestep_and_error_code.first;
  int error_code = new_timestep_and_error_code.second;
  double evolved_energy = test_satellite.get_total_energy();

  EXPECT_TRUE(abs(initial_energy - evolved_energy)/initial_energy < energy_cons_relative_tolerance)
      << "Total energy not preserved within relative tolerance. Relative difference: "
      << abs(initial_energy - evolved_energy)/initial_energy << "\n";
}

TEST(CircularOrbitTests, EvolvedOrbitalRadius1) {
  Satellite test_satellite("../tests/circular_orbit_test_2_input.json");
  double calculated_initial_radius = test_satellite.get_radius();
  double test_timestep = 1;
  bool perturbation_bool =
      false;  // While from what I can tell (see, e.g.,
              // https://ocw.tudelft.nl/wp-content/uploads/AE2104-Orbital-Mechanics-Slides_8.pdf)
              // there's no major effects on semimajor axis from J2
              // perturbation, not clear to me that this should be exactly
              // constant with J2 perturbation enabled

  std::pair<double, int> new_timestep_and_error_code =
      test_satellite.evolve_RK45(epsilon, test_timestep, perturbation_bool);
  double next_timestep = new_timestep_and_error_code.first;
  int error_code = new_timestep_and_error_code.second;
  double calculated_evolved_radius = test_satellite.get_radius();

  EXPECT_TRUE(abs(calculated_initial_radius - calculated_evolved_radius) <
              length_tolerance)
      << "Orbital radius not constant within tolerance. Difference: "
      << calculated_initial_radius - calculated_evolved_radius << "\n";
}

TEST(CircularOrbitTests, CircEvolvedOrbitalSpeed1) {
  Satellite test_satellite("../tests/circular_orbit_test_2_input.json");
  double calculated_initial_speed = test_satellite.get_speed();
  double test_timestep = 1;
  bool perturbation_bool =
      false;  // While from what I can tell (see, e.g.,
              // https://ocw.tudelft.nl/wp-content/uploads/AE2104-Orbital-Mechanics-Slides_8.pdf)
              // there's no major effects on semimajor axis from J2
              // perturbation, not clear to me that this should be exactly
              // constant with J2 perturbation enabled

  std::pair<double, int> new_timestep_and_error_code =
      test_satellite.evolve_RK45(epsilon, test_timestep, perturbation_bool);
  double next_timestep = new_timestep_and_error_code.first;
  int error_code = new_timestep_and_error_code.second;
  double calculated_evolved_speed = test_satellite.get_speed();

  EXPECT_TRUE(abs(calculated_initial_speed - calculated_evolved_speed) <
              tolerance)
      << "Orbital speed not constant within tolerance. Difference: "
      << calculated_initial_speed - calculated_evolved_speed << "\n";
}

TEST(CircularOrbitTests, BasicOrbitalElementsTest) {
  Satellite test_satellite("../tests/circular_orbit_test_2_input.json");
  std::array<double, 6> initial_orbit_elements =
      test_satellite.get_orbital_elements();

  test_satellite.update_orbital_elements_from_position_and_velocity();
  std::array<double, 6> recalculated_orbit_elements =
      test_satellite.get_orbital_elements();

  std::array<std::string, 6> orbital_element_name_array;
  orbital_element_name_array.at(0) = "Semimajor Axis";
  orbital_element_name_array.at(1) = "Eccentricity";
  orbital_element_name_array.at(2) = "Inclination";
  orbital_element_name_array.at(3) = "RAAN";
  orbital_element_name_array.at(4) = "Argument of Periapsis";
  orbital_element_name_array.at(5) = "True Anomaly";

  for (size_t orbital_elem_index = 0; orbital_elem_index < 6;
       orbital_elem_index++) {
    if (orbital_elem_index == 0) {
      EXPECT_TRUE(abs(initial_orbit_elements.at(orbital_elem_index) -
                      recalculated_orbit_elements.at(orbital_elem_index)) <
                  length_tolerance)
          << orbital_element_name_array.at(orbital_elem_index)
          << " was not constant within tolerance. Diff:"
          << initial_orbit_elements.at(orbital_elem_index) -
                 recalculated_orbit_elements.at(orbital_elem_index)
          << "\n";
    } else {
      EXPECT_TRUE(abs(initial_orbit_elements.at(orbital_elem_index) -
                      recalculated_orbit_elements.at(orbital_elem_index)) <
                  tolerance)
          << orbital_element_name_array.at(orbital_elem_index)
          << " was not constant within tolerance. Diff:"
          << initial_orbit_elements.at(orbital_elem_index) -
                 recalculated_orbit_elements.at(orbital_elem_index)
          << "\n";
    }
  }
}

TEST(CircularOrbitTests, ConstantEvolvedOrbitalElementsTest) {
  // The idea behind this test is that after evolving a timestep, orbital
  // elements besides true anomaly should be constant (when J2 perturbation is
  // not taken into account)
  Satellite test_satellite("../tests/circular_orbit_test_2_input.json");
  std::array<double, 6> initial_orbit_elements =
      test_satellite.get_orbital_elements();

  double test_timestep = 1;
  bool perturbation_bool = false;

  std::pair<double, int> new_timestep_and_error_code =
      test_satellite.evolve_RK45(epsilon, test_timestep, perturbation_bool);
  double next_timestep = new_timestep_and_error_code.first;
  int error_code = new_timestep_and_error_code.second;
  std::array<double, 6> evolved_orbit_elements =
      test_satellite.get_orbital_elements();

  std::array<std::string, 6> orbital_element_name_array;
  orbital_element_name_array.at(0) = "Semimajor Axis";
  orbital_element_name_array.at(1) = "Eccentricity";
  orbital_element_name_array.at(2) = "Inclination";
  orbital_element_name_array.at(3) = "RAAN";
  orbital_element_name_array.at(4) = "Argument of Periapsis";
  orbital_element_name_array.at(5) = "True Anomaly";

  for (size_t orbital_elem_index = 0; orbital_elem_index < 5;
       orbital_elem_index++) {
    // True anomaly shouldn't be constant over evolution
    if (orbital_elem_index == 0) {
      EXPECT_TRUE(abs(initial_orbit_elements.at(orbital_elem_index) -
                      evolved_orbit_elements.at(orbital_elem_index)) <
                  length_tolerance)
          << orbital_element_name_array.at(orbital_elem_index)
          << " was not constant within tolerance. Diff:"
          << initial_orbit_elements.at(orbital_elem_index) -
                 evolved_orbit_elements.at(orbital_elem_index)
          << "\n";
    } else {
      EXPECT_TRUE(abs(initial_orbit_elements.at(orbital_elem_index) -
                      evolved_orbit_elements.at(orbital_elem_index)) <
                  tolerance)
          << orbital_element_name_array.at(orbital_elem_index)
          << " was not constant within tolerance. Diff:"
          << initial_orbit_elements.at(orbital_elem_index) -
                 evolved_orbit_elements.at(orbital_elem_index)
          << "\n";
    }
  }
}

TEST(CircularOrbitTests, Thruster_Eccentricity_Change) {
  Satellite test_satellite("../tests/circular_orbit_test_2_input.json");
  std::array<double, 3> LVLH_thrust_direction = {1, 0, 0};
  double thrust_magnitude = 100;  // N
  double t_thrust_start = 1;
  double t_thrust_end = 100;

  test_satellite.add_LVLH_thrust_profile(
      LVLH_thrust_direction, thrust_magnitude, t_thrust_start, t_thrust_end);
  double test_timestep = 1;  // s
  double current_satellite_time = test_satellite.get_instantaneous_time();
  double sim_end_time = 110;
  while (current_satellite_time < sim_end_time) {
    std::pair<double, int> new_timestep_and_error_code =
        test_satellite.evolve_RK45(epsilon, test_timestep);
    double next_timestep = new_timestep_and_error_code.first;
    int error_code = new_timestep_and_error_code.second;
    test_timestep = next_timestep;
    current_satellite_time = test_satellite.get_instantaneous_time();
  }
  std::array<double, 6> evolved_orbit_elements =
      test_satellite.get_orbital_elements();
  double resulting_eccentricity = evolved_orbit_elements.at(1);

  EXPECT_TRUE(resulting_eccentricity > 0)
      << "Resulting eccentricity was not greater than 0. Calculated value: "
      << resulting_eccentricity << "\n";
}


// Attitude-related tests
const double pitch_tolerance = 5*pow(10.0, -3);

TEST(AttitudeTests, PassivePitchTest1) {
  // Without any external or initial torques, the satellite's pitch angle
  // w.r.t. the LVLH frame should progress 2*pi radians over one full orbit
  Satellite test_satellite("../tests/attitude_test_input_1.json");
  const double initial_true_anomaly =
      test_satellite.get_orbital_element("True Anomaly");
  const double initial_pitch = test_satellite.get_attitude_val("Pitch");
  const double orbital_period = test_satellite.calculate_orbital_period();
  double test_timestep = 0.1;  // s
  double current_true_anomaly = initial_true_anomaly;
  double next_timestep = 0;
  bool wrapped_around = false;
  while ((current_true_anomaly < initial_true_anomaly) || (wrapped_around == false)) {
    std::pair<double, int> new_timestep_and_error_code =
        test_satellite.evolve_RK45(epsilon, test_timestep);
    current_true_anomaly = test_satellite.get_orbital_element("True Anomaly");
    if ((!wrapped_around) && (0 <= current_true_anomaly) && (current_true_anomaly <= initial_true_anomaly)){
      wrapped_around = true;
    }
    next_timestep = new_timestep_and_error_code.first;
    int error_code = new_timestep_and_error_code.second;
    test_timestep = next_timestep;
  }
  double evolved_pitch = test_satellite.get_attitude_val("Pitch");
  EXPECT_TRUE(abs(evolved_pitch - initial_pitch) <
  pitch_tolerance)
      << "Pitch didn't progress 2*pi radians (within tolerance) "
      "over one orbit as expected. Difference: "
      << evolved_pitch - initial_pitch << "\n";
}

TEST(AttitudeTests, PassivePitchTest2) {
  // Let's make sure the pitch is progressing in the expected direction
  Satellite test_satellite("../tests/attitude_test_input_1.json");
  const double initial_true_anomaly =
      test_satellite.get_orbital_element("True Anomaly");
  const double initial_pitch = test_satellite.get_attitude_val("Pitch");
  const double time_to_sim = 1;  // s
  double test_timestep = 1;      // s
  double current_time = test_satellite.get_instantaneous_time();
  double next_timestep = 0;
  while (current_time < time_to_sim) {
    std::pair<double, int> new_timestep_and_error_code =
        test_satellite.evolve_RK45(epsilon, test_timestep);
    current_time = test_satellite.get_instantaneous_time();
    next_timestep = new_timestep_and_error_code.first;
    int error_code = new_timestep_and_error_code.second;
    test_timestep = next_timestep;
  }
  double evolved_pitch = test_satellite.get_attitude_val("Pitch");
  EXPECT_TRUE(evolved_pitch > initial_pitch)
      << "Pitch didn't increase as expected\n";
}


// Misc tests
TEST(MiscTests, ThrustProfileInitializationTest1) {
  double t_start = 1.0;
  double t_end = 9.0;
  std::array<double,3> thrust_vec = {1.0,101.2,-0.4};
  double magnitude = sqrt(pow(thrust_vec.at(0),2) + pow(thrust_vec.at(1),2) + pow(thrust_vec.at(2),2));
  std::array<double,3> thrust_vec_direction = {0.0,0.0,0.0};
  for (size_t ind=0;ind<thrust_vec.size();ind++) {
    thrust_vec_direction.at(ind) = thrust_vec.at(ind)/magnitude;
  }

  ThrustProfileLVLH thrust_profile_1(t_start,t_end,thrust_vec);
  ThrustProfileLVLH thrust_profile_2(t_start,t_end,thrust_vec_direction,magnitude);

  EXPECT_TRUE(thrust_profile_1 == thrust_profile_2)
      << "Thrust profiles initialized differently didn't agree.\n";
}

TEST(MiscTests, TorqueProfileInitializationTest1) {
  double t_start = 1.0;
  double t_end = 9.0;
  std::array<double,3> torque_vec = {1.0,101.2,-0.4};
  double magnitude = sqrt(pow(torque_vec.at(0),2) + pow(torque_vec.at(1),2) + pow(torque_vec.at(2),2));
  std::array<double,3> torque_vec_direction = {0.0,0.0,0.0};
  for (size_t ind=0;ind<torque_vec.size();ind++) {
    torque_vec_direction.at(ind) = torque_vec.at(ind)/magnitude;
  }

  BodyframeTorqueProfile torque_profile_1(t_start,t_end,torque_vec);
  BodyframeTorqueProfile torque_profile_2(t_start,t_end,torque_vec_direction,magnitude);

  EXPECT_TRUE(torque_profile_1 == torque_profile_2)
      << "Torque profiles initialized differently didn't agree.\n";
}

// Make sure zero-inclination orbit error is thrown correctly
TEST(MiscTests, ZeroInclinationTest) {
  bool caught = false;
  try {
    Satellite test_satellite("../tests/zero_inclination_test.json");
  }
  catch( const std::invalid_argument& invalid_arg_error) {
    caught = true;
  }

  EXPECT_TRUE(caught)
      << "Didn't catch the invalid argument error expected for zero inclination\n";
}

// Make sure function to get name of Satellite is working as expected
TEST(MiscTests, SatelliteNameTest) {
  Satellite test_satellite("../tests/circular_orbit_test_1_input.json");
  std::string expected_name = "Circ_Test_1";
  std::string recovered_name = test_satellite.get_name();


  EXPECT_TRUE(expected_name == recovered_name)
      << "Satellite name fetching didn't work as expected\n";
}
