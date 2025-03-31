#include <gtest/gtest.h>

#include <iostream>

#include "Satellite.h"
#include "utils.h"

const double tolerance = pow(10.0, -12);
const double epsilon = pow(10.0, -11);

TEST(AttitudeTests, PassivePitchTest1) {
  // Without any external or initial torques, the satellite's pitch angle
  // w.r.t. the LVLH frame should progress 2*pi radians over one full orbit
  Satellite test_satellite("../tests/attitude_test_input_1.json");
  const double initial_true_anomaly =
      test_satellite.get_orbital_element("True Anomaly");
  const double initial_pitch = test_satellite.get_attitude_val("Pitch");
  const double orbital_period = test_satellite.calculate_orbital_period();
  double test_timestep = 0.1;  // s
  double current_time = test_satellite.get_instantaneous_time();
  double next_timestep = 0;
  while (current_time < orbital_period) {
    std::pair<double, int> new_timestep_and_error_code =
        test_satellite.evolve_RK45(epsilon, test_timestep);
    current_time = test_satellite.get_instantaneous_time();
    next_timestep = new_timestep_and_error_code.first;
    int error_code = new_timestep_and_error_code.second;
    test_timestep = next_timestep;
  }
  double evolved_pitch = test_satellite.get_attitude_val("Pitch");
  double after_loop_time = test_satellite.get_instantaneous_time();
  EXPECT_TRUE(abs(evolved_pitch - initial_pitch) <=
              abs(after_loop_time - orbital_period))
      << "Pitch didn't progress 2*pi radians over one orbit as expected. "
         "Difference: "
      << evolved_pitch - initial_pitch
      << ", whereas the difference between the after-loop "
         "time and orbital period was "
      << after_loop_time - orbital_period << "\n";
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