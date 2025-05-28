#include <iostream>

#include "PhasedArrayGroundStation.h"
#include "Satellite.h"
#include "utils.h"

int main() {
  // This file demonstrates a few different ways you can run and
  // visualize data from satellite simulations.
  // First let's initialize the struct containing parameters for simulation and
  // plotting

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
  double test_timestep = 1.0;
  test_sat_1.add_LVLH_thrust_profile(LVLH_thrust_vec_2, t_thrust_2_start,
                                     t_thrust_2_end);
    size_t num_sim_steps = 100000000;
    const double temp_epsilon = sim_parameters.epsilon;
    for (size_t sim_step = 0; sim_step < num_sim_steps; sim_step++){
        std::pair<double, int> new_timestep_and_error_code =test_sat_1.evolve_RK45(temp_epsilon, test_timestep, false);
        double next_timestep = new_timestep_and_error_code.first;
        test_timestep = next_timestep;
        int error_code = new_timestep_and_error_code.second;
    }
  return 0;
}