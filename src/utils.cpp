#include "utils.h"

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>

#include "PhasedArrayGroundStation.h"
#include "Satellite.h"

using Eigen::Matrix3d;
using Eigen::Matrix4d;
using Eigen::Quaterniond;
using Eigen::Vector3d;
using Eigen::Vector4d;

Vector4d quaternion_multiplication(const Vector4d quaternion_1,
                                   const Vector4d quaternion_2) {
  // Ref: https://en.wikipedia.org/wiki/Quaternion#Scalar_and_vector_parts
  Vector3d quat_1_vec_component = {
      quaternion_1(1), quaternion_1(2),
      quaternion_1(3)};  // Here adopting convention that scalar part is first
                         // element of quaternion
  Vector3d quat_2_vec_component = {quaternion_2(1), quaternion_2(2),
                                   quaternion_2(3)};
  Vector3d vec_part_of_result =
      quaternion_1(0) * quat_2_vec_component +
      quaternion_2(0) * quat_1_vec_component +
      quat_1_vec_component.cross(quat_2_vec_component);
  Vector4d result;
  result << quaternion_1(0) * quaternion_2(0) -
                quat_1_vec_component.dot(quat_2_vec_component),
      vec_part_of_result(0), vec_part_of_result(1), vec_part_of_result(2);
  return result;
}

//"manual" version, via dot products and cross products with position and
// velocity vectors
std::array<double, 3> convert_LVLH_to_ECI_manual(
    const std::array<double, 3> input_LVLH_vec,
    const std::array<double, 3> input_position_vec,
    const std::array<double, 3> input_velocity_vec) {
  // LVLH x-axis is defined as in the direction of motion
  // LVLH z-axis is defined as pointing back towards Earth, so along the
  // reversed direction of the position vector from the center of the Earth
  // y-axis determined from a cross product

  Vector3d ECI_unit_vec_x = {1.0, 0.0, 0.0};
  Vector3d ECI_unit_vec_y = {0.0, 1.0, 0.0};
  Vector3d ECI_unit_vec_z = {0.0, 0.0, 1.0};

  std::array<double, 3> current_ECI_position_array = input_position_vec;
  Vector3d current_ECI_position_unit_vec;
  current_ECI_position_unit_vec << current_ECI_position_array.at(0),
      current_ECI_position_array.at(1), current_ECI_position_array.at(2);
  current_ECI_position_unit_vec
      .normalize();  // to make it actually a unit vector

  std::array<double, 3> current_ECI_velocity_array = input_velocity_vec;
  Vector3d current_ECI_velocity_unit_vec;
  current_ECI_velocity_unit_vec << current_ECI_velocity_array.at(0),
      current_ECI_velocity_array.at(1), current_ECI_velocity_array.at(2);
  current_ECI_velocity_unit_vec.normalize();

  Vector3d LVLH_x_unit_vec = current_ECI_velocity_unit_vec;
  Vector3d LVLH_z_unit_vec = (-1) * current_ECI_position_unit_vec;

  Vector3d LVLH_y_unit_vec = LVLH_z_unit_vec.cross(LVLH_x_unit_vec);
  // Should already be normalized, just in case though
  LVLH_y_unit_vec.normalize();

  double v_x_ECI = input_LVLH_vec.at(0) * ECI_unit_vec_x.dot(LVLH_x_unit_vec) +
                   input_LVLH_vec.at(1) * ECI_unit_vec_x.dot(LVLH_y_unit_vec) +
                   input_LVLH_vec.at(2) * ECI_unit_vec_x.dot(LVLH_z_unit_vec);

  double v_y_ECI = input_LVLH_vec.at(0) * ECI_unit_vec_y.dot(LVLH_x_unit_vec) +
                   input_LVLH_vec.at(1) * ECI_unit_vec_y.dot(LVLH_y_unit_vec) +
                   input_LVLH_vec.at(2) * ECI_unit_vec_y.dot(LVLH_z_unit_vec);

  double v_z_ECI = input_LVLH_vec.at(0) * ECI_unit_vec_z.dot(LVLH_x_unit_vec) +
                   input_LVLH_vec.at(1) * ECI_unit_vec_z.dot(LVLH_y_unit_vec) +
                   input_LVLH_vec.at(2) * ECI_unit_vec_z.dot(LVLH_z_unit_vec);

  std::array<double, 3> output_ECI_arr = {v_x_ECI, v_y_ECI, v_z_ECI};
  return output_ECI_arr;
}

// baselining cartesian coordinates in ECI frame

// std::array<double, 3> calculate_orbital_acceleration(
//     const std::array<double, 3> input_r_vec, const double
//     input_spacecraft_mass, const std::vector<std::array<double, 3>>
//         input_vec_of_force_vectors_in_ECI) {
//   // Note: this is the version used in the RK4 solver
//   // orbital acceleration = -G m_Earth/distance^3 * r_vec (just based on
//   // rearranging F=ma with a the acceleration due to gravitational attraction
//   // between satellite and Earth
//   // https://en.wikipedia.org/wiki/Newton%27s_law_of_universal_gravitation)
//   // going to be assuming Earth's position doesn't change for now
//   // also assuming Earth is spherical, can loosen this assumption in the
//   future
//   // note: this is in ECI frame

//   std::array<double, 3> acceleration_vec_due_to_gravity =
//       input_r_vec;  // shouldn't need to explicitly call copy function
//       because
//                     // input_r_vec is passed by value not ref

//   // F=ma
//   // a=F/m = (F_grav + F_ext)/m = (F_grav/m) + (F_ext/m) =
//   -G*M_Earth/distance^3
//   // + ...

//   const double distance = sqrt(input_r_vec.at(0) * input_r_vec.at(0) +
//                                input_r_vec.at(1) * input_r_vec.at(1) +
//                                input_r_vec.at(2) * input_r_vec.at(2));
//   const double overall_factor = -G * mass_Earth / pow(distance, 3);

//   for (size_t ind = 0; ind < input_r_vec.size(); ind++) {
//     acceleration_vec_due_to_gravity.at(ind) *= overall_factor;
//   }

//   // now add effects from externally-applied forces, e.g., thrusters, if any
//   std::array<double, 3> acceleration_vec = acceleration_vec_due_to_gravity;

//   for (std::array<double, 3> external_force_vec_in_ECI :
//        input_vec_of_force_vectors_in_ECI) {
//     acceleration_vec.at(0) +=
//         (external_force_vec_in_ECI.at(0) / input_spacecraft_mass);
//     acceleration_vec.at(1) +=
//         (external_force_vec_in_ECI.at(1) / input_spacecraft_mass);
//     acceleration_vec.at(2) +=
//         (external_force_vec_in_ECI.at(2) / input_spacecraft_mass);
//   }
//   return acceleration_vec;
// }

std::array<double, 3> convert_cylindrical_to_cartesian(
    const double input_r_comp, const double input_theta_comp,
    const double input_z_comp, const double input_theta) {
  std::array<double, 3> output_cartesian_vec = {0, 0, 0};

  // Dot product method
  // See
  // https://en.wikipedia.org/wiki/Vector_fields_in_cylindrical_and_spherical_coordinates
  // for relation between unit vectors
  double x_comp =
      input_r_comp * cos(input_theta) + input_theta_comp * (-sin(input_theta));
  double y_comp =
      input_r_comp * sin(input_theta) + input_theta_comp * (cos(input_theta));
  double z_comp = input_z_comp;

  output_cartesian_vec.at(0) = x_comp;
  output_cartesian_vec.at(1) = y_comp;
  output_cartesian_vec.at(2) = z_comp;
  return output_cartesian_vec;
}

std::array<double, 3> calculate_orbital_acceleration(
    const std::array<double, 3> input_r_vec, const double input_spacecraft_mass,
    const std::vector<ThrustProfileLVLH> input_list_of_thrust_profiles_LVLH,
    const double input_evaluation_time,
    const std::array<double, 3> input_velocity_vec,
    const double input_inclination, const double input_arg_of_periapsis,
    const double input_true_anomaly, const double input_F_10,
    const double input_A_p, const double input_A_s,
    const double input_satellite_mass, const bool perturbation,
    const bool atmospheric_drag) {
  // Note: this is the version used in the RK45 solver (this has a more updated
  // workflow) orbital acceleration = -G m_Earth/distance^3 * r_vec (just based
  // on rearranging F=ma with a the acceleration due to gravitational attraction
  // between satellite and Earth
  // https://en.wikipedia.org/wiki/Newton%27s_law_of_universal_gravitation)
  // going to be assuming Earth's position doesn't change for now
  // also assuming Earth is spherical, can loosen this assumption in the future
  // note: this is in ECI frame (r_vec and velocity vec should also be in ECI
  // frame)

  std::array<double, 3> acceleration_vec_due_to_gravity =
      input_r_vec;  // shouldn't need to explicitly call copy function because
                    // input_r_vec is passed by value not ref

  // F=ma
  // a=F/m = (F_grav + F_ext)/m = (F_grav/m) + (F_ext/m) = -G*M_Earth/distance^3
  // + ...

  const double distance = sqrt(input_r_vec.at(0) * input_r_vec.at(0) +
                               input_r_vec.at(1) * input_r_vec.at(1) +
                               input_r_vec.at(2) * input_r_vec.at(2));
  const double overall_factor = -G * mass_Earth / pow(distance, 3);

  for (size_t ind = 0; ind < input_r_vec.size(); ind++) {
    acceleration_vec_due_to_gravity.at(ind) *= overall_factor;
  }

  std::array<double, 3> acceleration_vec = acceleration_vec_due_to_gravity;

  // now add effects from externally-applied forces, e.g., thrusters, if any

  std::vector<std::array<double, 3>> list_of_LVLH_forces_at_evaluation_time =
      {};
  std::vector<std::array<double, 3>> list_of_ECI_forces_at_evaluation_time = {};

  for (const ThrustProfileLVLH thrust_profile :
       input_list_of_thrust_profiles_LVLH) {
    if ((input_evaluation_time >= thrust_profile.t_start_) &&
        (input_evaluation_time <= thrust_profile.t_end_) &&
        (thrust_profile.arg_of_periapsis_change_thrust_profile == false)) {
      list_of_LVLH_forces_at_evaluation_time.push_back(
          thrust_profile.LVLH_force_vec_);
      std::array<double, 3> ECI_thrust_vector = convert_LVLH_to_ECI_manual(
          thrust_profile.LVLH_force_vec_, input_r_vec, input_velocity_vec);
      list_of_ECI_forces_at_evaluation_time.push_back(ECI_thrust_vector);
    }
  }

  for (std::array<double, 3> external_force_vec_in_ECI :
       list_of_ECI_forces_at_evaluation_time) {
    acceleration_vec.at(0) +=
        (external_force_vec_in_ECI.at(0) / input_spacecraft_mass);
    acceleration_vec.at(1) +=
        (external_force_vec_in_ECI.at(1) / input_spacecraft_mass);
    acceleration_vec.at(2) +=
        (external_force_vec_in_ECI.at(2) / input_spacecraft_mass);
  }

  if (perturbation) {
    // If accounting for J2 perturbation

    // Now let's add the additional acceleration components due to the J2
    // perturbation Ref:
    // https://vatankhahghadim.github.io/AER506/Notes/6%20-%20Orbital%20Perturbations.pdf
    double J2 = 1.083 * pow(10, -3);
    const double mu_Earth = G * mass_Earth;
    double C = 3 * mu_Earth * J2 * radius_Earth * radius_Earth /
               (2 * pow(distance, 4));
    double x = input_r_vec.at(0);
    double y = input_r_vec.at(1);
    double rho = sqrt(pow(x, 2) + pow(y, 2));
    double theta;
    // Ref:
    // https://en.wikipedia.org/wiki/Cylindrical_coordinate_system#Line_and_volume_elements
    if (x >= 0) {
      theta = asin(
          y / rho);  // Note: here setting theta=0 even if both x and y are
                     // zero, whereas it should technically be indeterminate,
                     // but I'm going to assume this edge condition won't be hit
                     // and I don't want undefined behavior
    } else {
      if (y >= 0) {
        theta = -asin(y / rho) + M_PI;
      } else {
        theta = -asin(y / rho) - M_PI;
      }
    }

    double a_r =
        C * (3 * pow(sin(input_inclination), 2) *
                 pow(sin(input_arg_of_periapsis + input_true_anomaly), 2) -
             1);
    double a_theta = -C * pow(sin(input_inclination), 2) *
                     sin(2 * (input_arg_of_periapsis + input_true_anomaly));
    double a_z = -C * sin(2 * input_inclination) *
                 sin(input_arg_of_periapsis + input_true_anomaly);
    std::array<double, 3> cartesian_acceleration_components =
        convert_cylindrical_to_cartesian(a_r, a_theta, a_z, theta);

    for (size_t ind = 0; ind < 3; ind++) {
      acceleration_vec.at(ind) += cartesian_acceleration_components.at(ind);
    }
  }
  double altitude = (distance - radius_Earth) / 1000.0;  // km
  if ((atmospheric_drag) && (altitude >= 140) && (altitude <= 400)) {
    // Refs: https://angeo.copernicus.org/articles/39/397/2021/
    // https://www.spaceacademy.net.au/watch/debris/atmosmod.htm
    double speed = sqrt(pow(input_velocity_vec.at(0), 2) +
                        pow(input_velocity_vec.at(1), 2) +
                        pow(input_velocity_vec.at(2), 2));
    // First, esimate atmospheric density
    double rho = {0};
    if (altitude < 180) {
      double a0 = 7.001985 * pow(10, -2);
      double a1 = -4.336216 * pow(10, -3);
      double a2 = -5.009831 * pow(10, -3);
      double a3 = 1.621827 * pow(10, -4);
      double a4 = -2.471283 * pow(10, -6);
      double a5 = 1.904383 * pow(10, -8);
      double a6 = -7.189421 * pow(10, -11);
      double a7 = 1.060067 * pow(10, -13);
      double fit_val = a0 + a1 * altitude + a2 * pow(altitude, 2) +
                       a3 * pow(altitude, 3) + a4 * pow(altitude, 4) +
                       a5 * pow(altitude, 5) + a6 * pow(altitude, 6) +
                       a7 * pow(altitude, 7);

      rho = pow(10, fit_val);
    } else {
      double T = 900 + 2.5 * (input_F_10 - 70) + 1.5 * input_A_p;
      double new_mu = 27 - 0.012 * (altitude - 200);
      double H = T / new_mu;
      rho = 6 * pow(10, -10) * exp(-(altitude - 175) / H);
    }

    // Now estimate the satellite's ballistic coefficient B
    double C_d = 2.2;
    double B = C_d * input_A_s / input_satellite_mass;
    double drag_deceleration = (1.0 / 2.0) * rho * B * pow(speed, 2);
    // Should act in direction directly opposite to velocity
    std::array<double, 3> velocity_unit_vec = {0.0, 0.0, 0.0};
    for (size_t ind = 0; ind < 3; ind++) {
      velocity_unit_vec.at(ind) = input_velocity_vec.at(ind) / speed;
    }
    std::array<double, 3> drag_acceleration_vec = {0.0, 0.0, 0.0};
    for (size_t ind = 0; ind < 3; ind++) {
      drag_acceleration_vec.at(ind) =
          drag_deceleration * (-1) * velocity_unit_vec.at(ind);
      // Factor of (-1) because this acceleration acts in direction opposite to
      // velocity
    }
    for (size_t ind = 0; ind < 3; ind++) {
      acceleration_vec.at(ind) += drag_acceleration_vec.at(ind);
    }
  }

  return acceleration_vec;
}

// std::array<double, 6> RK4_deriv_function_orbit_position_and_velocity(
//     const std::array<double, 6> input_position_and_velocity,
//     const double input_spacecraft_mass,
//     const std::vector<std::array<double, 3>>
//         input_vec_of_force_vectors_in_ECI) {
//   std::array<double, 6> derivative_of_input_y = {};
//   std::array<double, 3> position_array = {};

//   for (size_t ind = 0; ind < 3; ind++) {
//     derivative_of_input_y.at(ind) = input_position_and_velocity.at(ind + 3);
//     position_array.at(ind) = input_position_and_velocity.at(ind);
//   }

//   std::array<double, 3> calculated_orbital_acceleration =
//       calculate_orbital_acceleration(position_array, input_spacecraft_mass,
//                                      input_vec_of_force_vectors_in_ECI);

//   for (size_t ind = 3; ind < 6; ind++) {
//     derivative_of_input_y.at(ind) = calculated_orbital_acceleration.at(ind -
//     3);
//   }

//   return derivative_of_input_y;
// }

std::array<double, 6> RK45_deriv_function_orbit_position_and_velocity(
    const std::array<double, 6> input_position_and_velocity,
    const double input_spacecraft_mass,
    const std::vector<ThrustProfileLVLH> input_list_of_thrust_profiles_LVLH,
    const double input_evaluation_time, const double input_inclination,
    const double input_arg_of_periapsis, const double input_true_anomaly,
    const double input_F_10, const double input_A_p, const double input_A_s,
    const double input_satellite_mass, const bool perturbation,
    const bool atmospheric_drag) {
  // std::cout << "In RK45 deriv function orbit position and velocity, received
  // perturbation bool: " << perturbation << "\n";
  std::array<double, 6> derivative_of_input_y = {};
  std::array<double, 3> position_array = {};
  std::array<double, 3> velocity_array = {};

  for (size_t ind = 0; ind < 3; ind++) {
    derivative_of_input_y.at(ind) = input_position_and_velocity.at(ind + 3);
    velocity_array.at(ind) = input_position_and_velocity.at(ind + 3);
    position_array.at(ind) = input_position_and_velocity.at(ind);
  }

  std::array<double, 3> calculated_orbital_acceleration =
      calculate_orbital_acceleration(
          position_array, input_spacecraft_mass,
          input_list_of_thrust_profiles_LVLH, input_evaluation_time,
          velocity_array, input_inclination, input_arg_of_periapsis,
          input_true_anomaly, input_F_10, input_A_p, input_A_s,
          input_satellite_mass, perturbation, atmospheric_drag);

  for (size_t ind = 3; ind < 6; ind++) {
    derivative_of_input_y.at(ind) = calculated_orbital_acceleration.at(ind - 3);
  }

  return derivative_of_input_y;
}

/// @brief Runs the requested simulation and writes output data to a text file
/// @param input_satellite_vector: The vector of satellite objects to be simulated
/// @param input_sim_parameters: A SimParameters struct containing simulation parameters
/// @param output_file_name_prefix: The prefix of the text file containing simulation data. Individual satellite data files will start with this prefix, followed by the satellite name.
int sim_and_write_to_file(std::vector<Satellite> input_satellite_vector,
                                const SimParameters& input_sim_parameters,
                                const std::string output_file_name_prefix,
                                const double size_limit) {
  // Estimate of bytes of output file (NOTE: size on disk may be different): 372.67*n - 189.97 with n the number of timesteps at which there exists data
  // This would change if, e.g., the number of values stored at each timestep changed, or the precision at which these values are written are changed
  double size_estimate_slope = 372.67;
  double size_estimate_y_intercept = -189.97;
  double size_estimate = size_estimate_y_intercept; // bytes (this isn't physical at n = 0, but need to account for possibility that no data is written)
  // One satellite at a time
  for (Satellite current_satellite : input_satellite_vector) {
    std::string satellite_name = current_satellite.get_name();
    std::string file_name = "../data/" + output_file_name_prefix + satellite_name;
    std::ofstream output_file(file_name.c_str());
    // Get initial values of all stored parameters
    double current_satellite_time =
        current_satellite.get_instantaneous_time();
    // ECI Position
    std::array<double,3> current_satellite_position_ECI = current_satellite.get_ECI_position();
    double current_satellite_position_ECI_x = current_satellite_position_ECI.at(0);
    double current_satellite_position_ECI_y = current_satellite_position_ECI.at(1);
    double current_satellite_position_ECI_z = current_satellite_position_ECI.at(2);
    // ECI Velocity
    std::array<double,3> current_satellite_velocity_ECI = current_satellite.get_ECI_velocity();
    double current_satellite_velocity_ECI_x = current_satellite_velocity_ECI.at(0);
    double current_satellite_velocity_ECI_y = current_satellite_velocity_ECI.at(1);
    double current_satellite_velocity_ECI_z = current_satellite_velocity_ECI.at(2);
    // Attitude-related values
    double current_satellite_roll = current_satellite.get_attitude_val("Roll");
    double current_satellite_pitch = current_satellite.get_attitude_val("Pitch");
    double current_satellite_yaw = current_satellite.get_attitude_val("Yaw");
    double current_satellite_omega_x = current_satellite.get_attitude_val("omega_x");
    double current_satellite_omega_y = current_satellite.get_attitude_val("omega_y");
    double current_satellite_omega_z = current_satellite.get_attitude_val("omega_z");
    double current_satellite_q_0 = current_satellite.get_attitude_val("q_0");
    double current_satellite_q_1 = current_satellite.get_attitude_val("q_1");
    double current_satellite_q_2 = current_satellite.get_attitude_val("q_2");
    double current_satellite_q_3 = current_satellite.get_attitude_val("q_3");
    // Orbit-related parameters
    double current_satellite_a = current_satellite.get_orbital_parameter("Semimajor Axis");
    double current_satellite_eccentricity = current_satellite.get_orbital_parameter("Eccentricity");
    double current_satellite_inclination = current_satellite.get_orbital_parameter("Inclination");
    double current_satellite_raan = current_satellite.get_orbital_parameter("RAAN");
    double current_satellite_arg_of_periapsis = current_satellite.get_orbital_parameter("Argument of Periapsis");
    double current_satellite_true_anomaly = current_satellite.get_orbital_parameter("True Anomaly");
    double current_satellite_orbital_rate = current_satellite.get_orbital_parameter("Orbital Rate");
    double current_satellite_orbital_ang_acc = current_satellite.get_orbital_parameter("Orbital Angular Acceleration");
    double current_satellite_total_energy = current_satellite.get_orbital_parameter("Total Energy");

    int num_digits = 15;

    // Write initial datapoints
    output_file << std::setprecision(num_digits) 
    << current_satellite_time << " "
    << current_satellite_position_ECI_x << " "
    << current_satellite_position_ECI_y << " "
    << current_satellite_position_ECI_z << " "
    << current_satellite_velocity_ECI_x << " "
    << current_satellite_velocity_ECI_y << " "
    << current_satellite_velocity_ECI_z << " "
    << current_satellite_roll << " "
    << current_satellite_pitch << " "
    << current_satellite_yaw << " "
    << current_satellite_omega_x << " "
    << current_satellite_omega_y << " "
    << current_satellite_omega_z << " "
    << current_satellite_q_0 << " "
    << current_satellite_q_1 << " "
    << current_satellite_q_2 << " "
    << current_satellite_q_3 << " "
    << current_satellite_a << " "
    << current_satellite_eccentricity << " "
    << current_satellite_inclination << " "
    << current_satellite_raan << " "
    << current_satellite_arg_of_periapsis << " "
    << current_satellite_true_anomaly << " "
    << current_satellite_orbital_rate << " "
    << current_satellite_orbital_ang_acc << " "
    << current_satellite_total_energy << "\n";
    size_estimate += size_estimate_slope;
    if (size_estimate >= size_limit) {
      std::cout << "Hit estimated number of bytes above input limit, exiting early\n";
      output_file.close();
      return 2;
    }

    double timestep_to_use = input_sim_parameters.initial_timestep_guess;
    while (current_satellite_time < input_sim_parameters.total_sim_time) {

      std::pair<double, double> drag_elements = {input_sim_parameters.F_10,
                                                  input_sim_parameters.A_p};
      std::pair<double, int> new_timestep_and_error_code =
          current_satellite.evolve_RK45(
              input_sim_parameters.epsilon, timestep_to_use,
              input_sim_parameters.perturbation_bool,
              input_sim_parameters.drag_bool, drag_elements);
      double new_timestep = new_timestep_and_error_code.first;
      int error_code = new_timestep_and_error_code.second;

      if (error_code != 0) {
        std::cout << "Error code " << error_code
                  << " detected, halting simulation\n";
        output_file << "Error code " << error_code
                  << " detected, halting simulation\n";
        output_file.close();
        return 1;
      }
      timestep_to_use = new_timestep;
      current_satellite_time = current_satellite.get_instantaneous_time();
      // Get updated parameters
      // ECI Position
        current_satellite_position_ECI = current_satellite.get_ECI_position();
        current_satellite_position_ECI_x = current_satellite_position_ECI.at(0);
        current_satellite_position_ECI_y = current_satellite_position_ECI.at(1);
        current_satellite_position_ECI_z = current_satellite_position_ECI.at(2);
        // ECI Velocity
        current_satellite_velocity_ECI = current_satellite.get_ECI_velocity();
        current_satellite_velocity_ECI_x = current_satellite_velocity_ECI.at(0);
        current_satellite_velocity_ECI_y = current_satellite_velocity_ECI.at(1);
        current_satellite_velocity_ECI_z = current_satellite_velocity_ECI.at(2);
        // Attitude-related values
        current_satellite_roll = current_satellite.get_attitude_val("Roll");
        current_satellite_pitch = current_satellite.get_attitude_val("Pitch");
        current_satellite_yaw = current_satellite.get_attitude_val("Yaw");
        current_satellite_omega_x = current_satellite.get_attitude_val("omega_x");
        current_satellite_omega_y = current_satellite.get_attitude_val("omega_y");
        current_satellite_omega_z = current_satellite.get_attitude_val("omega_z");
        current_satellite_q_0 = current_satellite.get_attitude_val("q_0");
        current_satellite_q_1 = current_satellite.get_attitude_val("q_1");
        current_satellite_q_2 = current_satellite.get_attitude_val("q_2");
        current_satellite_q_3 = current_satellite.get_attitude_val("q_3");
        // Orbit-related parameters
        current_satellite_a = current_satellite.get_orbital_parameter("Semimajor Axis");
        current_satellite_eccentricity = current_satellite.get_orbital_parameter("Eccentricity");
        current_satellite_inclination = current_satellite.get_orbital_parameter("Inclination");
        current_satellite_raan = current_satellite.get_orbital_parameter("RAAN");
        current_satellite_arg_of_periapsis = current_satellite.get_orbital_parameter("Argument of Periapsis");
        current_satellite_true_anomaly = current_satellite.get_orbital_parameter("True Anomaly");
        current_satellite_orbital_rate = current_satellite.get_orbital_parameter("Orbital Rate");
        current_satellite_orbital_ang_acc = current_satellite.get_orbital_parameter("Orbital Angular Acceleration");
        current_satellite_total_energy = current_satellite.get_orbital_parameter("Total Energy");

        // Write updated parameters to file
        output_file << std::setprecision(num_digits) 
        << current_satellite_time << " "
        << current_satellite_position_ECI_x << " "
        << current_satellite_position_ECI_y << " "
        << current_satellite_position_ECI_z << " "
        << current_satellite_velocity_ECI_x << " "
        << current_satellite_velocity_ECI_y << " "
        << current_satellite_velocity_ECI_z << " "
        << current_satellite_roll << " "
        << current_satellite_pitch << " "
        << current_satellite_yaw << " "
        << current_satellite_omega_x << " "
        << current_satellite_omega_y << " "
        << current_satellite_omega_z << " "
        << current_satellite_q_0 << " "
        << current_satellite_q_1 << " "
        << current_satellite_q_2 << " "
        << current_satellite_q_3 << " "
        << current_satellite_a << " "
        << current_satellite_eccentricity << " "
        << current_satellite_inclination << " "
        << current_satellite_raan << " "
        << current_satellite_arg_of_periapsis << " "
        << current_satellite_true_anomaly << " "
        << current_satellite_orbital_rate << " "
        << current_satellite_orbital_ang_acc << " "
        << current_satellite_total_energy << "\n";
        size_estimate += size_estimate_slope;
        if (size_estimate >= size_limit) {
          std::cout << "Hit estimated number of bytes above input limit, exiting early\n";
          output_file.close();
          return 2;
        }
    }
    output_file.close();
  }
  return 0;
}


/// @brief Uses Gnuplot to plot a specified parameter from a specified simulation data file.
/// @param data_file_name: Name of the data file to use
/// @param parameter_to_plot: Parameter you would like to plot. See comments above implementation for details.

/*
Possible values of "parameter_to_plot":
"3D Position" - makes an interactive 3D plot of the satellite's ECI position [m]. This is the only 3D plot currently supported.alignas
The following are available as 2D plots (the specified parameter will be plotted with respect to simulation time)
"X Velocity" - X component of ECI velocity [m/s]
"Y Velocity" - Y component of ECI velocity [m/s]
"Z Velocity" - Z component of ECI velocity [m/s]
"Roll" - Satellite's roll angle with respect to LVLH frame, represented in the satellite body frame [deg]
"Pitch" - Satellite's pitch angle with respect to LVLH frame, represented in the satellite body frame [deg]
"Yaw" - Satellite's yaw angle with respect to LVLH frame, represented in the satellite body frame [deg]
"omega_x" - x component of satellite's body-frame angular velocity with respect to the LVLH frame [deg/s]
"omega_y" - y component of satellite's body-frame angular velocity with respect to the LVLH frame [deg/s]
"omega_z" - z component of satellite's body-frame angular velocity with respect to the LVLH frame [deg/s]
"q_0" - First (scalar) component of quaternion representing satellite's attitude with respect to the LVLH frame
"q_1" - Second component of quaternion representing satellite's attitude with respect to the LVLH frame
"q_2" - Third component of quaternion representing satellite's attitude with respect to the LVLH frame
"q_3" - Fourth component of quaternion representing satellite's attitude with respect to the LVLH frame
"Semimajor Axis" - Semimajor axis of satellite's orbit [m]
"Eccentricity" - Eccentricity of satellite's orbit
"Inclination" - Inclination of satellite's orbit [deg]
"RAAN" - Right ascension of the ascending node of satellite's orbit [deg]
        (assumed to be equivalent to the longitude of the ascending node for Earth-orbiting satellites simulated here)
"Argument of Periapsis" - Argument of periapsis of the satellite's orbit [deg]
"True Anomaly" - True anomaly of the satellite's orbit [deg]
"Orbital Rate" - Rate of change of the true anomaly
"Orbital Angular Acceleration" - Second time derivative of the true anomaly
"Total Energy" - Satellite's total energy 
                (here taken to be gravitational potential energy + translational kinetic energy + rotational kinetic energy)
*/
void plot_3D_from_datafile(std::vector<std::string> input_datafile_name_vector,
                                const std::string output_file_name,
                                const std::string terminal_name,
                                const double x_increment,
                                const double y_increment,
                                const double z_increment) {
  if (input_datafile_name_vector.size() < 1) {
    throw std::invalid_argument("Empty datafile name vector");
  }
  // first, open "pipe" to gnuplot
  std::string gnuplot_arg_string = "gnuplot";
  if (terminal_name == "qt") {
    gnuplot_arg_string += " -persist";
  }
  FILE* gnuplot_pipe = popen(gnuplot_arg_string.c_str(), "w");

  // if it exists
  if (gnuplot_pipe) {
    fprintf(gnuplot_pipe, "set terminal '%s' size 900,700 font ',14'\n",
            terminal_name.c_str());
    if (terminal_name == "png") {
      fprintf(gnuplot_pipe, "set output '../%s.png'\n",
              output_file_name.c_str());
    }
    // formatting
    fprintf(gnuplot_pipe, "set xlabel 'x [m]' offset 0,-2\n");
    fprintf(gnuplot_pipe, "set ylabel 'y [m]' offset -2,0\n");
    fprintf(gnuplot_pipe, "set zlabel 'z [m]'\n");
    // fprintf(gnuplot_pipe,"set view 70,1,1,1\n");
    fprintf(gnuplot_pipe, "set view equal xyz\n");
    if (x_increment != 0) {
      fprintf(gnuplot_pipe, "set xtics %e offset 0,-1.5\n",
              x_increment);
    } else {
      fprintf(gnuplot_pipe, "set xtics offset 0,-1.5\n");
    }
    if (y_increment != 0) {
      fprintf(gnuplot_pipe, "set ytics %e offset -1.5,0\n",
              y_increment);
    } else {
      fprintf(gnuplot_pipe, "set ytics offset -1.5,0\n");
    }
    if (z_increment != 0) {
      fprintf(gnuplot_pipe, "set ztics %e\n", z_increment);
    }
    fprintf(gnuplot_pipe, "unset colorbox\n");
    fprintf(gnuplot_pipe, "set style fill transparent solid 1.0\n");

    fprintf(gnuplot_pipe, "set key offset 0,-10\n");
    fprintf(gnuplot_pipe, "set hidden3d front\n");

    // plotting
    // first let's set the stage for plotting the Earth
    fprintf(gnuplot_pipe, "R_Earth=%.17g\n", radius_Earth);
    fprintf(gnuplot_pipe, "set isosamples 50,50\n");
    fprintf(gnuplot_pipe, "set parametric\n");
    fprintf(gnuplot_pipe, "set urange [-pi/2:pi/2]\n");
    fprintf(gnuplot_pipe, "set vrange [0:2*pi]\n");

    // first satellite
    std::string current_file_name = input_datafile_name_vector.at(0);
    if (input_datafile_name_vector.size() == 1) {
      fprintf(gnuplot_pipe, "splot '../data/%s' using 2:3:4 with lines lw 1 title '%s' \\\n",
              current_file_name.c_str(),
              current_file_name.c_str());
    }

    else {
      fprintf(gnuplot_pipe, "splot '../data/%s' using 2:3:4 with lines lw 1 title '%s'\\\n",
              current_file_name.c_str(),
              current_file_name.c_str());
    }

    for (size_t datafile_index = 1;
         datafile_index < input_datafile_name_vector.size(); datafile_index++) {
      current_file_name = input_datafile_name_vector.at(datafile_index);
      if (datafile_index < input_datafile_name_vector.size() - 1) {
        fprintf(gnuplot_pipe, ",'../data/%s' using 2:3:4 with lines lw 1 title '%s' \\\n",
                current_file_name.c_str(),
                current_file_name.c_str());
      }
      else {
        fprintf(gnuplot_pipe, ",'../data/%s' using 2:3:4 with lines lw 1 title '%s' \\\n",
                current_file_name.c_str(),
                current_file_name.c_str());
        }
    }
    fprintf(gnuplot_pipe,
            ",R_Earth*cos(u)*cos(v),R_Earth*cos(u)*sin(v),R_Earth*sin(u) "
            "notitle with pm3d fillcolor rgbcolor 'navy'\n");

    if (terminal_name == "qt") {
      fprintf(gnuplot_pipe, "pause mouse keypress\n");
    }
    fprintf(gnuplot_pipe, "exit \n");
    pclose(gnuplot_pipe);
    std::cout << "Done\n";

  } else {
    std::cout << "gnuplot not found";
  }

  return;
}



// Now the 2D version for plotting parameters over time
void plot_2D_from_datafile(std::vector<std::string> input_datafile_name_vector,
                            const std::string plotted_parameter,
                            const std::string output_file_name) {
  if (input_datafile_name_vector.size() < 1) {
    throw std::invalid_argument("Empty datafile name vector");
  }

  std::unordered_map<std::string,std::string> name_to_units_map;
  name_to_units_map["X Velocity"] = "[m/s]";
  name_to_units_map["Y Velocity"] = "[m/s]";
  name_to_units_map["Z Velocity"] = "[m/s]";
  name_to_units_map["Semimajor Axis"] = "[m]";
  name_to_units_map["Eccentricity"] = "";
  name_to_units_map["Inclination"] = "[deg]";
  name_to_units_map["Argument of Periapsis"] = "[deg]";
  name_to_units_map["RAAN"] = "[deg]";
  name_to_units_map["True Anomaly"] = "[deg]";
  name_to_units_map["Roll"] = "[deg]";
  name_to_units_map["Pitch"] = "[deg]";
  name_to_units_map["Yaw"] = "[deg]";
  name_to_units_map["omega_x"] = "[deg/s]";
  name_to_units_map["omega_y"] = "[deg/s]";
  name_to_units_map["omega_z"] = "[deg/s]";
  name_to_units_map["q_0"] = "";
  name_to_units_map["q_1"] = "";
  name_to_units_map["q_2"] = "";
  name_to_units_map["q_3"] = "";
  name_to_units_map["Total Energy"] = "[J]";
  name_to_units_map["Orbital Rate"] = "[rad/s]";
  name_to_units_map["Orbital Angular Acceleration"] = "[rad/s^2]";

  std::unordered_map<std::string,int> name_to_col_map;
  name_to_col_map["X Velocity"] = 5;
  name_to_col_map["Y Velocity"] = 6;
  name_to_col_map["Z Velocity"] = 7;
  name_to_col_map["Roll"] = 8;
  name_to_col_map["Pitch"] = 9;
  name_to_col_map["Yaw"] = 10;
  name_to_col_map["omega_x"] = 11;
  name_to_col_map["omega_y"] = 12;
  name_to_col_map["omega_z"] = 13;
  name_to_col_map["q_0"] = 14;
  name_to_col_map["q_1"] = 15;
  name_to_col_map["q_2"] = 16;
  name_to_col_map["q_3"] = 17;
  name_to_col_map["Semimajor Axis"] = 18;
  name_to_col_map["Eccentricity"] = 19;
  name_to_col_map["Inclination"] = 20;
  name_to_col_map["RAAN"] = 21;
  name_to_col_map["Argument of Periapsis"] = 22;
  name_to_col_map["True Anomaly"] = 23;
  name_to_col_map["Orbital Rate"] = 24;
  name_to_col_map["Orbital Angular Acceleration"] = 25;
  name_to_col_map["Total Energy"] = 26;

  // first, open "pipe" to gnuplot
  FILE* gnuplot_pipe = popen("gnuplot", "w");
  // if it exists
  if (gnuplot_pipe) {
    fprintf(gnuplot_pipe,
            "set terminal png size 800,500 font ',14' linewidth 2\n");
    // formatting
    fprintf(gnuplot_pipe, "set output '../%s.png'\n", output_file_name.c_str());
    fprintf(gnuplot_pipe, "set offsets graph 0, 0, 0.05, 0.05\n");
    fprintf(gnuplot_pipe, "set xlabel 'Time [s]'\n");
    if (name_to_units_map.count(plotted_parameter) == 0) {
      throw std::invalid_argument("Unrecognized parameter to plot");
    }
    std::string plotted_element_units = name_to_units_map[plotted_parameter];
    int plotted_element_column = name_to_col_map[plotted_parameter];
    fprintf(gnuplot_pipe, "set ylabel '%s %s'\n",
              plotted_parameter.c_str(),
              plotted_element_units.c_str());
    fprintf(gnuplot_pipe, "set key top right\n");

    // plotting

    // first satellite
    std::string current_datafile_name = input_datafile_name_vector.at(0);
    if (input_datafile_name_vector.size() == 1) {
      fprintf(gnuplot_pipe,
                "plot '../data/%s' using 1:%d with lines lw 1 title '%s' \n",
                current_datafile_name.c_str(),
                plotted_element_column, 
                current_datafile_name.c_str());
    }

    else {
      fprintf(gnuplot_pipe,
                "plot '../data/%s' using 1:%d with lines lw 1 title '%s'\\\n",
                current_datafile_name.c_str(),
                plotted_element_column, 
                current_datafile_name.c_str());

    }

    for (size_t datafile_index = 1;
         datafile_index < input_datafile_name_vector.size(); datafile_index++) {
      current_datafile_name = input_datafile_name_vector.at(datafile_index);
      if (datafile_index < input_datafile_name_vector.size() - 1) {
          fprintf(gnuplot_pipe,
                    ",'../data/%s' using 1:%d with lines lw 1 title '%s'\\\n",
                    current_datafile_name.c_str(),
                    plotted_element_column, 
                    current_datafile_name.c_str());
      }

      else {
          fprintf(gnuplot_pipe,
                    ",'../data/%s' using 1:%d with lines lw 1 title '%s'\n",
                    current_datafile_name.c_str(),
                    plotted_element_column, 
                    current_datafile_name.c_str());
      }
    }

    fprintf(gnuplot_pipe, "exit \n");
    pclose(gnuplot_pipe);
    std::cout << "Done\n";

  } else {
    std::cout << "gnuplot not found";
  }

  return;
}


// Objective: simulate the input satellites over the specified total sim time,
// and visualize the resulting orbits in an interactive 3D plot using gnuplot
void sim_and_draw_orbit_gnuplot(std::vector<Satellite> input_satellite_vector,
                                const SimParameters& input_sim_parameters,
                                const std::string output_file_name) {
  if (input_satellite_vector.size() < 1) {
    std::cout << "No input Satellite objects\n";
    return;
  }
  // first, open "pipe" to gnuplot
  std::string gnuplot_arg_string = "gnuplot";
  if (input_sim_parameters.terminal_name_3D == "qt") {
    gnuplot_arg_string += " -persist";
  }
  FILE* gnuplot_pipe = popen(gnuplot_arg_string.c_str(), "w");

  // if it exists
  if (gnuplot_pipe) {
    fprintf(gnuplot_pipe, "set terminal '%s' size 900,700 font ',14'\n",
            input_sim_parameters.terminal_name_3D.c_str());
    if (input_sim_parameters.terminal_name_3D == "png") {
      fprintf(gnuplot_pipe, "set output '../%s.png'\n",
              output_file_name.c_str());
    }
    // formatting
    fprintf(gnuplot_pipe, "set xlabel 'x [m]' offset 0,-2\n");
    fprintf(gnuplot_pipe, "set ylabel 'y [m]' offset -2,0\n");
    fprintf(gnuplot_pipe, "set zlabel 'z [m]'\n");
    fprintf(gnuplot_pipe,
            "set title 'Simulated orbits up to time %.2f s' offset 0,-7.5\n",
            input_sim_parameters.total_sim_time);
    // fprintf(gnuplot_pipe,"set view 70,1,1,1\n");
    fprintf(gnuplot_pipe, "set view equal xyz\n");
    if (input_sim_parameters.x_increment != 0) {
      fprintf(gnuplot_pipe, "set xtics %e offset 0,-1.5\n",
              input_sim_parameters.x_increment);
    } else {
      fprintf(gnuplot_pipe, "set xtics offset 0,-1.5\n");
    }
    if (input_sim_parameters.y_increment != 0) {
      fprintf(gnuplot_pipe, "set ytics %e offset -1.5,0\n",
              input_sim_parameters.y_increment);
    } else {
      fprintf(gnuplot_pipe, "set ytics offset -1.5,0\n");
    }
    if (input_sim_parameters.z_increment != 0) {
      fprintf(gnuplot_pipe, "set ztics %e\n", input_sim_parameters.z_increment);
    }
    fprintf(gnuplot_pipe, "unset colorbox\n");
    fprintf(gnuplot_pipe, "set style fill transparent solid 1.0\n");

    fprintf(gnuplot_pipe, "set key offset 0,-10\n");
    fprintf(gnuplot_pipe, "set hidden3d front\n");

    // plotting
    // first let's set the stage for plotting the Earth
    fprintf(gnuplot_pipe, "R_Earth=%.17g\n", radius_Earth);
    fprintf(gnuplot_pipe, "set isosamples 50,50\n");
    fprintf(gnuplot_pipe, "set parametric\n");
    fprintf(gnuplot_pipe, "set urange [-pi/2:pi/2]\n");
    fprintf(gnuplot_pipe, "set vrange [0:2*pi]\n");

    // first satellite
    Satellite current_satellite = input_satellite_vector.at(0);
    if (input_satellite_vector.size() == 1) {
      if (current_satellite.plotting_color_.size() > 0) {
        fprintf(gnuplot_pipe,
                "splot '-' with lines lw 1 lc rgb '%s' title '%s' \\\n",
                current_satellite.plotting_color_.c_str(),
                current_satellite.get_name().c_str());
      } else {
        fprintf(gnuplot_pipe, "splot '-' with lines lw 1 title '%s' \\\n",
                current_satellite.get_name().c_str());
      }

    }

    else {
      if (current_satellite.plotting_color_.size() > 0) {
        fprintf(gnuplot_pipe,
                "splot '-' with lines lw 1 lc rgb '%s' title '%s'\\\n",
                current_satellite.plotting_color_.c_str(),
                current_satellite.get_name().c_str());
      } else {
        fprintf(gnuplot_pipe, "splot '-' with lines lw 1 title '%s'\\\n",
                current_satellite.get_name().c_str());
      }
    }

    for (size_t satellite_index = 1;
         satellite_index < input_satellite_vector.size(); satellite_index++) {
      current_satellite = input_satellite_vector.at(satellite_index);
      if (satellite_index < input_satellite_vector.size() - 1) {
        if (current_satellite.plotting_color_.size() > 0) {
          fprintf(gnuplot_pipe,
                  ",'-' with lines lw 1 lc rgb '%s' title '%s' \\\n",
                  current_satellite.plotting_color_.c_str(),
                  current_satellite.get_name().c_str());
        } else {
          fprintf(gnuplot_pipe, ",'-' with lines lw 1 title '%s' \\\n",
                  current_satellite.get_name().c_str());
        }

      }

      else {
        if (current_satellite.plotting_color_.size() > 0) {
          fprintf(gnuplot_pipe,
                  ",'-' with lines lw 1 lc rgb '%s' title '%s'\\\n",
                  current_satellite.plotting_color_.c_str(),
                  current_satellite.get_name().c_str());
        } else {
          fprintf(gnuplot_pipe, ",'-' with lines lw 1 title '%s'\\\n",
                  current_satellite.get_name().c_str());
        }
      }
    }
    fprintf(gnuplot_pipe,
            ",R_Earth*cos(u)*cos(v),R_Earth*cos(u)*sin(v),R_Earth*sin(u) "
            "notitle with pm3d fillcolor rgbcolor 'navy'\n");

    // now the orbit data, inline, one satellite at a time
    std::cout << "Running simulation...\n";
    for (size_t satellite_index = 0;
         satellite_index < input_satellite_vector.size(); satellite_index++) {
      Satellite current_satellite = input_satellite_vector.at(satellite_index);
      std::array<double, 3> initial_position =
          current_satellite.get_ECI_position();
      fprintf(gnuplot_pipe, "%.17g %.17g %.17g\n", initial_position.at(0),
              initial_position.at(1), initial_position.at(2));

      std::array<double, 3> evolved_position = {};

      double timestep_to_use = input_sim_parameters.initial_timestep_guess;
      double current_satellite_time =
          current_satellite.get_instantaneous_time();
      while (current_satellite_time < input_sim_parameters.total_sim_time) {
        std::pair<double, double> drag_elements = {input_sim_parameters.F_10,
                                                   input_sim_parameters.A_p};
        std::pair<double, int> new_timestep_and_error_code =
            current_satellite.evolve_RK45(
                input_sim_parameters.epsilon, timestep_to_use,
                input_sim_parameters.perturbation_bool,
                input_sim_parameters.drag_bool, drag_elements);
        double new_timestep = new_timestep_and_error_code.first;
        int error_code = new_timestep_and_error_code.second;
        if (error_code != 0) {
          std::cout << "Error code " << error_code
                    << " detected, halting visualization\n";
          fprintf(gnuplot_pipe, "e\n");
          fprintf(gnuplot_pipe, "exit \n");
          pclose(gnuplot_pipe);
          return;
        }
        timestep_to_use = new_timestep;
        evolved_position = current_satellite.get_ECI_position();
        current_satellite_time = current_satellite.get_instantaneous_time();
        fprintf(gnuplot_pipe, "%.17g %.17g %.17g\n", evolved_position.at(0),
                evolved_position.at(1), evolved_position.at(2));
      }
      fprintf(gnuplot_pipe, "e\n");
    }

    if (input_sim_parameters.terminal_name_3D == "qt") {
      fprintf(gnuplot_pipe, "pause mouse keypress\n");
    }
    fprintf(gnuplot_pipe, "exit \n");
    pclose(gnuplot_pipe);
    std::cout << "Done\n";

  } else {
    std::cout << "gnuplot not found";
  }

  return;
}

// Objective: simulate the input satellites over the specified total sim time,
// and plot a specific orbital element over time
void sim_and_plot_orbital_elem_gnuplot(
    std::vector<Satellite> input_satellite_vector,
    const SimParameters& input_sim_parameters,
    std::string input_orbital_element_name, const std::string file_name) {
  if (input_satellite_vector.size() < 1) {
    std::cout << "No input Satellite objects\n";
    return;
  }

  // first, open "pipe" to gnuplot
  FILE* gnuplot_pipe = popen("gnuplot", "w");
  // if it exists
  if (gnuplot_pipe) {
    fprintf(gnuplot_pipe,
            "set terminal png size 800,500 font ',14' linewidth 2\n");
    // formatting
    fprintf(gnuplot_pipe, "set output '../%s.png'\n", file_name.c_str());
    fprintf(gnuplot_pipe, "set offsets graph 0, 0, 0.05, 0.05\n");
    fprintf(gnuplot_pipe, "set xlabel 'Time [s]'\n");
    if (input_orbital_element_name == "Semimajor Axis") {
      fprintf(gnuplot_pipe, "set ylabel '%s [m]'\n",
              input_orbital_element_name.c_str());
    } else if (input_orbital_element_name == "Eccentricity") {
      fprintf(gnuplot_pipe, "set ylabel '%s'\n",
              input_orbital_element_name.c_str());
    } else {
      fprintf(gnuplot_pipe, "set ylabel '%s [deg]'\n",
              input_orbital_element_name.c_str());
    }
    fprintf(gnuplot_pipe, "set title '%s simulated up to time %.2f s'\n",
            input_orbital_element_name.c_str(),
            input_sim_parameters.total_sim_time);
    fprintf(gnuplot_pipe, "set key top right\n");

    // plotting

    // first satellite
    Satellite current_satellite = input_satellite_vector.at(0);
    if (input_satellite_vector.size() == 1) {
      if (current_satellite.plotting_color_.size() > 0) {
        fprintf(gnuplot_pipe,
                "plot '-' using 1:2 with lines lw 1 lc rgb '%s' title '%s' \n",
                current_satellite.plotting_color_.c_str(),
                current_satellite.get_name().c_str());
      } else {
        fprintf(gnuplot_pipe,
                "plot '-' using 1:2 with lines lw 1 title '%s' \n",
                current_satellite.get_name().c_str());
      }

    }

    else {
      if (current_satellite.plotting_color_.size() > 0) {
        fprintf(gnuplot_pipe,
                "plot '-' using 1:2 with lines lw 1 lc rgb '%s' title '%s'\\\n",
                current_satellite.plotting_color_.c_str(),
                current_satellite.get_name().c_str());
      } else {
        fprintf(gnuplot_pipe,
                "plot '-' using 1:2 with lines lw 1 title '%s'\\\n",
                current_satellite.get_name().c_str());
      }
    }

    for (size_t satellite_index = 1;
         satellite_index < input_satellite_vector.size(); satellite_index++) {
      current_satellite = input_satellite_vector.at(satellite_index);
      if (satellite_index < input_satellite_vector.size() - 1) {
        if (current_satellite.plotting_color_.size() > 0) {
          fprintf(gnuplot_pipe,
                  ",'-' using 1:2 with lines lw 1 lc rgb '%s' title '%s' \\\n",
                  current_satellite.plotting_color_.c_str(),
                  current_satellite.get_name().c_str());
        } else {
          fprintf(gnuplot_pipe,
                  ",'-' using 1:2 with lines lw 1 title '%s' \\\n",
                  current_satellite.get_name().c_str());
        }

      }

      else {
        if (current_satellite.plotting_color_.size() > 0) {
          fprintf(gnuplot_pipe,
                  ",'-' using 1:2 with lines lw 1 lc rgb '%s' title '%s'\n",
                  current_satellite.plotting_color_.c_str(),
                  current_satellite.get_name().c_str());
        } else {
          fprintf(gnuplot_pipe, ",'-' using 1:2 with lines lw 1 title '%s'\n",
                  current_satellite.get_name().c_str());
        }
      }
    }

    // now the orbit data, inline, one satellite at a time
    std::cout << "Running simulation...\n";
    for (size_t satellite_index = 0;
         satellite_index < input_satellite_vector.size(); satellite_index++) {
      Satellite current_satellite = input_satellite_vector.at(satellite_index);
      double val =
          current_satellite.get_orbital_parameter(input_orbital_element_name);
      double current_satellite_time =
          current_satellite.get_instantaneous_time();
      fprintf(gnuplot_pipe, "%.17g %.17g\n", current_satellite_time, val);

      double evolved_val = {0};

      double timestep_to_use = input_sim_parameters.initial_timestep_guess;
      current_satellite_time = current_satellite.get_instantaneous_time();
      while (current_satellite_time < input_sim_parameters.total_sim_time) {
        // std::cout <<
        // "========================================================\n";
        // std::cout << "Running an evolve step at satellite time " <<
        // current_satellite_time << "\n";
        std::pair<double, double> drag_elements = {input_sim_parameters.F_10,
                                                   input_sim_parameters.A_p};
        std::pair<double, int> new_timestep_and_error_code =
            current_satellite.evolve_RK45(
                input_sim_parameters.epsilon, timestep_to_use,
                input_sim_parameters.perturbation_bool,
                input_sim_parameters.drag_bool, drag_elements);
        double new_timestep = new_timestep_and_error_code.first;
        int error_code = new_timestep_and_error_code.second;

        if (error_code != 0) {
          std::cout << "Error code " << error_code
                    << " detected, halting visualization\n";
          fprintf(gnuplot_pipe, "e\n");
          fprintf(gnuplot_pipe, "exit \n");
          pclose(gnuplot_pipe);
          return;
        }
        timestep_to_use = new_timestep;
        evolved_val =
            current_satellite.get_orbital_parameter(input_orbital_element_name);
        current_satellite_time = current_satellite.get_instantaneous_time();
        // printf("%.17g %.17g\n", current_satellite_time,
        //   evolved_val);
        fprintf(gnuplot_pipe, "%.17g %.17g\n", current_satellite_time,
                evolved_val);
      }
      fprintf(gnuplot_pipe, "e\n");
    }

    fprintf(gnuplot_pipe, "exit \n");
    pclose(gnuplot_pipe);
    std::cout << "Done\n";

  } else {
    std::cout << "gnuplot not found";
  }

  return;
}

Matrix3d z_rot_matrix(const double input_angle) {
  Matrix3d z_rotation_matrix;
  z_rotation_matrix << cos(input_angle), -sin(input_angle), 0, sin(input_angle),
      cos(input_angle), 0, 0, 0, 1;

  return z_rotation_matrix;
}

Matrix3d y_rot_matrix(const double input_angle) {
  Matrix3d y_rotation_matrix;
  y_rotation_matrix << cos(input_angle), 0, sin(input_angle), 0, 1, 0,
      -sin(input_angle), 0, cos(input_angle);

  return y_rotation_matrix;
}

Matrix3d x_rot_matrix(const double input_angle) {
  Matrix3d x_rotation_matrix;
  x_rotation_matrix << 1, 0, 0, 0, cos(input_angle), -sin(input_angle), 0,
      sin(input_angle), cos(input_angle);

  return x_rotation_matrix;
}

// Objective: simulate the input satellites over the specified total sim time,
// and plot a specific attitude-related value over time
void sim_and_plot_attitude_evolution_gnuplot(
    std::vector<Satellite> input_satellite_vector,
    const SimParameters& input_sim_parameters,
    const std::string input_plotted_val_name, const std::string file_name) {
  if (input_satellite_vector.size() < 1) {
    std::cout << "No input Satellite objects\n";
    return;
  }

  // first, open "pipe" to gnuplot
  FILE* gnuplot_pipe = popen("gnuplot", "w");
  // if it exists
  if (gnuplot_pipe) {
    fprintf(gnuplot_pipe,
            "set terminal png size 800,500 font ',14' linewidth 2\n");
    // formatting
    fprintf(gnuplot_pipe, "set output '../%s.png'\n", file_name.c_str());
    fprintf(gnuplot_pipe, "set offsets graph 0, 0, 0.05, 0.05\n");
    fprintf(gnuplot_pipe, "set xlabel 'Time [s]'\n");
    if ((input_plotted_val_name == "omega_x") ||
        (input_plotted_val_name == "omega_y") ||
        (input_plotted_val_name == "omega_z")) {
      fprintf(gnuplot_pipe, "set ylabel '%s [deg/s]'\n",
              input_plotted_val_name.c_str());
    } else if ((input_plotted_val_name == "Roll") ||
               (input_plotted_val_name == "Pitch") ||
               (input_plotted_val_name == "Yaw")) {
      fprintf(gnuplot_pipe, "set ylabel '%s [deg]'\n",
              input_plotted_val_name.c_str());
    } else {
      fprintf(gnuplot_pipe, "set ylabel '%s'\n",
              input_plotted_val_name.c_str());
    }
    fprintf(gnuplot_pipe, "set title '%s simulated up to time %.2f s'\n",
            input_plotted_val_name.c_str(),
            input_sim_parameters.total_sim_time);
    fprintf(gnuplot_pipe, "set key top right\n");

    // plotting

    // first satellite
    Satellite current_satellite = input_satellite_vector.at(0);
    if (input_satellite_vector.size() == 1) {
      if (current_satellite.plotting_color_.size() > 0) {
        fprintf(gnuplot_pipe,
                "plot '-' using 1:2 with lines lw 1 lc rgb '%s' title '%s' \n",
                current_satellite.plotting_color_.c_str(),
                current_satellite.get_name().c_str());
      } else {
        fprintf(gnuplot_pipe,
                "plot '-' using 1:2 with lines lw 1 title '%s' \n",
                current_satellite.get_name().c_str());
      }

    }

    else {
      if (current_satellite.plotting_color_.size() > 0) {
        fprintf(gnuplot_pipe,
                "plot '-' using 1:2 with lines lw 1 lc rgb '%s' title '%s'\\\n",
                current_satellite.plotting_color_.c_str(),
                current_satellite.get_name().c_str());
      } else {
        fprintf(gnuplot_pipe,
                "plot '-' using 1:2 with lines lw 1 title '%s'\\\n",
                current_satellite.get_name().c_str());
      }
    }

    for (size_t satellite_index = 1;
         satellite_index < input_satellite_vector.size(); satellite_index++) {
      current_satellite = input_satellite_vector.at(satellite_index);
      if (satellite_index < input_satellite_vector.size() - 1) {
        if (current_satellite.plotting_color_.size() > 0) {
          fprintf(gnuplot_pipe,
                  ",'-' using 1:2 with lines lw 1 lc rgb '%s' title '%s' \\\n",
                  current_satellite.plotting_color_.c_str(),
                  current_satellite.get_name().c_str());
        } else {
          fprintf(gnuplot_pipe,
                  ",'-' using 1:2 with lines lw 1 title '%s' \\\n",
                  current_satellite.get_name().c_str());
        }

      }

      else {
        if (current_satellite.plotting_color_.size() > 0) {
          fprintf(gnuplot_pipe,
                  ",'-' using 1:2 with lines lw 1 lc rgb '%s' title '%s'\n",
                  current_satellite.plotting_color_.c_str(),
                  current_satellite.get_name().c_str());
        } else {
          fprintf(gnuplot_pipe, ",'-' using 1:2 with lines lw 1 title '%s'\n",
                  current_satellite.get_name().c_str());
        }
      }
    }

    // now the orbit data, inline, one satellite at a time
    std::cout << "Running simulation...\n";
    for (size_t satellite_index = 0;
         satellite_index < input_satellite_vector.size(); satellite_index++) {
      Satellite current_satellite = input_satellite_vector.at(satellite_index);
      double val = current_satellite.get_attitude_val(input_plotted_val_name);
      double current_satellite_time =
          current_satellite.get_instantaneous_time();
      fprintf(gnuplot_pipe, "%.17g %.17g\n", current_satellite_time, val);

      double evolved_val = {0};

      double timestep_to_use = input_sim_parameters.initial_timestep_guess;
      current_satellite_time = current_satellite.get_instantaneous_time();
      while (current_satellite_time < input_sim_parameters.total_sim_time) {
        std::pair<double, double> drag_elements = {input_sim_parameters.F_10,
                                                   input_sim_parameters.A_p};
        std::pair<double, int> new_timestep_and_error_code =
            current_satellite.evolve_RK45(
                input_sim_parameters.epsilon, timestep_to_use,
                input_sim_parameters.perturbation_bool,
                input_sim_parameters.drag_bool, drag_elements);
        double new_timestep = new_timestep_and_error_code.first;
        int error_code = new_timestep_and_error_code.second;

        if (error_code != 0) {
          std::cout << "Error code " << error_code
                    << " detected, halting visualization\n";
          fprintf(gnuplot_pipe, "e\n");
          fprintf(gnuplot_pipe, "exit \n");
          pclose(gnuplot_pipe);
          return;
        }
        timestep_to_use = new_timestep;
        evolved_val =
            current_satellite.get_attitude_val(input_plotted_val_name);
        current_satellite_time = current_satellite.get_instantaneous_time();
        fprintf(gnuplot_pipe, "%.17g %.17g\n", current_satellite_time,
                evolved_val);
      }
      fprintf(gnuplot_pipe, "e\n");
    }

    fprintf(gnuplot_pipe, "exit \n");
    pclose(gnuplot_pipe);
    std::cout << "Done\n";

  } else {
    std::cout << "gnuplot not found";
  }

  return;
}

Matrix3d rollyawpitch_bodyframe_to_LVLH_matrix(const double input_roll,
                                               const double input_pitch,
                                               const double input_yaw) {
  // Going off convention in
  // https://ntrs.nasa.gov/api/citations/19770024112/downloads/19770024112.pdf
  // which appears to be first rotating by pitch, then yaw, then roll
  //  This is a convention where the rotations are performed as: v_body =
  //  R_x(roll)R_z(yaw)R_y(pitch)v_LVLH That's an intrinsic Tait-Bryan xzy
  //  sequence, following convention described in
  //  https://en.wikipedia.org/wiki/Rotation_matrix To my understanding, this is
  //  equivalent to an extrinsic Tait-Bryan yzx sequence

  // First going to construct LVLH_to_bodyframe matrix, then take transpose to
  // make it the bodyframe_to_LVLH matrix
  Matrix3d bodyframe_to_LVLH_matrix = x_rot_matrix(input_roll) *
                                      z_rot_matrix(input_yaw) *
                                      y_rot_matrix(input_pitch);
  bodyframe_to_LVLH_matrix.transpose();

  return bodyframe_to_LVLH_matrix;
}

std::array<double, 4> rollyawpitch_angles_to_quaternion(
    const double input_roll, const double input_pitch, const double input_yaw) {
  // Refs :
  // https://www.euclideanspace.com/maths/geometry/rotations/conversions/eulerToQuaternion/Euler%20to%20quat.pdf
  //  and
  //  https://ntrs.nasa.gov/api/citations/19770024112/downloads/19770024112.pdf
  Vector4d q_pitch = {cos(input_pitch / 2), 0.0, sin(input_pitch / 2), 0.0};
  Vector4d q_yaw = {cos(input_yaw / 2), 0.0, 0.0, sin(input_yaw / 2)};
  Vector4d q_roll = {cos(input_roll / 2), sin(input_roll / 2), 0.0, 0.0};
  Vector4d q_tot = quaternion_multiplication(
      q_roll, quaternion_multiplication(q_yaw, q_pitch));
  std::array<double, 4> output_quaternion = {q_tot(0), q_tot(1), q_tot(2),
                                             q_tot(3)};
  return output_quaternion;
}

Matrix3d LVLH_to_body_transformation_matrix_from_quaternion(
    const std::array<double, 4> input_bodyframe_quaternion_relative_to_LVLH) {
  // Ref:
  // https://ntrs.nasa.gov/api/citations/20240009554/downloads/Space%20Attitude%20Development%20Control.pdf
  // section 4.3.1
  double q0 = input_bodyframe_quaternion_relative_to_LVLH.at(0);
  double q1 = input_bodyframe_quaternion_relative_to_LVLH.at(1);
  double q2 = input_bodyframe_quaternion_relative_to_LVLH.at(2);
  double q3 = input_bodyframe_quaternion_relative_to_LVLH.at(3);

  Matrix3d LVLH_to_body_mat;
  LVLH_to_body_mat << 2 * q0 * q0 - 1 + 2 * q1 * q1, 2 * q1 * q2 + 2 * q0 * q3,
      2 * q1 * q3 - 2 * q0 * q2, 2 * q1 * q2 - 2 * q0 * q3,
      2 * q0 * q0 - 1 + 2 * q2 * q2, 2 * q2 * q3 + 2 * q0 * q1,
      2 * q1 * q3 + 2 * q0 * q2, 2 * q2 * q3 - 2 * q0 * q1,
      2 * q0 * q0 - 1 + 2 * q3 * q3;

  return LVLH_to_body_mat;
}

// Objective: compute time derivatives of bodyframe angular velocities
std::array<double, 3> calculate_spacecraft_bodyframe_angular_acceleration(
    const Matrix3d J_matrix,
    const std::vector<BodyframeTorqueProfile>
        input_bodyframe_torque_profile_list,
    const Vector3d input_omega_I,
    const double input_orbital_angular_acceleration,
    const Vector3d input_omega_bodyframe_wrt_LVLH_in_body_frame,
    const Matrix3d input_LVLH_to_bodyframe_transformation_matrix,
    const Vector3d input_omega_LVLH_wrt_inertial_in_LVLH,
    const double input_evaluation_time) {
  // Using approach from
  // https://ntrs.nasa.gov/api/citations/20240009554/downloads/Space%20Attitude%20Development%20Control.pdf
  // , Ch. 4 especially Objective: calculate omega_dot in spacecraft body frame
  // w.r.t. the LVLH frame, expressed in the spacecraft body frame Disturbance
  // torques, if eventually added in, should just be more torque profiles in the
  // satellite object's list of torque profiles So the input_torques vector
  // includes both disturbance and control torques
  Vector3d bodyframe_torque_vec = {0, 0, 0};

  for (const BodyframeTorqueProfile bodyframe_torque_profile :
       input_bodyframe_torque_profile_list) {
    if ((input_evaluation_time >= bodyframe_torque_profile.t_start_) &&
        (input_evaluation_time <= bodyframe_torque_profile.t_end_)) {
      for (size_t ind = 0; ind < 3; ind++) {
        bodyframe_torque_vec(ind) +=
            bodyframe_torque_profile.bodyframe_torque_list.at(ind);
      }
    }
  }

  Vector3d omega_lvlh_dot = {0, -input_orbital_angular_acceleration, 0};
  Matrix3d inverted_J_matrix = J_matrix;
  inverted_J_matrix.inverse();

  Vector3d comp1 =
      -inverted_J_matrix * (input_omega_I.cross(J_matrix * input_omega_I));
  Vector3d comp2 = inverted_J_matrix * bodyframe_torque_vec;
  Vector3d comp3 = input_omega_bodyframe_wrt_LVLH_in_body_frame.cross(
      input_LVLH_to_bodyframe_transformation_matrix *
      input_omega_LVLH_wrt_inertial_in_LVLH);
  Vector3d omega_dot_vec_LVLH_wrt_inertial_in_ECI = {
      0, -input_orbital_angular_acceleration, 0};
  Vector3d comp4 = -input_LVLH_to_bodyframe_transformation_matrix *
                   omega_dot_vec_LVLH_wrt_inertial_in_ECI;

  Vector3d angular_acceleration_bodyframe_wrt_LVLH_in_bodyframe =
      comp1 + comp2 + comp3 + comp4;
  std::array<double, 3>
      angular_acceleration_bodyframe_wrt_LVLH_in_bodyframe_array = {
          angular_acceleration_bodyframe_wrt_LVLH_in_bodyframe(0),
          angular_acceleration_bodyframe_wrt_LVLH_in_bodyframe(1),
          angular_acceleration_bodyframe_wrt_LVLH_in_bodyframe(2)};

  return angular_acceleration_bodyframe_wrt_LVLH_in_bodyframe_array;
}

// Objective: compute time derivatives of components of quaternion describing
// orientation of satellite body frame with respect to LVLH frame
Vector4d quaternion_kinematics_equation(
    const Vector4d quaternion_of_bodyframe_relative_to_ref_frame,
    const Vector3d angular_velocity_vec_wrt_ref_frame_in_body_frame) {
  // Ref:
  // https://ntrs.nasa.gov/api/citations/20240009554/downloads/Space%20Attitude%20Development%20Control.pdf
  // Section 4.1.2

  // Here assuming scalar component is first in quaternion
  Vector3d vec_component_of_quaternion = {
      quaternion_of_bodyframe_relative_to_ref_frame(1),
      quaternion_of_bodyframe_relative_to_ref_frame(2),
      quaternion_of_bodyframe_relative_to_ref_frame(3)};
  double q_0_dot =
      (-0.5) * angular_velocity_vec_wrt_ref_frame_in_body_frame.dot(
                   vec_component_of_quaternion);
  Vector3d vec_q_dot =
      (-0.5) * angular_velocity_vec_wrt_ref_frame_in_body_frame.cross(
                   vec_component_of_quaternion) +
      (0.5) * quaternion_of_bodyframe_relative_to_ref_frame(0) *
          angular_velocity_vec_wrt_ref_frame_in_body_frame;

  Vector4d quaternion_derivative;
  quaternion_derivative << q_0_dot, vec_q_dot(0), vec_q_dot(1), vec_q_dot(2);

  return quaternion_derivative;
}

// Computes time derivative of combined satellite quaternion + angular velocity
// vector
std::array<double, 7> RK45_satellite_body_angular_deriv_function(
    const std::array<double, 7> combined_bodyframe_angular_array,
    const Matrix3d J_matrix,
    const std::vector<BodyframeTorqueProfile>
        input_bodyframe_torque_profile_list,
    const Vector3d input_omega_I, double input_orbital_angular_acceleration,
    const Matrix3d input_LVLH_to_bodyframe_transformation_matrix,
    const Vector3d input_omega_LVLH_wrt_inertial_in_LVLH,
    const double input_evaluation_time) {
  // Setting this up for an RK45 step with
  // y={q_0,q_1,q_2,q_3,omega_1,omega_2,omega_3} Objective is to produce dy/dt
  // Input quaternion should be quaternion of bodyframe relative to LVLH
  std::array<double, 7> combined_angular_derivative_array = {0, 0, 0, 0,
                                                             0, 0, 0};
  Vector4d quaternion = {combined_bodyframe_angular_array.at(0),
                         combined_bodyframe_angular_array.at(1),
                         combined_bodyframe_angular_array.at(2),
                         combined_bodyframe_angular_array.at(3)};
  Vector3d body_angular_velocity_vec_wrt_LVLH_in_body_frame = {
      combined_bodyframe_angular_array.at(4),
      combined_bodyframe_angular_array.at(5),
      combined_bodyframe_angular_array.at(6)};
  Vector4d quaternion_derivative = quaternion_kinematics_equation(
      quaternion, body_angular_velocity_vec_wrt_LVLH_in_body_frame);
  std::array<double, 3> body_angular_acceleration_vec_wrt_LVLH_in_body_frame =
      calculate_spacecraft_bodyframe_angular_acceleration(
          J_matrix, input_bodyframe_torque_profile_list, input_omega_I,
          input_orbital_angular_acceleration,
          body_angular_velocity_vec_wrt_LVLH_in_body_frame,
          input_LVLH_to_bodyframe_transformation_matrix,
          input_omega_LVLH_wrt_inertial_in_LVLH, input_evaluation_time);

  for (size_t ind = 0; ind < 4; ind++) {
    combined_angular_derivative_array.at(ind) = quaternion_derivative(ind);
  }
  for (size_t ind = 4; ind < 7; ind++) {
    combined_angular_derivative_array.at(ind) =
        body_angular_acceleration_vec_wrt_LVLH_in_body_frame.at(ind - 4);
  }

  return combined_angular_derivative_array;
}

// Calculate angular velocity of spacecraft body frame with respect to
// inertial (ECI) frame
Vector3d calculate_omega_I(
    const Vector3d input_bodyframe_ang_vel_vector_wrt_lvlh,
    const Matrix3d input_LVLH_to_bodyframe_transformation_matrix,
    const double input_orbital_rate) {
  // Eq. 4.17 in
  // https://ntrs.nasa.gov/api/citations/20240009554/downloads/Space%20Attitude%20Development%20Control.pdf
  Vector3d omega_LVLH_wrt_inertial_in_LVLH = {0, -input_orbital_rate, 0};
  Vector3d omega_I = input_bodyframe_ang_vel_vector_wrt_lvlh +
                     input_LVLH_to_bodyframe_transformation_matrix *
                         omega_LVLH_wrt_inertial_in_LVLH;
  return omega_I;
}

Matrix3d construct_J_matrix(const double input_Jxx, const double input_Jyy,
                            const double input_Jzz) {
  Matrix3d J_matrix = Matrix3d::Zero();
  J_matrix(0, 0) = input_Jxx;
  J_matrix(1, 1) = input_Jyy;
  J_matrix(2, 2) = input_Jzz;
  return J_matrix;
}

std::array<double, 13>
RK45_combined_orbit_position_velocity_attitude_deriv_function(
    const std::array<double, 13>
        combined_position_velocity_bodyframe_angular_array,
    const Matrix3d J_matrix,
    const std::vector<BodyframeTorqueProfile>
        input_bodyframe_torque_profile_list,
    const Vector3d input_omega_I, double input_orbital_angular_acceleration,
    const Matrix3d input_LVLH_to_bodyframe_transformation_matrix,
    const Vector3d input_omega_LVLH_wrt_inertial_in_LVLH,
    const double input_spacecraft_mass,
    const std::vector<ThrustProfileLVLH> input_list_of_thrust_profiles_LVLH,
    const double input_evaluation_time, const double input_inclination,
    const double input_arg_of_periapsis, const double input_true_anomaly,
    const double input_F_10, const double input_A_p, const double input_A_s,
    const bool perturbation, const bool atmospheric_drag) {
  // Input vector is in the form of {ECI_position,
  // ECI_velocity,bodyframe_quaternion_to_LVLH,bodyframe_omega_wrt_LVLH}
  // Objective is to output derivative of that vector
  std::array<double, 13> derivative_vec;
  std::array<double, 6> combined_position_and_velocity_array;
  std::array<double, 7> combined_angular_array;

  for (size_t ind = 0; ind < combined_position_and_velocity_array.size();
       ind++) {
    combined_position_and_velocity_array.at(ind) =
        combined_position_velocity_bodyframe_angular_array.at(ind);
  }
  for (size_t ind = 0; ind < combined_angular_array.size(); ind++) {
    combined_angular_array.at(ind) =
        combined_position_velocity_bodyframe_angular_array.at(
            ind + combined_position_and_velocity_array.size());
  }
  std::array<double, 6> orbital_position_and_velocity_derivative_array =
      RK45_deriv_function_orbit_position_and_velocity(
          combined_position_and_velocity_array, input_spacecraft_mass,
          input_list_of_thrust_profiles_LVLH, input_evaluation_time,
          input_inclination, input_arg_of_periapsis, input_true_anomaly,
          input_F_10, input_A_p, input_A_s, input_spacecraft_mass, perturbation,
          atmospheric_drag);

  std::array<double, 7> angular_derivative_array =
      RK45_satellite_body_angular_deriv_function(
          combined_angular_array, J_matrix, input_bodyframe_torque_profile_list,
          input_omega_I, input_orbital_angular_acceleration,
          input_LVLH_to_bodyframe_transformation_matrix,
          input_omega_LVLH_wrt_inertial_in_LVLH, input_evaluation_time);

  for (size_t ind = 0;
       ind < orbital_position_and_velocity_derivative_array.size(); ind++) {
    derivative_vec.at(ind) =
        orbital_position_and_velocity_derivative_array.at(ind);
  }
  for (size_t ind = 0; ind < angular_derivative_array.size(); ind++) {
    derivative_vec.at(ind +
                      orbital_position_and_velocity_derivative_array.size()) =
        angular_derivative_array.at(ind);
  }

  return derivative_vec;
}

std::array<double, 4> normalize_quaternion(
    std::array<double, 4> input_quaternion) {
  double length = 0;
  for (size_t ind = 0; ind < input_quaternion.size(); ind++) {
    length += (input_quaternion.at(ind) * input_quaternion.at(ind));
  }
  for (size_t ind = 0; ind < input_quaternion.size(); ind++) {
    input_quaternion.at(ind) /= sqrt(length);
  }
  return input_quaternion;
}

std::array<double, 3> convert_quaternion_to_roll_yaw_pitch_angles(
    const std::array<double, 4> input_quaternion) {
  // Based on approach in
  // https://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToEuler/quat_2_euler_paper_ver2-1.pdf
  // Again, here I'm going off the pitch-yaw-roll sequence used in
  // https://ntrs.nasa.gov/api/citations/19770024112/downloads/19770024112.pdf
  // Where the angles which defined the orbiter attitude with respect to LVLH
  // are [A]_from_LVLH_to_BY = Rx(roll)*R_z(yaw)*R_y(pitch)

  Vector3d v3 = {0, 1,
                 0};  // unit vector in direction of last rotation performed,
                      // which in this sequence is the y unit vector
  Quaterniond eigen_quaternion(input_quaternion.at(0), input_quaternion.at(1),
                               input_quaternion.at(2), input_quaternion.at(3));
  eigen_quaternion.normalize();  // Just in case
  Matrix3d rotation_matrix = eigen_quaternion.toRotationMatrix();
  Vector3d euler_angles = rotation_matrix.eulerAngles(
      0, 2,
      1);  // Eigen uses intrinsic convention, so this is an intrinsic XZY
           // (denoted x-z'-y'' in the notation Wikipedia uses) or extrinsic YZX
           // sequence corresponding to R_x R_z R_y (e.g., see
           // https://arg.usask.ca/docs/skplatform/appendices/rotation_matrices.html)

  double roll = euler_angles(0);
  double yaw = euler_angles(1);
  double pitch = euler_angles(2);

  std::array<double, 3> output_angle_vec = {roll, yaw, pitch};
  return output_angle_vec;
}

std::array<double, 3> convert_array_from_LVLH_to_bodyframe(
    const std::array<double, 3> input_LVLH_frame_array, const double input_roll,
    const double input_yaw, const double input_pitch) {
  Matrix3d transformation_matrix =
      rollyawpitch_bodyframe_to_LVLH_matrix(input_roll, input_pitch, input_yaw);
  Vector3d input_LVLH_frame_vector = {input_LVLH_frame_array.at(0),
                                      input_LVLH_frame_array.at(1),
                                      input_LVLH_frame_array.at(2)};
  Vector3d bodyframe_vec = transformation_matrix * input_LVLH_frame_vector;
  std::array<double, 3> bodyframe_arr = {bodyframe_vec(0), bodyframe_vec(1),
                                         bodyframe_vec(2)};
  return bodyframe_arr;
}

Vector3d convert_lat_long_to_ECEF(const double latitude, const double longitude,
                                  const double height) {
  double a = 6378137;  // Equatorial radius [m]
  double b = 6356752;  // Polar radius [m]

  // Refs:
  // https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_geodetic_to_ECEF_coordinates
  // https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
  double e = 1 - ((b * b) / (a * a));
  double N =
      a / (sqrt(1 - ((e * e) / (1 + pow(cos(latitude) / sin(latitude), 2)))));
  double x = (N + height) * cos(latitude) * cos(longitude);
  double y = (N + height) * cos(latitude) * sin(longitude);
  double z = ((1 - e * e) * N + height) * sin(latitude);
  Vector3d output_ECEF_array = {x, y, z};
  return output_ECEF_array;
}

Vector3d convert_ECEF_to_ECI(const Vector3d input_ECEF_position,
                             const double input_time) {
  // Simple rotation-based conversion, not yet accounting for higher-fidelity
  // effects like changes to Earth axes Ref:
  // https://x-lumin.com/wp-content/uploads/2020/09/Coordinate_Transforms.pdf
  // https://space.stackexchange.com/questions/43187/is-this-commonly-attributed-eci-earth-centered-inertial-to-ecef-earth-centere
  // https://space.stackexchange.com/questions/38807/transform-eci-to-ecef
  double omega_Earth =
      0.261799387799149 * (1.0 / 3600.0);      // [radians / second]
  double ERA_at_J200 = 280.46 * (M_PI / 180);  // radians
  // Going to be assuming, for simplicity, that satellites start orbiting at the
  // J200 epoch
  double theta_g = ERA_at_J200 + omega_Earth * input_time;
  Matrix3d rotation_matrix = z_rot_matrix(theta_g);
  Vector3d ECI_array = rotation_matrix * input_ECEF_position;
  return ECI_array;
}

void sim_and_plot_gs_connectivity_distance_gnuplot(
    PhasedArrayGroundStation input_ground_station,
    std::vector<Satellite> input_satellite_vector,
    const SimParameters& input_sim_parameters, const std::string file_name) {
  // Objective: given an input ground station and satellite vector,
  // plot contact distances between ground station and satellites over time

  if (input_satellite_vector.size() < 1) {
    std::cout << "No input Satellite objects\n";
    return;
  }

  // first, open "pipe" to gnuplot
  FILE* gnuplot_pipe = popen("gnuplot", "w");
  // if it exists
  if (gnuplot_pipe) {
    fprintf(gnuplot_pipe,
            "set terminal png size 800,500 font ',14' linewidth 2\n");
    // formatting
    fprintf(gnuplot_pipe, "set output '../%s.png'\n", file_name.c_str());
    fprintf(gnuplot_pipe, "set xlabel 'Time [s]'\n");
    fprintf(gnuplot_pipe, "set ylabel 'Distance [m]'\n");

    fprintf(gnuplot_pipe,
            "set title 'Phased array ground station connectivity "
            "with %d beams'\n",
            input_ground_station.num_beams_);
    fprintf(gnuplot_pipe, "set key outside\n");
    int x_max_plot_window =
        std::floor(input_sim_parameters.total_sim_time * 1.05);
    fprintf(gnuplot_pipe, "set xrange[0:%d]\n", x_max_plot_window);
    // plotting

    // first satellite
    Satellite current_satellite = input_satellite_vector.at(0);
    if (input_satellite_vector.size() == 1) {
      if (current_satellite.plotting_color_.size() > 0) {
        fprintf(gnuplot_pipe,
                "plot '-' using 1:2 with lines lw 1 lc rgb '%s' title '%s' \n",
                current_satellite.plotting_color_.c_str(),
                current_satellite.get_name().c_str());
      } else {
        fprintf(gnuplot_pipe,
                "plot '-' using 1:2 with lines lw 1 title '%s' \n",
                current_satellite.get_name().c_str());
      }

    }

    else {
      if (current_satellite.plotting_color_.size() > 0) {
        fprintf(gnuplot_pipe,
                "plot '-' using 1:2 with lines lw 1 lc rgb '%s' title '%s'\\\n",
                current_satellite.plotting_color_.c_str(),
                current_satellite.get_name().c_str());
      } else {
        fprintf(gnuplot_pipe,
                "plot '-' using 1:2 with lines lw 1 title '%s'\\\n",
                current_satellite.get_name().c_str());
      }
    }

    for (size_t satellite_index = 1;
         satellite_index < input_satellite_vector.size(); satellite_index++) {
      current_satellite = input_satellite_vector.at(satellite_index);
      if (satellite_index < input_satellite_vector.size() - 1) {
        if (current_satellite.plotting_color_.size() > 0) {
          fprintf(gnuplot_pipe,
                  ",'-' using 1:2 with lines lw 1 lc rgb '%s' title '%s' \\\n",
                  current_satellite.plotting_color_.c_str(),
                  current_satellite.get_name().c_str());
        } else {
          fprintf(gnuplot_pipe,
                  ",'-' using 1:2 with lines lw 1 title '%s' \\\n",
                  current_satellite.get_name().c_str());
        }

      }

      else {
        if (current_satellite.plotting_color_.size() > 0) {
          fprintf(gnuplot_pipe,
                  ",'-' using 1:2 with lines lw 1 lc rgb '%s' title '%s'\n",
                  current_satellite.plotting_color_.c_str(),
                  current_satellite.get_name().c_str());
        } else {
          fprintf(gnuplot_pipe, ",'-' using 1:2 with lines lw 1 title '%s'\n",
                  current_satellite.get_name().c_str());
        }
      }
    }

    // now the orbit data, inline, one satellite at a time
    std::cout << "Running simulation...\n";
    for (size_t satellite_index = 0;
         satellite_index < input_satellite_vector.size(); satellite_index++) {
      Satellite current_satellite = input_satellite_vector.at(satellite_index);
      double current_satellite_time =
          current_satellite.get_instantaneous_time();
      double previous_time = current_satellite_time - 1;
      double ground_station_beam_angle =
          input_ground_station.angle_to_satellite_from_normal(
              current_satellite);
      double initial_distance;
      if (ground_station_beam_angle >
          input_ground_station.max_beam_angle_from_normal_) {
        fprintf(gnuplot_pipe, "%.17g NaN\n", current_satellite_time);
      } else {
        initial_distance =
            input_ground_station.distance_to_satellite(current_satellite);
        fprintf(gnuplot_pipe, "%.17g %.17g\n", current_satellite_time,
                initial_distance);
      }

      double evolved_distance = {0};

      double timestep_to_use = input_sim_parameters.initial_timestep_guess;
      current_satellite_time = current_satellite.get_instantaneous_time();
      while (current_satellite_time < input_sim_parameters.total_sim_time) {
        std::pair<double, double> drag_elements = {input_sim_parameters.F_10,
                                                   input_sim_parameters.A_p};
        std::pair<double, int> new_timestep_and_error_code =
            current_satellite.evolve_RK45(
                input_sim_parameters.epsilon, timestep_to_use,
                input_sim_parameters.perturbation_bool,
                input_sim_parameters.drag_bool, drag_elements);
        double new_timestep = new_timestep_and_error_code.first;
        int error_code = new_timestep_and_error_code.second;

        if (error_code != 0) {
          std::cout << "Error code " << error_code
                    << " detected, halting visualization\n";
          fprintf(gnuplot_pipe, "e\n");
          fprintf(gnuplot_pipe, "exit \n");
          pclose(gnuplot_pipe);
          return;
        }
        timestep_to_use = new_timestep;
        evolved_distance =
            input_ground_station.distance_to_satellite(current_satellite);
        previous_time = current_satellite_time;
        current_satellite_time = current_satellite.get_instantaneous_time();
        ground_station_beam_angle =
            input_ground_station.angle_to_satellite_from_normal(
                current_satellite);
        // Check number of existing connections to ground station at this point
        int num_sats_at_this_time =
            input_ground_station.num_sats_connected_at_this_time(
                current_satellite_time);
        if ((ground_station_beam_angle >
             input_ground_station.max_beam_angle_from_normal_) ||
            (num_sats_at_this_time == input_ground_station.num_beams_)) {
          fprintf(gnuplot_pipe, "%.17g NaN\n", current_satellite_time);
        } else {
          evolved_distance =
              input_ground_station.distance_to_satellite(current_satellite);
          fprintf(gnuplot_pipe, "%.17g %.17g\n", current_satellite_time,
                  evolved_distance);
          input_ground_station.update_linked_sats_map(
              satellite_index, current_satellite_time, previous_time);
        }
      }
      fprintf(gnuplot_pipe, "e\n");
    }
    fprintf(gnuplot_pipe, "pause mouse keypress\n");

    fprintf(gnuplot_pipe, "exit \n");
    pclose(gnuplot_pipe);
    std::cout << "Done\n";

  } else {
    std::cout << "gnuplot not found";
  }

  return;
}

void sim_and_plot_gs_connectivity_gnuplot(
    PhasedArrayGroundStation input_ground_station,
    std::vector<Satellite> input_satellite_vector,
    const SimParameters& input_sim_parameters, const std::string file_name) {
  // Objective: given an input ground station and satellite vector,
  // plot contact distances between ground station and satellites over time

  if (input_satellite_vector.size() < 1) {
    std::cout << "No input Satellite objects\n";
    return;
  }

  // first, open "pipe" to gnuplot
  FILE* gnuplot_pipe = popen("gnuplot", "w");
  // if it exists
  if (gnuplot_pipe) {
    fprintf(gnuplot_pipe,
            "set terminal png size 800,500 font ',14' linewidth 2\n");
    // formatting
    fprintf(gnuplot_pipe, "set output '../%s.png'\n", file_name.c_str());
    fprintf(gnuplot_pipe, "set xlabel 'Time [s]'\n");
    fprintf(gnuplot_pipe, "set ylabel 'Satellite Index'\n");

    fprintf(gnuplot_pipe,
            "set title 'Phased array ground station connectivity "
            "with %d beams'\n",
            input_ground_station.num_beams_);
    fprintf(gnuplot_pipe, "set key outside\n");
    int x_max_plot_window =
        std::floor(input_sim_parameters.total_sim_time * 1.05);
    fprintf(gnuplot_pipe, "set xrange[0:%d]\n", x_max_plot_window);
    fprintf(gnuplot_pipe, "set yrange[-0.5:%lu]\n",
            input_satellite_vector.size());

    // plotting

    // first satellite
    Satellite current_satellite = input_satellite_vector.at(0);
    if (input_satellite_vector.size() == 1) {
      if (current_satellite.plotting_color_.size() > 0) {
        fprintf(gnuplot_pipe,
                "plot '-' using 1:2 with lines lw 1 lc rgb '%s' title '%s' \n",
                current_satellite.plotting_color_.c_str(),
                current_satellite.get_name().c_str());
      } else {
        fprintf(gnuplot_pipe,
                "plot '-' using 1:2 with lines lw 1 title '%s' \n",
                current_satellite.get_name().c_str());
      }

    }

    else {
      if (current_satellite.plotting_color_.size() > 0) {
        fprintf(gnuplot_pipe,
                "plot '-' using 1:2 with lines lw 1 lc rgb '%s' title '%s'\\\n",
                current_satellite.plotting_color_.c_str(),
                current_satellite.get_name().c_str());
      } else {
        fprintf(gnuplot_pipe,
                "plot '-' using 1:2 with lines lw 1 title '%s'\\\n",
                current_satellite.get_name().c_str());
      }
    }

    for (size_t satellite_index = 1;
         satellite_index < input_satellite_vector.size(); satellite_index++) {
      current_satellite = input_satellite_vector.at(satellite_index);
      if (satellite_index < input_satellite_vector.size() - 1) {
        if (current_satellite.plotting_color_.size() > 0) {
          fprintf(gnuplot_pipe,
                  ",'-' using 1:2 with lines lw 1 lc rgb '%s' title '%s' \\\n",
                  current_satellite.plotting_color_.c_str(),
                  current_satellite.get_name().c_str());
        } else {
          fprintf(gnuplot_pipe,
                  ",'-' using 1:2 with lines lw 1 title '%s' \\\n",
                  current_satellite.get_name().c_str());
        }

      }

      else {
        if (current_satellite.plotting_color_.size() > 0) {
          fprintf(gnuplot_pipe,
                  ",'-' using 1:2 with lines lw 1 lc rgb '%s' title '%s'\n",
                  current_satellite.plotting_color_.c_str(),
                  current_satellite.get_name().c_str());
        } else {
          fprintf(gnuplot_pipe, ",'-' using 1:2 with lines lw 1 title '%s'\n",
                  current_satellite.get_name().c_str());
        }
      }
    }

    // now the orbit data, inline, one satellite at a time
    std::cout << "Running simulation...\n";
    for (size_t satellite_index = 0;
         satellite_index < input_satellite_vector.size(); satellite_index++) {
      Satellite current_satellite = input_satellite_vector.at(satellite_index);
      double current_satellite_time =
          current_satellite.get_instantaneous_time();
      double previous_time = current_satellite_time - 1;
      double ground_station_beam_angle =
          input_ground_station.angle_to_satellite_from_normal(
              current_satellite);

      if (ground_station_beam_angle >
          input_ground_station.max_beam_angle_from_normal_) {
        fprintf(gnuplot_pipe, "%.17g NaN\n", current_satellite_time);
      } else {
        fprintf(gnuplot_pipe, "%.17g %lu\n", current_satellite_time,
                satellite_index);
      }

      double timestep_to_use = input_sim_parameters.initial_timestep_guess;
      current_satellite_time = current_satellite.get_instantaneous_time();
      while (current_satellite_time < input_sim_parameters.total_sim_time) {
        std::pair<double, double> drag_elements = {input_sim_parameters.F_10,
                                                   input_sim_parameters.A_p};
        std::pair<double, int> new_timestep_and_error_code =
            current_satellite.evolve_RK45(
                input_sim_parameters.epsilon, timestep_to_use,
                input_sim_parameters.perturbation_bool,
                input_sim_parameters.drag_bool, drag_elements);
        double new_timestep = new_timestep_and_error_code.first;
        int error_code = new_timestep_and_error_code.second;

        if (error_code != 0) {
          std::cout << "Error code " << error_code
                    << " detected, halting visualization\n";
          fprintf(gnuplot_pipe, "e\n");
          fprintf(gnuplot_pipe, "exit \n");
          pclose(gnuplot_pipe);
          return;
        }
        timestep_to_use = new_timestep;
        previous_time = current_satellite_time;
        current_satellite_time = current_satellite.get_instantaneous_time();
        ground_station_beam_angle =
            input_ground_station.angle_to_satellite_from_normal(
                current_satellite);
        // Check number of existing connections to ground station at this point
        int num_sats_at_this_time =
            input_ground_station.num_sats_connected_at_this_time(
                current_satellite_time);
        if ((ground_station_beam_angle >
             input_ground_station.max_beam_angle_from_normal_) ||
            (num_sats_at_this_time == input_ground_station.num_beams_)) {
          fprintf(gnuplot_pipe, "%.17g NaN\n", current_satellite_time);
        } else {
          fprintf(gnuplot_pipe, "%.17g %lu\n", current_satellite_time,
                  satellite_index);
          input_ground_station.update_linked_sats_map(
              satellite_index, current_satellite_time, previous_time);
        }
      }
      fprintf(gnuplot_pipe, "e\n");
    }
    // for (auto key_val_pair : input_ground_station.linked_sats_map_) {
    //   std::cout << "==========================================\n";
    //   std::cout << "Satellite index " << key_val_pair.first << "\n";
    //   for (std::pair<double,double> range : key_val_pair.second) {
    //     std::cout << "{" << range.first << "," << range.second << "}\n";
    //   }
    // }
    fprintf(gnuplot_pipe, "pause mouse keypress\n");

    fprintf(gnuplot_pipe, "exit \n");
    pclose(gnuplot_pipe);
    std::cout << "Done\n";

  } else {
    std::cout << "gnuplot not found";
  }

  return;
}

int add_lowthrust_orbit_transfer(Satellite& input_satellite_object,
                                 const double final_orbit_semimajor_axis_km,
                                 const double input_thrust_magnitude,
                                 const double transfer_initiation_time) {
  // Only transfers between circular orbits are supported at this time
  // Ref: https://prussing.ae.illinois.edu/AE402/low.thrust.pdf
  int error_code = 0;
  double initial_eccentricity =
      input_satellite_object.get_orbital_parameter("Eccentricity");
  if (initial_eccentricity != 0) {
    std::cout << "Satellite's initial orbit was not circular\n";
    error_code = 1;
  }
  double satellite_mass = input_satellite_object.get_mass();
  if (satellite_mass <= 0) {
    std::cout << "Error: satellite mass was <= 0\n";
    error_code = 2;
  }
  double thrust_acceleration = input_thrust_magnitude / satellite_mass;
  double semimajor_axis_final = 1000 * final_orbit_semimajor_axis_km;  // m
  double semimajor_axis_initial =
      input_satellite_object.get_orbital_parameter("Semimajor Axis");
  const double mu_Earth = G * mass_Earth;
  double comp1 = sqrt(mu_Earth / semimajor_axis_initial);
  double comp2 = sqrt(mu_Earth / semimajor_axis_final);
  double time_to_burn = (comp1 - comp2) / thrust_acceleration;

  // Thrust is purely co-linear with velocity vector, so in the +- x direction
  // of the LVLH frame (along +x if raising orbit, -x if lowering orbit)
  std::array<double, 3> LVLH_thrust_direction;
  if (semimajor_axis_initial < semimajor_axis_final) {
    LVLH_thrust_direction = {1, 0, 0};
  } else if (semimajor_axis_initial > semimajor_axis_final) {
    LVLH_thrust_direction = {-1, 0, 0};
    time_to_burn =
        (-1) * time_to_burn;  // Since this would have otherwise been negative
  } else {
    std::cout << "Error: initial and final semimajor axes were equal.\n";
    error_code = 3;  // Arbitrarily choose a burn direction in this scenario,
                     // since it shouldn't matter
  }
  if (error_code == 0) {
    input_satellite_object.add_LVLH_thrust_profile(
        LVLH_thrust_direction, input_thrust_magnitude, transfer_initiation_time,
        transfer_initiation_time + time_to_burn);
  }
  return error_code;
}