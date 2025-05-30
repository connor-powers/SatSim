
#include <iostream>
#define _USE_MATH_DEFINES
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cmath>

#include "Satellite.h"
#include "utils.h"
using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::Vector3d;

// Currently baselining use of cartesian coordinates in ECI frame (J2000
// specifically, if that comes up later)

// Objective: given the semimajor axis of an orbit, calculate its orbital period
double Satellite::calculate_orbital_period() {
  // https://en.wikipedia.org/wiki/Orbital_period
  double T = 2 * M_PI * sqrt(pow(a_, 3) / (G * mass_Earth));
  return T;
}

// Objective: calculate the perifocal position of the satellite
std::array<double, 3> Satellite::calculate_perifocal_position() {
  // Using approach from Fundamentals of Astrodynamics
  std::array<double, 3> calculated_perifocal_position;
  double p =
      a_ * (1 - eccentricity_) * (1 + eccentricity_);  // rearranging eq. 1-47

  double r = p / (1 + eccentricity_ * cos(true_anomaly_));

  double r_perifocal_P = r * cos(true_anomaly_);
  double r_perifocal_Q = r * sin(true_anomaly_);

  calculated_perifocal_position.at(0) = r_perifocal_P;
  calculated_perifocal_position.at(1) = r_perifocal_Q;
  calculated_perifocal_position.at(2) = 0;  // W component is 0

  return calculated_perifocal_position;
}

// Objective: calculate the perifocal velocity of the satellite
std::array<double, 3> Satellite::calculate_perifocal_velocity() {
  // Using approach from Fundamentals of Astrodynamics
  std::array<double, 3> calculated_perifocal_velocity;
  double p =
      a_ * (1 - eccentricity_) * (1 + eccentricity_);  // rearranging eq. 1-47
  const double mu_Earth = G * mass_Earth;
  double v_perifocal_P = sqrt(mu_Earth / p) * (-sin(true_anomaly_));
  double v_perifocal_Q =
      sqrt(mu_Earth / p) * (eccentricity_ + cos(true_anomaly_));

  calculated_perifocal_velocity.at(0) = v_perifocal_P;
  calculated_perifocal_velocity.at(1) = v_perifocal_Q;
  calculated_perifocal_velocity.at(2) = 0;  // 0 component in the W direction

  return calculated_perifocal_velocity;
}

// Objective: convert a vector from perifocal frame to ECI coordinates
std::array<double, 3> Satellite::convert_perifocal_to_ECI(
    const std::array<double, 3> input_perifocal_vec) {
  // Method from Fundamentals of Astrodynamics
  // Sounds like what they describe as their "Geocentric-Equatorial" coordinate
  // system is ECI
  Matrix3d R_tilde_matrix;

  R_tilde_matrix(0, 0) =
      cos(raan_) * cos(arg_of_periapsis_) -
      sin(raan_) * sin(arg_of_periapsis_) * cos(inclination_);

  R_tilde_matrix(0, 1) =
      -cos(raan_) * sin(arg_of_periapsis_) -
      sin(raan_) * cos(arg_of_periapsis_) * cos(inclination_);

  R_tilde_matrix(0, 2) = sin(raan_) * sin(inclination_);

  R_tilde_matrix(1, 0) =
      sin(raan_) * cos(arg_of_periapsis_) +
      cos(raan_) * sin(arg_of_periapsis_) * cos(inclination_);

  R_tilde_matrix(1, 1) =
      -sin(raan_) * sin(arg_of_periapsis_) +
      cos(raan_) * cos(arg_of_periapsis_) * cos(inclination_);

  R_tilde_matrix(1, 2) = -cos(raan_) * sin(inclination_);

  R_tilde_matrix(2, 0) = sin(arg_of_periapsis_) * sin(inclination_);

  R_tilde_matrix(2, 1) = cos(arg_of_periapsis_) * sin(inclination_);

  R_tilde_matrix(2, 2) = cos(inclination_);

  // convert to Eigen vectors just to make sure matrix-vector multiplication is
  // done correctly
  Vector3d vector_ijk(0, 0, 0);
  Vector3d vector_pqw(0, 0, 0);

  vector_pqw(0) = input_perifocal_vec.at(0);
  vector_pqw(1) = input_perifocal_vec.at(1);
  vector_pqw(2) = input_perifocal_vec.at(2);

  vector_ijk = R_tilde_matrix * vector_pqw;
  std::array<double, 3> output_vector_ijk = {vector_ijk(0), vector_ijk(1),
                                             vector_ijk(2)};

  return output_vector_ijk;
}

// Objective: convert vector from ECI frame to perifocal frame
std::array<double, 3> Satellite::convert_ECI_to_perifocal(
    const std::array<double, 3> input_ECI_vec) {
  // Method from Fundamentals of Astrodynamics, and using
  // https://en.wikipedia.org/wiki/Perifocal_coordinate_system Sounds like what
  // they describe as their "Geocentric-Equatorial" coordinate system is ECI
  Matrix3d R_tilde_matrix;

  R_tilde_matrix(0, 0) =
      cos(raan_) * cos(arg_of_periapsis_) -
      sin(raan_) * sin(arg_of_periapsis_) * cos(inclination_);

  R_tilde_matrix(0, 1) =
      -cos(raan_) * sin(arg_of_periapsis_) -
      sin(raan_) * cos(arg_of_periapsis_) * cos(inclination_);

  R_tilde_matrix(0, 2) = sin(raan_) * sin(inclination_);

  R_tilde_matrix(1, 0) =
      sin(raan_) * cos(arg_of_periapsis_) +
      cos(raan_) * sin(arg_of_periapsis_) * cos(inclination_);

  R_tilde_matrix(1, 1) =
      -sin(raan_) * sin(arg_of_periapsis_) +
      cos(raan_) * cos(arg_of_periapsis_) * cos(inclination_);

  R_tilde_matrix(1, 2) = -cos(raan_) * sin(inclination_);

  R_tilde_matrix(2, 0) = sin(arg_of_periapsis_) * sin(inclination_);

  R_tilde_matrix(2, 1) = cos(arg_of_periapsis_) * sin(inclination_);

  R_tilde_matrix(2, 2) = cos(inclination_);

  // convert to Eigen vectors just to make sure matrix-vector multiplication is
  // done correctly
  Vector3d vector_ijk(0, 0, 0);
  Vector3d vector_pqw(0, 0, 0);
  vector_ijk(0) = input_ECI_vec.at(0);
  vector_ijk(1) = input_ECI_vec.at(1);
  vector_ijk(2) = input_ECI_vec.at(2);

  vector_pqw = (R_tilde_matrix.inverse()) * vector_ijk;
  std::array<double, 3> output_vector_pqw = {vector_pqw(0), vector_pqw(1),
                                             vector_pqw(2)};

  return output_vector_pqw;
}

// Objective: evolve position and velocity (combined into one vector) one
// timestep via RK4 method
// void Satellite::evolve_RK4(const double input_step_size) {
//   // format input position and velocity arrays into single array for RK4 step
//   std::array<double, 6> combined_initial_position_and_velocity_array = {};
//   std::pair<std::array<double, 3>, std::array<double, 3>>
//       output_position_velocity_pair = {};

//   for (size_t ind = 0; ind < 3; ind++) {
//     combined_initial_position_and_velocity_array.at(ind) =
//         ECI_position_.at(ind);
//   }
//   for (size_t ind = 3; ind < 6; ind++) {
//     combined_initial_position_and_velocity_array.at(ind) =
//         ECI_velocity_.at(ind - 3);
//   }

//   // populate list of thrust forces at half a timestep past, for RK4
//   // calculations
//   std::vector<std::array<double, 3>>
//   list_of_LVLH_forces_at_half_timestep_past =
//       {};
//   std::vector<std::array<double, 3>> list_of_ECI_forces_at_half_timestep_past
//   =
//       {};

//   for (const ThrustProfileLVLH thrust_profile : thrust_profile_list_) {
//     if (((t_ + (input_step_size / 2.0)) >= thrust_profile.t_start_) &&
//         ((t_ + (input_step_size / 2.0)) <= thrust_profile.t_end_)) {
//       list_of_LVLH_forces_at_half_timestep_past.push_back(
//           thrust_profile.LVLH_force_vec_);
//       std::array<double, 3> ECI_thrust_vector = convert_LVLH_to_ECI_manual(
//           thrust_profile.LVLH_force_vec_, ECI_position_, ECI_velocity_);
//       list_of_ECI_forces_at_half_timestep_past.push_back(ECI_thrust_vector);
//     }
//   }

//   // populate list of thrust forces at a timestep past, for RK4 calculations
//   std::vector<std::array<double, 3>> list_of_LVLH_forces_at_one_timestep_past
//   =
//       {};
//   std::vector<std::array<double, 3>> list_of_ECI_forces_at_one_timestep_past
//   =
//       {};

//   for (const ThrustProfileLVLH thrust_profile : thrust_profile_list_) {
//     if (((t_ + input_step_size) >= thrust_profile.t_start_) &&
//         ((t_ + input_step_size) <= thrust_profile.t_end_)) {
//       list_of_LVLH_forces_at_one_timestep_past.push_back(
//           thrust_profile.LVLH_force_vec_);
//       std::array<double, 3> ECI_thrust_vector = convert_LVLH_to_ECI_manual(
//           thrust_profile.LVLH_force_vec_, ECI_position_, ECI_velocity_);
//       list_of_ECI_forces_at_one_timestep_past.push_back(ECI_thrust_vector);
//     }
//   }

//   std::array<double, 6> output_combined_position_and_velocity_array =
//       RK4_step<6>(combined_initial_position_and_velocity_array,
//       input_step_size,
//                   RK4_deriv_function_orbit_position_and_velocity, m_,
//                   list_of_ECI_forces_at_this_time_,
//                   list_of_ECI_forces_at_half_timestep_past,
//                   list_of_ECI_forces_at_one_timestep_past);
//   // std::array<double,6> output_combined_angular_array=
//   //
//   RK4_step<6>(combined_initial_angular_array,input_step_size,RK4_deriv_function_angular,I_,list_of_body_frame_torques_at_this_time_,list_of_body_frame_torques_at_half_timestep_past,list_of_body_frame_torques_at_one_timestep_past);

//   for (size_t ind = 0; ind < 3; ind++) {
//     ECI_position_.at(ind) =
//     output_combined_position_and_velocity_array.at(ind);
//     ECI_velocity_.at(ind) =
//         output_combined_position_and_velocity_array.at(ind + 3);
//     // also update the perifocal versions
//     perifocal_position_ = convert_ECI_to_perifocal(ECI_position_);
//     perifocal_velocity_ = convert_ECI_to_perifocal(ECI_velocity_);
//   }
//   t_ += input_step_size;

//   list_of_LVLH_forces_at_this_time_ =
//   list_of_LVLH_forces_at_one_timestep_past; list_of_ECI_forces_at_this_time_
//   = list_of_ECI_forces_at_one_timestep_past;

//   // Update orbital parameters
//   update_orbital_elements_from_position_and_velocity();

//   return;
// }

// Objective: add a LVLH frame thrust profile to the satellite
void Satellite::add_LVLH_thrust_profile(
    const std::array<double, 3> input_LVLH_thrust_vector,
    const double input_thrust_start_time, const double input_thrust_end_time) {
  ThrustProfileLVLH new_thrust_profile(
      input_thrust_start_time, input_thrust_end_time, input_LVLH_thrust_vector);
  thrust_profile_list_.push_back(new_thrust_profile);
  if (input_thrust_start_time == t_) {
    list_of_LVLH_forces_at_this_time_.push_back(input_LVLH_thrust_vector);
    std::array<double, 3> ECI_thrust_vector = convert_LVLH_to_ECI_manual(
        input_LVLH_thrust_vector, ECI_position_, ECI_velocity_);
    list_of_ECI_forces_at_this_time_.push_back(ECI_thrust_vector);
  }
}

// Alternate input argument style
void Satellite::add_LVLH_thrust_profile(
    const std::array<double, 3> input_LVLH_normalized_thrust_direction,
    const double input_LVLH_thrust_magnitude,
    const double input_thrust_start_time, const double input_thrust_end_time) {
  ThrustProfileLVLH new_thrust_profile(
      input_thrust_start_time, input_thrust_end_time,
      input_LVLH_normalized_thrust_direction, input_LVLH_thrust_magnitude);
  thrust_profile_list_.push_back(new_thrust_profile);

  std::array<double, 3> LVLH_thrust_vec = {0, 0, 0};

  for (size_t ind = 0; ind < 3; ind++) {
    LVLH_thrust_vec.at(ind) = input_LVLH_thrust_magnitude *
                              input_LVLH_normalized_thrust_direction.at(ind);
  }

  if (input_thrust_start_time == t_) {
    list_of_LVLH_forces_at_this_time_.push_back(LVLH_thrust_vec);
    std::array<double, 3> ECI_thrust_vector = convert_LVLH_to_ECI_manual(
        LVLH_thrust_vec, ECI_position_, ECI_velocity_);
    list_of_ECI_forces_at_this_time_.push_back(ECI_thrust_vector);
  }
}

// Now the version for argument of periapsis change maneuvers
void Satellite::add_LVLH_thrust_profile(const double input_thrust_start_time,
                                        const double final_arg_of_periapsis,
                                        const double input_thrust_magnitude) {
  ThrustProfileLVLH thrust_profile(
      input_thrust_start_time, final_arg_of_periapsis, input_thrust_magnitude,
      arg_of_periapsis_, eccentricity_, a_, m_);
  thrust_profile_list_.push_back(thrust_profile);
}

int Satellite::update_orbital_elements_from_position_and_velocity() {
  // Anytime the orbit is changed via external forces, need to update the
  // orbital parameters of the satellite. True anomaly should change over time
  // even in absence of external forces Using approach from Fundamentals of
  // Astrodynamics
  int error_code = 0;  // 0 represents nominal operation
  double mu_Earth = G * mass_Earth;
  Vector3d position_vector = {ECI_position_.at(0), ECI_position_.at(1),
                              ECI_position_.at(2)};
  Vector3d velocity_vector = {ECI_velocity_.at(0), ECI_velocity_.at(1),
                              ECI_velocity_.at(2)};
  Vector3d h_vector = position_vector.cross(velocity_vector);
  double h = h_vector.norm();

  Vector3d n_vector = {-h_vector(1), h_vector(0), 0};
  double n = n_vector.norm();

  double v_magnitude = get_speed();
  double r_magnitude = get_radius();

  Vector3d e_vec_component_1 =
      (1.0 / mu_Earth) *
      (v_magnitude * v_magnitude - (mu_Earth / r_magnitude)) * position_vector;
  Vector3d e_vec_component_2 = (1.0 / mu_Earth) *
                               (position_vector.dot(velocity_vector)) *
                               velocity_vector;
  Vector3d e_vec = e_vec_component_1 - e_vec_component_2;
  double calculated_eccentricity = e_vec.norm();

  double calculated_p = h * h / mu_Earth;

  double calculated_i = acos(h_vector(2) / h);

  double calculated_RAAN = acos(n_vector(0) / n);
  if (n_vector(1) < 0) {
    calculated_RAAN = 2 * M_PI - calculated_RAAN;
  }

  // Need to treat the e \approx 0 case (circular orbits) specially
  double calculated_arg_of_periapsis;
  double calculated_true_anomaly;
  if (calculated_eccentricity > pow(10, -15)) {
    calculated_arg_of_periapsis =
        acos(n_vector.dot(e_vec) / (n * calculated_eccentricity));
    if (e_vec(2) < 0) {
      calculated_arg_of_periapsis = 2 * M_PI - calculated_arg_of_periapsis;
    }
    calculated_true_anomaly = acos(e_vec.dot(position_vector) /
                                   (calculated_eccentricity * r_magnitude));
    if (std::isnan(calculated_true_anomaly)) {
      // First try arg of the following acos call to a float to avoid bug
      // where it ended up being represented by a value with magnitude
      // slightly larger than one, causing NaN in output of acos
      float true_anomaly_acos_arg_float =
          e_vec.dot(position_vector) / (calculated_eccentricity * r_magnitude);
      calculated_true_anomaly = acos(true_anomaly_acos_arg_float);
      if (std::isnan(calculated_true_anomaly)) {
        std::cout << "NaN true anomaly encountered.\n";
        error_code = 10;
      }
    }
    if (position_vector.dot(velocity_vector) < 0) {
      calculated_true_anomaly = 2 * M_PI - calculated_true_anomaly;
    }

  } else {
    // Approximately circular orbits
    // For this case, I'll set the true anomaly to be the argument of latitude
    // (which is valid when arg of periapsis is 0, which is what I'm setting it
    // to) Refs: https://en.wikipedia.org/wiki/True_anomaly#Circular_orbit ,
    // https://en.wikipedia.org/wiki/Argument_of_latitude ,
    // https://en.wikipedia.org/wiki/Argument_of_periapsis
    calculated_arg_of_periapsis = 0;  // Setting this

    double calculated_arg_of_latitude =
        acos(n_vector.dot(position_vector) / (n * r_magnitude));
    if (position_vector(2) < 0) {
      calculated_arg_of_latitude = (2 * M_PI - calculated_arg_of_latitude);
    }
    calculated_true_anomaly = calculated_arg_of_latitude;  // For this case
  }

  double calculated_a =
      calculated_p / (1.0 - calculated_eccentricity * calculated_eccentricity);

  // Update stored values of these orbital elements

  a_ = calculated_a;
  eccentricity_ = calculated_eccentricity;
  inclination_ = calculated_i;
  raan_ = calculated_RAAN;
  arg_of_periapsis_ = calculated_arg_of_periapsis;
  true_anomaly_ = calculated_true_anomaly;

  return error_code;
}

std::array<double, 6> Satellite::get_orbital_elements() {
  std::array<double, 6> orbit_elems_array;
  orbit_elems_array.at(0) = a_;
  orbit_elems_array.at(1) = eccentricity_;
  orbit_elems_array.at(2) = inclination_;
  orbit_elems_array.at(3) = raan_;
  orbit_elems_array.at(4) = arg_of_periapsis_;
  orbit_elems_array.at(5) = true_anomaly_;
  return orbit_elems_array;
}

std::pair<double, int> Satellite::evolve_RK45(
    const double input_epsilon, const double input_step_size,
    const bool perturbation, const bool atmospheric_drag,
    const std::pair<double, double> drag_elements) {
  // std::cout << "In EvolveRK45, received perturbation bool: " << perturbation
  // << "\n";
  // perturbation is a flag which, when set to true, currently accounts for J2
  // perturbation.

  // Let's do a single RK45_step call with y_n combined between orbital motion
  // and attitude variables,
  //  y_n = {ECI_position_x, ECI_position_y, ECI_position_z, ECI_velocity_x,
  //  ECI_velocity_y, ECI_velocity_z, q_0, q_1, q_2, q_3, omega_x, omega_y,
  //  omega_z} Where the quaternion represents the attitude of the spacecraft
  //  body frame with respect to the LVLH frame, and omega_i is the angular
  //  velocity around the ith axis of the body frame with respect to the LVLH
  //  frame, represented in the body frame

  // The tuple drag_elements contains the F_10 value and the A_p value used for
  // atmospheric drag calculations, if applicable
  // F_10 is the first element, A_p is the second element
  double input_F_10 = drag_elements.first;
  double input_A_p = drag_elements.second;

  // Initialize arg of periapsis change thrust profile from current orbital
  // parameters, if there is one that hasn't already been initialized
  check_for_maneuvers_to_initialize();
  // Add any thrusts from an argument of periapsis change maneuver, if
  // applicable
  int arg_of_periapsis_code = add_arg_of_periapsis_change_thrust();

  std::array<double, 13>
      combined_initial_position_velocity_quaternion_angular_velocity_array = {};
  std::pair<std::array<double, 3>, std::array<double, 3>>
      output_position_velocity_pair = {};

  for (size_t ind = 0; ind < 3; ind++) {
    combined_initial_position_velocity_quaternion_angular_velocity_array.at(
        ind) = ECI_position_.at(ind);
  }
  for (size_t ind = 3; ind < 6; ind++) {
    combined_initial_position_velocity_quaternion_angular_velocity_array.at(
        ind) = ECI_velocity_.at(ind - 3);
  }
  for (size_t ind = 6; ind < 10; ind++) {
    combined_initial_position_velocity_quaternion_angular_velocity_array.at(
        ind) = quaternion_satellite_bodyframe_wrt_LVLH_.at(ind - 6);
  }

  for (size_t ind = 10; ind < 13; ind++) {
    combined_initial_position_velocity_quaternion_angular_velocity_array.at(
        ind) = body_angular_velocity_vec_wrt_LVLH_in_body_frame_.at(ind - 10);
  }

  Matrix3d LVLH_to_body_transformation_matrix =
      LVLH_to_body_transformation_matrix_from_quaternion(
          quaternion_satellite_bodyframe_wrt_LVLH_);

  Vector3d body_angular_velocity_vec_wrt_LVLH_in_body_frame_eigenform = {
      body_angular_velocity_vec_wrt_LVLH_in_body_frame_.at(0),
      body_angular_velocity_vec_wrt_LVLH_in_body_frame_.at(1),
      body_angular_velocity_vec_wrt_LVLH_in_body_frame_.at(2)};

  Vector3d omega_I = calculate_omega_I(
      body_angular_velocity_vec_wrt_LVLH_in_body_frame_eigenform,
      LVLH_to_body_transformation_matrix, orbital_rate_);

  Vector3d omega_LVLH_wrt_inertial_in_LVLH = {0, -orbital_rate_, 0};

  Matrix3d J_matrix = construct_J_matrix(J_11_, J_22_, J_33_);

  std::pair<std::array<double, 13>, std::pair<double, double>> output_pair =
      RK45_step<13>(
          combined_initial_position_velocity_quaternion_angular_velocity_array,
          input_step_size,
          RK45_combined_orbit_position_velocity_attitude_deriv_function,
          J_matrix, bodyframe_torque_profile_list_, omega_I,
          orbital_angular_acceleration_, LVLH_to_body_transformation_matrix,
          omega_LVLH_wrt_inertial_in_LVLH, m_, thrust_profile_list_,
          inclination_, arg_of_periapsis_, true_anomaly_, input_F_10, input_A_p,
          A_s_, perturbation, atmospheric_drag, t_, input_epsilon);

  std::array<double, 13>
      output_combined_position_velocity_quaternion_angular_velocity_array =
          output_pair.first;
  double step_size_successfully_used_here = output_pair.second.first;
  double new_step_size = output_pair.second.second;

  for (size_t ind = 0; ind < 3; ind++) {
    ECI_position_.at(ind) =
        output_combined_position_velocity_quaternion_angular_velocity_array.at(
            ind);
    ECI_velocity_.at(ind) =
        output_combined_position_velocity_quaternion_angular_velocity_array.at(
            ind + 3);
    // Also update the perifocal versions
    perifocal_position_ = convert_ECI_to_perifocal(ECI_position_);
    perifocal_velocity_ = convert_ECI_to_perifocal(ECI_velocity_);
  }
  t_ += step_size_successfully_used_here;

  for (size_t ind = 0; ind < quaternion_satellite_bodyframe_wrt_LVLH_.size();
       ind++) {
    quaternion_satellite_bodyframe_wrt_LVLH_.at(ind) =
        output_combined_position_velocity_quaternion_angular_velocity_array.at(
            ind + 6);
  }
  quaternion_satellite_bodyframe_wrt_LVLH_ =
      normalize_quaternion(quaternion_satellite_bodyframe_wrt_LVLH_);
  std::array<double, 3> updated_roll_yaw_pitch =
      convert_quaternion_to_roll_yaw_pitch_angles(
          quaternion_satellite_bodyframe_wrt_LVLH_);

  roll_angle_ = updated_roll_yaw_pitch.at(0);
  yaw_angle_ = updated_roll_yaw_pitch.at(1);
  pitch_angle_ = updated_roll_yaw_pitch.at(2);

  for (size_t ind = 0;
       ind < body_angular_velocity_vec_wrt_LVLH_in_body_frame_.size(); ind++) {
    body_angular_velocity_vec_wrt_LVLH_in_body_frame_.at(ind) =
        output_combined_position_velocity_quaternion_angular_velocity_array.at(
            ind + 10);
  }

  // Update orbital parameters
  int orbit_elems_error_code =
      update_orbital_elements_from_position_and_velocity();
  orbital_rate_ = calculate_instantaneous_orbit_rate();
  orbital_angular_acceleration_ =
      calculate_instantaneous_orbit_angular_acceleration();
  std::pair<double, int> evolve_RK45_output_pair;
  // Orbital radius shouldn't be less than or equal to Earth's radius
  double new_orbital_radius = get_radius();
  if (new_orbital_radius <= radius_Earth) {
    orbit_elems_error_code = 2;
  }
  evolve_RK45_output_pair.first = new_step_size;
  evolve_RK45_output_pair.second = orbit_elems_error_code;

  // If you just added a temporally localized thrust profile stemming from an
  // argument of periapsis change thrust profile, remove it now that it's no
  // longer needed to avoid bloat of the thrust profile list
  if (arg_of_periapsis_code == 1) {
    thrust_profile_list_.pop_back();
  }
  return evolve_RK45_output_pair;
}

// Returns a specific orbital element
double Satellite::get_orbital_parameter(
    const std::string orbital_parameter_name) {
  if (orbital_parameter_name == "Semimajor Axis") {
    return a_;
  } else if (orbital_parameter_name == "Eccentricity") {
    return eccentricity_;
  } else if (orbital_parameter_name == "Inclination") {
    return inclination_ * (180 / M_PI);  // Returns val in degrees
  } else if (orbital_parameter_name == "RAAN") {
    return raan_ * (180.0 / M_PI);  // Returns val in degrees
  } else if (orbital_parameter_name == "Argument of Periapsis") {
    return arg_of_periapsis_ * (180.0 / M_PI);  // Returns val in degrees
  } else if (orbital_parameter_name == "True Anomaly") {
    return true_anomaly_ * (180.0 / M_PI);  // Returns val in degrees
  } else if (orbital_parameter_name == "Orbital Rate") {
    return orbital_rate_;
  } else if (orbital_parameter_name == "Orbital Angular Acceleration") {
    return orbital_angular_acceleration_;
  } else if (orbital_parameter_name == "Total Energy") {
    return get_total_energy();
  } else {
    throw std::invalid_argument("Value name not recognized");
  }
}

// Objective: add a body frame torque profile to the satellite
void Satellite::add_bodyframe_torque_profile(
    const std::array<double, 3> input_bodyframe_torque_vector,
    const double input_torque_start_time, const double input_torque_end_time) {
  BodyframeTorqueProfile new_torque_profile(input_torque_start_time,
                                            input_torque_end_time,
                                            input_bodyframe_torque_vector);
  bodyframe_torque_profile_list_.push_back(new_torque_profile);
  if (input_torque_start_time == 0) {
    list_of_body_frame_torques_at_this_time_.push_back(
        input_bodyframe_torque_vector);
  }
}
// Alternate input argument style
void Satellite::add_bodyframe_torque_profile(
    const std::array<double, 3> input_bodyframe_direction_unit_vec,
    const double input_bodyframe_torque_magnitude,
    const double input_torque_start_time, const double input_torque_end_time) {
  std::array<double, 3> input_bodyframe_torque_vector = {0, 0, 0};

  for (size_t ind = 0; ind < 3; ind++) {
    input_bodyframe_torque_vector.at(ind) =
        input_bodyframe_torque_magnitude *
        input_bodyframe_direction_unit_vec.at(ind);
  }
  BodyframeTorqueProfile new_torque_profile(input_torque_start_time,
                                            input_torque_end_time,
                                            input_bodyframe_torque_vector);
  bodyframe_torque_profile_list_.push_back(new_torque_profile);
  if (input_torque_start_time == 0) {
    list_of_body_frame_torques_at_this_time_.push_back(
        input_bodyframe_torque_vector);
  }
}

double Satellite::calculate_instantaneous_orbit_rate() {
  // Need to figure out how to calculate the equivalent of w_0 in section 4.3.1
  // of
  // https://ntrs.nasa.gov/api/citations/20240009554/downloads/Space%20Attitude%20Development%20Control.pdf
  // but for general elliptical orbits Is this proportional to the rate of
  // change of the true anomaly? Proof sketched out here:
  // https://space.stackexchange.com/questions/21615/how-to-calculate-the-complete-true-anomaly-motion-in-an-elliptical-orbit
  Vector3d position_vector = {perifocal_position_.at(0),
                              perifocal_position_.at(1),
                              perifocal_position_.at(2)};
  Vector3d velocity_vector = {perifocal_velocity_.at(0),
                              perifocal_velocity_.at(1),
                              perifocal_velocity_.at(2)};
  Vector3d h_vector = position_vector.cross(velocity_vector);
  double h = h_vector.norm();
  double r = get_radius();
  double orbit_rate = h / (r * r);
  return orbit_rate;
}

double Satellite::calculate_instantaneous_orbit_angular_acceleration() {
  // Need to figure out how to calculate the equivalent of w_0^dot to get the
  // w_lvlh^dot term in Ch.4 of
  // https://ntrs.nasa.gov/api/citations/20240009554/downloads/Space%20Attitude%20Development%20Control.pdf
  // but for general elliptical orbits Proof sketched out here:
  // https://space.stackexchange.com/questions/21615/how-to-calculate-the-complete-true-anomaly-motion-in-an-elliptical-orbit
  // This assumes angular momentum is constant
  Vector3d position_vector = {perifocal_position_.at(0),
                              perifocal_position_.at(1),
                              perifocal_position_.at(2)};
  Vector3d velocity_vector = {perifocal_velocity_.at(0),
                              perifocal_velocity_.at(1),
                              perifocal_velocity_.at(2)};
  Vector3d h_vector = position_vector.cross(velocity_vector);
  double h = h_vector.norm();
  double r = get_radius();
  Vector3d unit_position_vector = position_vector;
  unit_position_vector.normalize();
  double orbit_angular_acceleration =
      -2 * h / (r * r * r) * velocity_vector.dot(unit_position_vector);

  // Neglecting term relating to derivative of h for now (would be an h^dot /
  // r^2 term added to the above term)
  return orbit_angular_acceleration;
}

void Satellite::initialize_and_normalize_body_quaternion(
    const double roll_angle, const double pitch_angle, const double yaw_angle) {
  std::array<double, 4> quaternion =
      rollyawpitch_angles_to_quaternion(roll_angle, pitch_angle, yaw_angle);
  std::array<double, 4> normalized_quaternion =
      normalize_quaternion(quaternion);
  quaternion_satellite_bodyframe_wrt_LVLH_ = normalized_quaternion;
  return;
}

// Return a specific attitude-related value
double Satellite::get_attitude_val(std::string input_attitude_val_name) {
  if (input_attitude_val_name == "Roll") {
    return roll_angle_ * (180.0 / M_PI);  // Returns val in degrees
  } else if (input_attitude_val_name == "Pitch") {
    return pitch_angle_ * (180.0 / M_PI);  // Returns val in degrees
  } else if (input_attitude_val_name == "Yaw") {
    return yaw_angle_ * (180.0 / M_PI);  // Returns val in degrees
  } else if (input_attitude_val_name == "omega_x") {
    return body_angular_velocity_vec_wrt_LVLH_in_body_frame_.at(0) *
           (180.0 / M_PI);  // Returns val in degrees/s
  } else if (input_attitude_val_name == "omega_y") {
    return body_angular_velocity_vec_wrt_LVLH_in_body_frame_.at(1) *
           (180.0 / M_PI);  // Returns val in degrees/s
  } else if (input_attitude_val_name == "omega_z") {
    return body_angular_velocity_vec_wrt_LVLH_in_body_frame_.at(2) *
           (180.0 / M_PI);  // Returns val in degrees/s
  } else if (input_attitude_val_name == "q_0") {
    return quaternion_satellite_bodyframe_wrt_LVLH_.at(0);
  } else if (input_attitude_val_name == "q_1") {
    return quaternion_satellite_bodyframe_wrt_LVLH_.at(1);
  } else if (input_attitude_val_name == "q_2") {
    return quaternion_satellite_bodyframe_wrt_LVLH_.at(2);
  } else if (input_attitude_val_name == "q_3") {
    return quaternion_satellite_bodyframe_wrt_LVLH_.at(3);
  } else {
    throw std::invalid_argument("Value name not recognized");
  }
}

void Satellite::initialize_body_angular_velocity_vec_wrt_LVLH_in_body_frame() {
  std::array<double, 3> initial_body_angular_velocity_in_LVLH_frame = {
      0, orbital_rate_, 0};
  body_angular_velocity_vec_wrt_LVLH_in_body_frame_ =
      convert_array_from_LVLH_to_bodyframe(
          initial_body_angular_velocity_in_LVLH_frame, roll_angle_, yaw_angle_,
          pitch_angle_);  // Because LVLH is always rotating around y axis to
                          // keep z pointed towards Earth
}
int Satellite::add_arg_of_periapsis_change_thrust() {
  // Returning 0 means didn't find any arg of periapsis change thrust profiles
  // Returning 1 means it did find one and added a temporally localized thrust
  // at this time
  for (ThrustProfileLVLH thrust_profile : thrust_profile_list_) {
    if ((thrust_profile.arg_of_periapsis_change_thrust_profile) &&
        (thrust_profile.t_start_ <= t_) && (thrust_profile.t_end_ >= t_)) {
      std::array<double, 3> thrust_direction_vec = {sin(true_anomaly_), 0,
                                                    cos(true_anomaly_)};

      add_LVLH_thrust_profile(
          thrust_direction_vec,
          thrust_profile.sign_of_delta_omega * thrust_profile.thrust_magnitude_,
          t_, thrust_profile.t_end_);
      // can add one with t_end stretching all the way to the t_end of the total
      // arg of periapsis change maneuver because this temporary thrust profile
      // will be deleted after this RK45 step. This end time is given so that
      // whatever timestep is chosen for this step, all intermediate calculated
      // steps will also see this thrust.
      return 1;
    }
  }
  return 0;
}
void Satellite::check_for_maneuvers_to_initialize() {
  for (size_t maneuver_index = 0;
       maneuver_index < maneuvers_awaiting_initiation_.size();
       maneuver_index++) {
    std::pair<std::string, std::array<double, 3>> maneuver_info =
        maneuvers_awaiting_initiation_.at(maneuver_index);
    const double maneuver_start_time = maneuver_info.second.at(0);
    if ((maneuver_start_time <= t_) &&
        (maneuver_info.first == "Argument of Periapsis Change")) {
      // Need to initialize this maneuver
      const double target_arg_of_periapsis_deg = maneuver_info.second.at(1);
      const double target_arg_of_periapsis_rad =
          target_arg_of_periapsis_deg * (M_PI / 180.0);
      const double thrust_magnitude = maneuver_info.second.at(2);
      ThrustProfileLVLH thrust_profile(
          maneuver_start_time, target_arg_of_periapsis_rad, thrust_magnitude,
          arg_of_periapsis_, eccentricity_, a_, m_);
      thrust_profile_list_.push_back(thrust_profile);
      maneuvers_awaiting_initiation_.erase(
          maneuvers_awaiting_initiation_.begin() + maneuver_index);
    }
  }
}

void Satellite::add_maneuver(const std::string maneuver_type,
                             const double maneuver_start_time,
                             const double final_parameter_val,
                             const double thrust_magnitude) {
  std::array<double, 3> maneuver_vals = {maneuver_start_time,
                                         final_parameter_val, thrust_magnitude};
  std::pair<std::string, std::array<double, 3>> maneuver_info = {maneuver_type,
                                                                 maneuver_vals};
  maneuvers_awaiting_initiation_.push_back(maneuver_info);
}

double Satellite::get_total_energy() {
  // First, gravitational potential energy
  double orbital_radius = get_radius();
  double gravitational_potential_energy =
      -G * mass_Earth * m_ / orbital_radius;

  // Next, translational kinetic energy
  double orbital_speed = get_speed();
  double translational_kinetic_energy = (1.0 / 2.0) * m_ * (orbital_speed * orbital_speed);

  // Next, rotational kinetic energy
  double rotational_kinetic_energy = 0;
  std::array<double,3> body_angular_velocity_vec_wrt_LVLH_in_body_frame = body_angular_velocity_vec_wrt_LVLH_in_body_frame_;
  rotational_kinetic_energy += (0.5 * J_11_ * body_angular_velocity_vec_wrt_LVLH_in_body_frame.at(0)*body_angular_velocity_vec_wrt_LVLH_in_body_frame.at(0));
  rotational_kinetic_energy += (0.5 * J_22_ * body_angular_velocity_vec_wrt_LVLH_in_body_frame.at(1)*body_angular_velocity_vec_wrt_LVLH_in_body_frame.at(1));
  rotational_kinetic_energy += (0.5 * J_33_ * body_angular_velocity_vec_wrt_LVLH_in_body_frame.at(2)*body_angular_velocity_vec_wrt_LVLH_in_body_frame.at(2));

  return (gravitational_potential_energy + translational_kinetic_energy + rotational_kinetic_energy);
}