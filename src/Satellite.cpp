
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

// implementations of methods for Satellite class
double Satellite::calculate_orbital_period(double input_semimajor_axis) {
  // https://en.wikipedia.org/wiki/Orbital_period
  double T = 2 * M_PI * sqrt(pow(input_semimajor_axis, 3) / (G * mass_Earth));
  return T;
}

std::pair<double, double> Satellite::calculate_eccentric_anomaly(
    const double input_eccentricity, const double input_true_anomaly,
    const double input_semimajor_axis) {
  // https://en.wikipedia.org/wiki/Eccentric_anomaly
  double numerator =
      std::sqrt(1 - pow(input_eccentricity, 2)) * sin(input_true_anomaly);
  double denominator = input_eccentricity + cos(input_true_anomaly);
  double eccentric_anomaly = atan2(numerator, denominator);

  std::pair<double, double> output_pair;
  output_pair.first = eccentric_anomaly;
  double computed_distance_from_Earth =
      input_semimajor_axis * (1 - input_eccentricity * cos(eccentric_anomaly));
  output_pair.second = computed_distance_from_Earth;
  return output_pair;
}

std::array<double, 3> Satellite::calculate_perifocal_position() {
  // Using approach from Fundamentals of Astrodynamics
  std::array<double, 3> calculated_perifocal_position;
  double mu = G * mass_Earth;

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

std::array<double, 3> Satellite::calculate_perifocal_velocity() {
  // Using approach from Fundamentals of Astrodynamics
  std::array<double, 3> calculated_perifocal_velocity;
  double mu = G * mass_Earth;
  double p =
      a_ * (1 - eccentricity_) * (1 + eccentricity_);  // rearranging eq. 1-47

  double v_perifocal_P = sqrt(mu / p) * (-sin(true_anomaly_));
  double v_perifocal_Q = sqrt(mu / p) * (eccentricity_ + cos(true_anomaly_));

  calculated_perifocal_velocity.at(0) = v_perifocal_P;
  calculated_perifocal_velocity.at(1) = v_perifocal_Q;
  calculated_perifocal_velocity.at(2) = 0;  // 0 component in the W direction

  return calculated_perifocal_velocity;
}

std::array<double, 3> Satellite::convert_perifocal_to_ECI(
    std::array<double, 3> input_perifocal_vec) {
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

std::array<double, 3> Satellite::convert_ECI_to_perifocal(
    std::array<double, 3> input_ECI_vec) {
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

void Satellite::evolve_RK4(double input_step_size) {
  // format input position and velocity arrays into single array for RK4 step
  std::array<double, 6> combined_initial_position_and_velocity_array = {};
  std::pair<std::array<double, 3>, std::array<double, 3>>
      output_position_velocity_pair = {};

  for (size_t ind = 0; ind < 3; ind++) {
    combined_initial_position_and_velocity_array.at(ind) =
        ECI_position_.at(ind);
  }
  for (size_t ind = 3; ind < 6; ind++) {
    combined_initial_position_and_velocity_array.at(ind) =
        ECI_velocity_.at(ind - 3);
  }

  // populate list of thrust forces at half a timestep past, for RK4
  // calculations
  std::vector<std::array<double, 3>> list_of_LVLH_forces_at_half_timestep_past =
      {};
  std::vector<std::array<double, 3>> list_of_ECI_forces_at_half_timestep_past =
      {};

  for (ThrustProfileLVLH thrust_profile : thrust_profile_list_) {
    if (((t_ + (input_step_size / 2)) >= thrust_profile.t_start_) &&
        ((t_ + (input_step_size / 2)) <= thrust_profile.t_end_)) {
      list_of_LVLH_forces_at_half_timestep_past.push_back(
          thrust_profile.LVLH_force_vec_);
      std::array<double, 3> ECI_thrust_vector = convert_LVLH_to_ECI_manual(
          thrust_profile.LVLH_force_vec_, ECI_position_, ECI_velocity_);
      list_of_ECI_forces_at_half_timestep_past.push_back(ECI_thrust_vector);
    }
  }

  // populate list of thrust forces at a timestep past, for RK4 calculations
  std::vector<std::array<double, 3>> list_of_LVLH_forces_at_one_timestep_past =
      {};
  std::vector<std::array<double, 3>> list_of_ECI_forces_at_one_timestep_past =
      {};

  for (ThrustProfileLVLH thrust_profile : thrust_profile_list_) {
    if (((t_ + input_step_size) >= thrust_profile.t_start_) &&
        ((t_ + input_step_size) <= thrust_profile.t_end_)) {
      list_of_LVLH_forces_at_one_timestep_past.push_back(
          thrust_profile.LVLH_force_vec_);
      std::array<double, 3> ECI_thrust_vector = convert_LVLH_to_ECI_manual(
          thrust_profile.LVLH_force_vec_, ECI_position_, ECI_velocity_);
      list_of_ECI_forces_at_one_timestep_past.push_back(ECI_thrust_vector);
    }
  }

  std::array<double, 6> output_combined_position_and_velocity_array =
      RK4_step<6>(combined_initial_position_and_velocity_array, input_step_size,
                  RK4_deriv_function_orbit_position_and_velocity, m_,
                  list_of_ECI_forces_at_this_time_,
                  list_of_ECI_forces_at_half_timestep_past,
                  list_of_ECI_forces_at_one_timestep_past);
  // std::array<double,6> output_combined_angular_array=
  // RK4_step<6>(combined_initial_angular_array,input_step_size,RK4_deriv_function_angular,I_,list_of_body_frame_torques_at_this_time_,list_of_body_frame_torques_at_half_timestep_past,list_of_body_frame_torques_at_one_timestep_past);

  for (size_t ind = 0; ind < 3; ind++) {
    ECI_position_.at(ind) = output_combined_position_and_velocity_array.at(ind);
    ECI_velocity_.at(ind) =
        output_combined_position_and_velocity_array.at(ind + 3);
    // also update the perifocal versions
    perifocal_position_ = convert_ECI_to_perifocal(ECI_position_);
    perifocal_velocity_ = convert_ECI_to_perifocal(ECI_velocity_);
  }
  t_ += input_step_size;

  list_of_LVLH_forces_at_this_time_ = list_of_LVLH_forces_at_one_timestep_past;
  list_of_ECI_forces_at_this_time_ = list_of_ECI_forces_at_one_timestep_past;

  // Update orbital parameters
  update_orbital_elements_from_position_and_velocity();

  return;
}

// std::array<double,3> Satellite::convert_LVLH_to_ECI(std::array<double,3>
// input_LVLH_vec){
//     Vector3d input_LVLH_vec_eigen;
//     input_LVLH_vec_eigen <<
//     input_LVLH_vec.at(0),input_LVLH_vec.at(1),input_LVLH_vec.at(2);
//     //ref:
//     https://ntrs.nasa.gov/api/citations/20205003902/downloads/Introduction%20to%20Orbital%20Mechanics%20and%20Spacecraft%20Attitudes%20for%20Thermal%20Engineers%20CHARTS%20PDF.pdf
//     Matrix3d ref_mat;
//     ref_mat << 0,0,-1,
//                1,0,0,
//                0,-1,0;

//     Matrix3d omega_mat=z_rot_matrix(raan_);
//     Matrix3d inclination_mat=y_rot_matrix(inclination_);
//     Matrix3d arg_periapsis_mat=z_rot_matrix(arg_of_periapsis_);
//     Matrix3d true_anomaly_mat=z_rot_matrix(true_anomaly_);

//     Matrix3d pitch_angle_mat=y_rot_matrix(pitch_angle_);
//     Matrix3d yaw_angle_mat=z_rot_matrix(yaw_angle_);
//     Matrix3d roll_angle_mat=x_rot_matrix(roll_angle_);

//     Vector3d
//     ECI_vec_eigen=omega_mat*inclination_mat*arg_periapsis_mat*true_anomaly_mat*ref_mat*pitch_angle_mat*yaw_angle_mat*roll_angle_mat*input_LVLH_vec_eigen;

//     std::array<double,3>
//     ECI_vec_output={ECI_vec_eigen(0),ECI_vec_eigen(1),ECI_vec_eigen(2)};
//     return ECI_vec_output;
// }

void Satellite::add_LVLH_thrust_profile(
    std::array<double, 3> input_LVLH_thrust_vector,
    double input_thrust_start_time, double input_thrust_end_time) {
  ThrustProfileLVLH new_thrust_profile(
      input_thrust_start_time, input_thrust_end_time, input_LVLH_thrust_vector);
  thrust_profile_list_.push_back(new_thrust_profile);
  if (input_thrust_start_time == 0) {
    list_of_LVLH_forces_at_this_time_.push_back(input_LVLH_thrust_vector);
    std::array<double, 3> ECI_thrust_vector = convert_LVLH_to_ECI_manual(
        input_LVLH_thrust_vector, ECI_position_, ECI_velocity_);
    list_of_ECI_forces_at_this_time_.push_back(ECI_thrust_vector);
  }
}

void Satellite::add_LVLH_thrust_profile(
    std::array<double, 3> input_LVLH_normalized_thrust_direction,
    double input_LVLH_thrust_magnitude, double input_thrust_start_time,
    double input_thrust_end_time) {
  ThrustProfileLVLH new_thrust_profile(
      input_thrust_start_time, input_thrust_end_time,
      input_LVLH_normalized_thrust_direction, input_LVLH_thrust_magnitude);
  thrust_profile_list_.push_back(new_thrust_profile);

  std::array<double, 3> LVLH_thrust_vec = {0, 0, 0};

  for (size_t ind = 0; ind < 3; ind++) {
    LVLH_thrust_vec.at(ind) = input_LVLH_thrust_magnitude *
                              input_LVLH_normalized_thrust_direction.at(ind);
  }

  if (input_thrust_start_time == 0) {
    list_of_LVLH_forces_at_this_time_.push_back(LVLH_thrust_vec);
    std::array<double, 3> ECI_thrust_vector = convert_LVLH_to_ECI_manual(
        LVLH_thrust_vec, ECI_position_, ECI_velocity_);
    list_of_ECI_forces_at_this_time_.push_back(ECI_thrust_vector);
  }
}

// std::array<double,3>
// Satellite::convert_body_frame_to_LVLH(std::array<double,3>
// input_body_frame_vec){
//     //taking LVLH to be the "unrotated" body frame

//     std::array<double,3>
//     LVLH_vec=convert_rotated_body_frame_to_unrotated_body_frame(input_body_frame_vec,
//     theta_, phi_, psi_); return LVLH_vec;
// }

// std::array<double,3>
// Satellite::convert_body_frame_to_ECI(std::array<double,3>
// input_body_frame_vec){
//     //Current strategy: first convert rotated body frame to unrotated body
//     frame, which I'm taking to be the LVLH frame
//     //Then convert the unrotated LVLH frame into ECI
//     std::array<double,3>
//     LVLH_vec=convert_body_frame_to_LVLH(input_body_frame_vec);
//     std::array<double,3> ECI_vec=convert_LVLH_to_ECI_manual(LVLH_vec);
// }

int Satellite::update_orbital_elements_from_position_and_velocity() {
  // Anytime the orbit is changed via external forces, need to update the
  // orbital parameters of the satellite. True anomaly should change over time
  // even in absence of external forces Using approach from Fundamentals of
  // Astrodynamics
  int error_code = 0;  // 0 represents nominal operation
  double mu = G * mass_Earth;

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
      (1.0 / mu) * (v_magnitude * v_magnitude - (mu / r_magnitude)) *
      position_vector;
  Vector3d e_vec_component_2 =
      (1.0 / mu) * (position_vector.dot(velocity_vector)) * velocity_vector;
  Vector3d e_vec = e_vec_component_1 - e_vec_component_2;
  double calculated_eccentricity = e_vec.norm();

  double calculated_p = h * h / mu;

  double calculated_i = acos(h_vector(2) / h);

  double calculated_RAAN = acos(n_vector(0) / n);
  if (n_vector(1) < 0) {
    calculated_RAAN = 2 * M_PI - calculated_RAAN;
  }
  if (std::isnan(calculated_RAAN)) {
    std::cout
        << "Calculated RAAN was undefined. One possible cause of this is an "
           "orbit with zero inclination. Current magnitude of line of nodes: "
        << n << "\n";
    error_code = 1;
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
    if (position_vector.dot(velocity_vector) < 0) {
      calculated_true_anomaly = 2 * M_PI - calculated_true_anomaly;
    }

    if (std::isnan(calculated_arg_of_periapsis)) {
      std::cout << "Calculated argument of periapsis was undefined. One "
                   "possible cause of this is an orbit with zero inclination. "
                   "Current magnitude of line of nodes: "
                << n << "\n";
      error_code = 1;
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
    if (std::isnan(calculated_arg_of_latitude)) {
      std::cout << "Calculated argument of latitude was undefined. One "
                   "possible cause of this is an orbit with zero inclination. "
                   "Current magnitude of line of nodes: "
                << n << "\n";
      error_code = 1;
    }
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

std::pair<double, int> Satellite::evolve_RK45(double input_epsilon,
                                              double input_step_size,
                                              bool perturbation) {
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
          inclination_, arg_of_periapsis_, true_anomaly_, perturbation, t_,
          input_epsilon);

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

  evolve_RK45_output_pair.first = new_step_size;
  evolve_RK45_output_pair.second = orbit_elems_error_code;

  return evolve_RK45_output_pair;
}

double Satellite::get_orbital_element(std::string orbital_element_name) {
  if (orbital_element_name == "Semimajor Axis") {
    return a_;
  } else if (orbital_element_name == "Eccentricity") {
    return eccentricity_;
  } else if (orbital_element_name == "Inclination") {
    return inclination_;
  } else if (orbital_element_name == "RAAN") {
    return raan_;
  } else if (orbital_element_name == "Argument of Periapsis") {
    return arg_of_periapsis_;
  } else if (orbital_element_name == "True Anomaly") {
    return true_anomaly_;
  } else if (orbital_element_name == "Orbital Rate") {
    return orbital_rate_;
  } else if (orbital_element_name == "Orbital Angular Acceleration") {
    return orbital_angular_acceleration_;
  } else {
    std::cout << "Unrecognized argument, returning -1";
    return -1;
  }
}

void Satellite::add_bodyframe_torque_profile(
    std::array<double, 3> input_bodyframe_torque_vector,
    double input_torque_start_time, double input_torque_end_time) {
  BodyframeTorqueProfile new_torque_profile(input_torque_start_time,
                                            input_torque_end_time,
                                            input_bodyframe_torque_vector);
  bodyframe_torque_profile_list_.push_back(new_torque_profile);
  if (input_torque_start_time == 0) {
    list_of_body_frame_torques_at_this_time_.push_back(
        input_bodyframe_torque_vector);
  }
}

void Satellite::add_bodyframe_torque_profile(
    std::array<double, 3> input_bodyframe_direction_unit_vec,
    double input_bodyframe_torque_magnitude, double input_torque_start_time,
    double input_torque_end_time) {
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

void Satellite::initialize_and_normalize_body_quaternion(double roll_angle,
                                                         double pitch_angle,
                                                         double yaw_angle) {
  std::array<double, 4> quaternion =
      rollyawpitch_angles_to_quaternion(roll_angle, pitch_angle, yaw_angle);
  std::array<double, 4> normalized_quaternion =
      normalize_quaternion(quaternion);
  quaternion_satellite_bodyframe_wrt_LVLH_ = normalized_quaternion;
  return;
}

double Satellite::get_attitude_val(std::string input_attitude_val_name) {
  if (input_attitude_val_name == "Roll") {
    return roll_angle_;
  } else if (input_attitude_val_name == "Pitch") {
    return pitch_angle_;
  } else if (input_attitude_val_name == "Yaw") {
    return yaw_angle_;
  } else if (input_attitude_val_name == "omega_x") {
    return body_angular_velocity_vec_wrt_LVLH_in_body_frame_.at(0);
  } else if (input_attitude_val_name == "omega_y") {
    return body_angular_velocity_vec_wrt_LVLH_in_body_frame_.at(1);
  } else if (input_attitude_val_name == "omega_z") {
    return body_angular_velocity_vec_wrt_LVLH_in_body_frame_.at(2);
  } else if (input_attitude_val_name == "q_0") {
    return quaternion_satellite_bodyframe_wrt_LVLH_.at(0);
  } else if (input_attitude_val_name == "q_1") {
    return quaternion_satellite_bodyframe_wrt_LVLH_.at(1);
  } else if (input_attitude_val_name == "q_2") {
    return quaternion_satellite_bodyframe_wrt_LVLH_.at(2);
  } else if (input_attitude_val_name == "q_3") {
    return quaternion_satellite_bodyframe_wrt_LVLH_.at(3);
  } else {
    std::cout << "Invalid attitude val choice " << input_attitude_val_name
              << ", returning -1\n";
    return -1.0;
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