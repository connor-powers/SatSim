#ifndef UTILS_HEADER
#define UTILS_HEADER

#include <Eigen/Dense>
#include <functional>
#include <iostream>

#include "Satellite.h"

using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::Vector3d;
using Eigen::Vector4d;

std::array<double, 3> calculate_orbital_acceleration(
    const std::array<double, 3> input_r_vec, const double input_spacecraft_mass,
    const std::vector<std::array<double, 3>> input_vec_of_force_vectors_in_ECI =
        {});
std::array<double, 3> calculate_orbital_acceleration(
    const std::array<double, 3> input_r_vec, const double input_spacecraft_mass,
    const std::vector<ThrustProfileLVLH> input_list_of_thrust_profiles_LVLH,
    const double input_evaluation_time,
    const std::array<double, 3> input_velocity_vec,
    const double input_inclination, const double input_arg_of_periapsis,
    const double input_true_anomaly, const double input_F_10,
    const double input_A_p, const double input_A_s,
    const double input_satellite_mass, const bool perturbation,
    const bool atmospheric_drag);

std::array<double, 6> RK4_deriv_function_orbit_position_and_velocity(
    const std::array<double, 6> input_position_and_velocity,
    const double input_spacecraft_mass,
    const std::vector<std::array<double, 3>> input_vec_of_force_vectors_in_ECI =
        {});

// template <int T>
// std::array<double, T> RK4_step(
//     const std::array<double, T> y_n, const double input_step_size,
//     std::function<std::array<double, T>(const std::array<double, T> input_y_vec,
//                                         const double input_spacecraft_mass,
//                                         const std::vector<std::array<double, 3>>
//                                             input_vec_of_force_vectors_in_ECI)>
//         input_derivative_function,
//     const double input_spacecraft_mass,
//     const std::vector<std::array<double, 3>>
//         input_vec_of_force_vectors_in_ECI_at_t = {},
//     const std::vector<std::array<double, 3>>
//         input_vec_of_force_vectors_in_ECI_at_t_and_halfstep = {},
//     const std::vector<std::array<double, 3>>
//         input_vec_of_force_vectors_in_ECI_at_t_and_step = {}) {
//   // ref:
//   // https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#The_Runge%E2%80%93Kutta_method
//   std::array<double, T> y_nplus1 = y_n;

//   // first, k=1;
//   // going to be assuming derivative function does not have explicit time
//   // dependence
//   std::array<double, T> k_1 = input_derivative_function(
//       y_n, input_spacecraft_mass, input_vec_of_force_vectors_in_ECI_at_t);

//   // now k_2
//   std::array<double, T> y_vec_for_k_2 = y_n;
//   for (size_t ind = 0; ind < T; ind++) {
//     y_vec_for_k_2.at(ind) += ((input_step_size / 2) * k_1.at(ind));
//   }

//   std::array<double, T> k_2 = input_derivative_function(
//       y_vec_for_k_2, input_spacecraft_mass,
//       input_vec_of_force_vectors_in_ECI_at_t_and_halfstep);
//   // now k_3
//   std::array<double, T> y_vec_for_k_3 = y_n;
//   for (size_t ind = 0; ind < T; ind++) {
//     y_vec_for_k_3.at(ind) += ((input_step_size / 2) * k_2.at(ind));
//   }

//   std::array<double, T> k_3 = input_derivative_function(
//       y_vec_for_k_3, input_spacecraft_mass,
//       input_vec_of_force_vectors_in_ECI_at_t_and_halfstep);

//   // now k_4
//   std::array<double, T> y_vec_for_k_4 = y_n;
//   for (size_t ind = 0; ind < T; ind++) {
//     y_vec_for_k_4.at(ind) += (input_step_size * k_3.at(ind));
//   }

//   std::array<double, T> k_4 = input_derivative_function(
//       y_vec_for_k_4, input_spacecraft_mass,
//       input_vec_of_force_vectors_in_ECI_at_t_and_step);

//   for (size_t ind = 0; ind < T; ind++) {
//     y_nplus1.at(ind) +=
//         ((input_step_size / 6) *
//          (k_1.at(ind) + 2 * k_2.at(ind) + 2 * k_3.at(ind) + k_4.at(ind)));
//   }

//   return y_nplus1;
// }

void sim_and_draw_orbit_gnuplot(
    std::vector<Satellite> input_satellite_vector, const double input_timestep,
    const double input_total_sim_time, const double input_epsilon,
    const bool perturbation = true, const bool atmospheric_drag = false,
    const std::pair<double, double> drag_elements = {},
    const std::string input_terminal = "qt", //Currently, "qt" and "png" are supported
    const std::string output_file_name = "output");

template <int T>
std::pair<std::array<double, T>, std::pair<double, double>> RK45_step(
    std::array<double, T> y_n, double input_step_size,
    std::function<std::array<double, T>(
        const std::array<double, T>, const double,
        std::vector<ThrustProfileLVLH>, double, double, double, double, bool)>
        input_derivative_function,
    const double input_spacecraft_mass,
    std::vector<ThrustProfileLVLH> input_list_of_thrust_profiles_LVLH,
    double input_t_n, double input_epsilon, double input_inclination,
    double input_arg_of_periapsis, double input_true_anomaly,
    bool perturbation) {
  // Version for satellite orbital motion time evolution
  // Implementing RK4(5) method for its adaptive step size
  // Refs:https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method
  // ,
  // https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#The_Runge%E2%80%93Kutta_method
  std::array<double, 6> nodes = {0.0, 1.0 / 4, 3.0 / 8, 12.0 / 13,
                                 1.0, 1.0 / 2};  // c coefficients
  std::array<double, 6> CH_vec = {16.0 / 135,      0.0,       6656.0 / 12825,
                                  28561.0 / 56430, -9.0 / 50, 2.0 / 55};
  std::array<double, 6> CT_vec = {-1.0 / 360,     0.0,       128.0 / 4275,
                                  2197.0 / 75240, -1.0 / 50, -2.0 / 55};
  int s = 6;

  MatrixXd RK_matrix = MatrixXd::Zero(s, s);
  RK_matrix(1, 0) = 1.0 / 4;
  RK_matrix(2, 0) = 3.0 / 32;
  RK_matrix(2, 1) = 9.0 / 32;
  RK_matrix(3, 0) = 1932.0 / 2197;
  RK_matrix(3, 1) = -7200.0 / 2197;
  RK_matrix(3, 2) = 7296.0 / 2197;
  RK_matrix(4, 0) = 439.0 / 216;
  RK_matrix(4, 1) = -8.0;
  RK_matrix(4, 2) = 3680.0 / 513;
  RK_matrix(4, 3) = -845.0 / 4104;
  RK_matrix(5, 0) = -8.0 / 27;
  RK_matrix(5, 1) = 2.0;
  RK_matrix(5, 2) = -3544.0 / 2565;
  RK_matrix(5, 3) = 1859.0 / 4104;
  RK_matrix(5, 4) = -11.0 / 40;

  std::vector<std::array<double, T>> k_vec_vec = {};

  // Need to populate k_vec (compute vector k_i for i=0...s-1)
  for (size_t k_ind = 0; k_ind < s; k_ind++) {
    double evaluation_time = input_t_n + (nodes.at(k_ind) * input_step_size);
    std::array<double, T> y_n_evaluated_value = y_n;
    std::array<double, T> k_vec_at_this_s;
    for (size_t s_ind = 0; s_ind < k_ind; s_ind++) {
      for (size_t y_val_ind = 0; y_val_ind < y_n.size(); y_val_ind++) {
        y_n_evaluated_value.at(y_val_ind) += RK_matrix(k_ind, s_ind) *
                                             k_vec_vec.at(s_ind).at(y_val_ind);
      }
    }
    std::array<double, 6> derivative_function_output =
        input_derivative_function(y_n_evaluated_value, input_spacecraft_mass,
                                  input_list_of_thrust_profiles_LVLH,
                                  evaluation_time, input_inclination,
                                  input_arg_of_periapsis, input_true_anomaly,
                                  perturbation);
    for (size_t y_val_ind = 0; y_val_ind < y_n.size(); y_val_ind++) {
      k_vec_at_this_s.at(y_val_ind) =
          input_step_size * derivative_function_output.at(y_val_ind);
    }
    k_vec_vec.push_back(k_vec_at_this_s);
  }

  std::array<double, T> y_nplusone = y_n;

  for (size_t s_ind = 0; s_ind < s; s_ind++) {
    for (size_t y_ind = 0; y_ind < y_n.size(); y_ind++) {
      double tmp = CH_vec.at(s_ind) * k_vec_vec.at(s_ind).at(y_ind);
      y_nplusone.at(y_ind) += tmp;
    }
  }

  std::array<double, T> TE_vec;
  for (size_t s_ind = 0; s_ind < s; s_ind++) {
    for (size_t y_ind = 0; y_ind < std::size(TE_vec); y_ind++) {
      TE_vec.at(y_ind) += CT_vec.at(s_ind) * k_vec_vec.at(s_ind).at(y_ind);
    }
  }

  for (size_t y_ind = 0; y_ind < std::size(TE_vec); y_ind++) {
    TE_vec.at(y_ind) = abs(TE_vec.at(y_ind));
  }
  // I'm going to use the max TE found across the whole position+velocity vec as
  // the TE in the calculation of the next stepsize
  double max_TE = TE_vec.at(
      std::distance(std::begin(TE_vec),
                    std::max_element(std::begin(TE_vec), std::end(TE_vec))));
  double epsilon_ratio = input_epsilon / max_TE;

  double h_new = 0.9 * input_step_size * std::pow(epsilon_ratio, 1.0 / 5);

  if (max_TE <= input_epsilon) {
    std::pair<double, double> output_timestep_pair;
    output_timestep_pair.first =
        input_step_size;  // First timestep size in this pair is the one
                          // successfully used in this calculation
    output_timestep_pair.second =
        h_new;  // Second timestep size in this pair is the one to be used in
                // the next step
    std::pair<std::array<double, T>, std::pair<double, double>> output_pair;
    output_pair.first = y_nplusone;
    output_pair.second = output_timestep_pair;
    return output_pair;
  } else {
    return RK45_step<T>(
        y_n, h_new, input_derivative_function, input_spacecraft_mass,
        input_list_of_thrust_profiles_LVLH, input_t_n, input_epsilon,
        input_inclination, input_arg_of_periapsis, input_true_anomaly,
        perturbation);
  }
}

std::array<double, 3> convert_LVLH_to_ECI_manual(
    const std::array<double, 3> input_LVLH_vec,
    const std::array<double, 3> input_position_vec,
    const std::array<double, 3> input_velocity_vec);

std::array<double, 6> RK45_deriv_function_orbit_position_and_velocity(
    const std::array<double, 6> input_position_and_velocity,
    const double input_spacecraft_mass,
    const std::vector<ThrustProfileLVLH> input_list_of_thrust_profiles_LVLH,
    const double input_evaluation_time, const double input_inclination,
    const double input_arg_of_periapsis, const double input_true_anomaly,
    const double input_F_10, const double input_A_p, const double input_A_s,
    const double input_satellite_mass, const bool perturbation,
    const bool atmospheric_drag);

std::array<double, 3> convert_cylindrical_to_cartesian(
    const double input_r_comp, const double input_theta_comp,
    const double input_z_comp, const double input_theta);
void sim_and_plot_orbital_elem_gnuplot(
    std::vector<Satellite> input_satellite_vector, const double input_timestep,
    const double input_total_sim_time, const double input_epsilon,
    const std::string input_orbital_element_name,
    const std::string file_name = "output",
    const bool perturbation = true, const bool atmospheric_drag = false,
    const std::pair<double, double> drag_elements = {});
void sim_and_plot_attitude_evolution_gnuplot(
    std::vector<Satellite> input_satellite_vector, const double input_timestep,
    const double input_total_sim_time, const double input_epsilon,
    const std::string input_plotted_val_name, 
    const std::string file_name = "output",
    const bool perturbation = true,
    const bool atmospheric_drag = false,
    const std::pair<double, double> drag_elements = {});

Matrix3d rollyawpitch_bodyframe_to_LVLH(
    const std::array<double, 3> input_bodyframe_vec, const double input_roll,
    const double input_pitch, const double input_yaw);
std::array<double, 4> rollyawpitch_angles_to_quaternion(
    const double input_roll, const double input_pitch, const double input_yaw);
std::array<double, 4> rollpitchyaw_angles_to_quaternion(
    const double input_roll, const double input_pitch, const double input_yaw);

Matrix3d LVLH_to_body_transformation_matrix_from_quaternion(
    const std::array<double, 4> input_bodyframe_quaternion_relative_to_LVLH);
Vector4d quaternion_multiplication(const Vector4d quaternion_1,
                                   const Vector4d quaternion_2);
Vector4d quaternion_kinematics_equation(
    const Vector4d quaternion_of_bodyframe_relative_to_ref_frame,
    const Vector3d angular_velocity_vec_wrt_ref_frame_in_body_frame);
std::array<double, 7> RK45_satellite_body_angular_deriv_function(
    const std::array<double, 7> combined_bodyframe_angular_array,
    const Matrix3d J_matrix,
    const std::vector<BodyframeTorqueProfile>
        input_bodyframe_torque_profile_list,
    const Vector3d input_omega_I,
    const double input_orbital_angular_acceleration,
    const Matrix3d input_LVLH_to_bodyframe_transformation_matrix,
    const Vector3d input_omega_LVLH_wrt_inertial_in_LVLH,
    const double input_evaluation_time);
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
    const bool perturbation, const bool atmospheric_drag);

template <int T>
std::pair<std::array<double, T>, std::pair<double, double>> RK45_step(
    const std::array<double, T> y_n, const double input_step_size,
    std::function<std::array<double, T>(
        const std::array<double, T>, const Matrix3d,
        const std::vector<BodyframeTorqueProfile>, const Vector3d, const double,
        const Matrix3d, const Vector3d, const double,
        const std::vector<ThrustProfileLVLH>, const double, const double,
        const double, const double, const double, const double, const double,
        const bool, const bool)>
        input_combined_derivative_function,
    const Matrix3d J_matrix,
    const std::vector<BodyframeTorqueProfile>
        input_bodyframe_torque_profile_list,
    const Vector3d input_omega_I, double input_orbital_angular_acceleration,
    const Matrix3d input_LVLH_to_bodyframe_transformation_matrix,
    const Vector3d input_omega_LVLH_wrt_inertial_in_LVLH,
    const double input_spacecraft_mass,
    const std::vector<ThrustProfileLVLH> input_list_of_thrust_profiles_LVLH,
    const double input_inclination, const double input_arg_of_periapsis,
    const double input_true_anomaly, const double input_F_10,
    const double input_A_p, const double input_A_s, const bool perturbation,
    const double atmospheric_drag, const double input_t_n,
    const double input_epsilon) {
  // Version for combined satellite orbital motion and attitude time evolution
  // Implementing RK4(5) method for its adaptive step size
  // Refs:https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method
  // ,
  // https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#The_Runge%E2%80%93Kutta_method

  std::array<double, 6> nodes = {0.0, 1.0 / 4, 3.0 / 8, 12.0 / 13,
                                 1.0, 1.0 / 2};  // c coefficients
  std::array<double, 6> CH_vec = {16.0 / 135,      0.0,       6656.0 / 12825,
                                  28561.0 / 56430, -9.0 / 50, 2.0 / 55};
  std::array<double, 6> CT_vec = {-1.0 / 360,     0.0,       128.0 / 4275,
                                  2197.0 / 75240, -1.0 / 50, -2.0 / 55};
  int s = 6;

  MatrixXd RK_matrix = MatrixXd::Zero(s, s);
  RK_matrix(1, 0) = 1.0 / 4;
  RK_matrix(2, 0) = 3.0 / 32;
  RK_matrix(2, 1) = 9.0 / 32;
  RK_matrix(3, 0) = 1932.0 / 2197;
  RK_matrix(3, 1) = -7200.0 / 2197;
  RK_matrix(3, 2) = 7296.0 / 2197;
  RK_matrix(4, 0) = 439.0 / 216;
  RK_matrix(4, 1) = -8.0;
  RK_matrix(4, 2) = 3680.0 / 513;
  RK_matrix(4, 3) = -845.0 / 4104;
  RK_matrix(5, 0) = -8.0 / 27;
  RK_matrix(5, 1) = 2.0;
  RK_matrix(5, 2) = -3544.0 / 2565;
  RK_matrix(5, 3) = 1859.0 / 4104;
  RK_matrix(5, 4) = -11.0 / 40;

  std::vector<std::array<double, T>> k_vec_vec = {};

  // Need to populate k_vec (compute vector k_i for i=0...s-1)
  for (size_t k_ind = 0; k_ind < s; k_ind++) {
    double evaluation_time = input_t_n + (nodes.at(k_ind) * input_step_size);
    std::array<double, T> y_n_evaluated_value = y_n;
    std::array<double, T> k_vec_at_this_s;
    for (size_t s_ind = 0; s_ind < k_ind; s_ind++) {
      for (size_t y_val_ind = 0; y_val_ind < y_n.size(); y_val_ind++) {
        y_n_evaluated_value.at(y_val_ind) += input_step_size *
                                             RK_matrix(k_ind, s_ind) *
                                             k_vec_vec.at(s_ind).at(y_val_ind);
      }
    }
    std::array<double, T> derivative_function_output =
        input_combined_derivative_function(
            y_n_evaluated_value, J_matrix, input_bodyframe_torque_profile_list,
            input_omega_I, input_orbital_angular_acceleration,
            input_LVLH_to_bodyframe_transformation_matrix,
            input_omega_LVLH_wrt_inertial_in_LVLH, input_spacecraft_mass,
            input_list_of_thrust_profiles_LVLH, evaluation_time,
            input_inclination, input_arg_of_periapsis, input_true_anomaly,
            input_F_10, input_A_p, input_A_s, perturbation, atmospheric_drag);
    for (size_t y_val_ind = 0; y_val_ind < y_n.size(); y_val_ind++) {
      k_vec_at_this_s.at(y_val_ind) =
          input_step_size * derivative_function_output.at(y_val_ind);
    }

    k_vec_vec.push_back(k_vec_at_this_s);
  }

  std::array<double, T> y_nplusone = y_n;

  for (size_t s_ind = 0; s_ind < s; s_ind++) {
    for (size_t y_ind = 0; y_ind < y_n.size(); y_ind++) {
      double tmp = CH_vec.at(s_ind) * k_vec_vec.at(s_ind).at(y_ind);
      y_nplusone.at(y_ind) += tmp;
    }
  }

  std::array<double, T> TE_vec = {0};
  for (size_t y_ind = 0; y_ind < y_n.size(); y_ind++) {
    for (size_t s_ind = 0; s_ind < s; s_ind++) {
      TE_vec.at(y_ind) += CT_vec.at(s_ind) * k_vec_vec.at(s_ind).at(y_ind);
    }
  }

  for (size_t ind = 0; ind < TE_vec.size(); ind++) {
    TE_vec.at(ind) = abs(TE_vec.at(ind));
  }

  // I'm going to use the max TE found across the whole position+velocity vec as
  // the TE in the calculation of the next stepsize
  double max_TE = TE_vec.at(
      std::distance(std::begin(TE_vec),
                    std::max_element(std::begin(TE_vec), std::end(TE_vec))));
  double epsilon_ratio = input_epsilon / max_TE;

  double h_new = 0.9 * input_step_size * std::pow(epsilon_ratio, 1.0 / 5);

  if (max_TE <= input_epsilon) {
    std::pair<double, double> output_timestep_pair;
    output_timestep_pair.first =
        input_step_size;  // First timestep size in this pair is the one
                          // successfully used in this calculation
    output_timestep_pair.second =
        h_new;  // Second timestep size in this pair is the one to be used in
                // the next step
    std::pair<std::array<double, T>, std::pair<double, double>> output_pair;
    output_pair.first = y_nplusone;
    output_pair.second = output_timestep_pair;
    return output_pair;
  } else {
    return RK45_step<T>(
        y_n, h_new, input_combined_derivative_function, J_matrix,
        input_bodyframe_torque_profile_list, input_omega_I,
        input_orbital_angular_acceleration,
        input_LVLH_to_bodyframe_transformation_matrix,
        input_omega_LVLH_wrt_inertial_in_LVLH, input_spacecraft_mass,
        input_list_of_thrust_profiles_LVLH, input_inclination,
        input_arg_of_periapsis, input_true_anomaly, input_F_10, input_A_p,
        input_A_s, perturbation, atmospheric_drag, input_t_n, input_epsilon);
  }
}
Vector3d calculate_omega_I(
    const Vector3d input_bodyframe_ang_vel_vector_wrt_lvlh,
    const Matrix3d input_LVLH_to_bodyframe_transformation_matrix,
    const double input_orbital_rate);

Matrix3d construct_J_matrix(const double input_Jxx, const double input_Jyy,
                            const double input_Jzz);

std::array<double, 3> convert_quaternion_to_roll_yaw_pitch_angles(
    const std::array<double, 4>);

std::array<double, 4> normalize_quaternion(
    std::array<double, 4> input_quaternion);

std::array<double, 3> convert_array_from_LVLH_to_bodyframe(
    const std::array<double, 3> input_LVLH_frame_array, const double input_roll,
    const double input_yaw, const double input_pitch);


std::array<double,3>  convert_lat_long_to_ECEF(double latitude, double longitude, double height);
std::array<double,3>  convert_ECEF_to_ECI(std::array<double,3> input_ECEF_position);
#endif