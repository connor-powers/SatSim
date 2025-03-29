#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cmath>
#include <iostream>

#include "Satellite.h"

using Eigen::Matrix3d;
using Eigen::Matrix4d;
using Eigen::Quaterniond;
using Eigen::Vector3d;
using Eigen::Vector4d;

Vector4d quaternion_multiplication(Vector4d quaternion_1,
                                   Vector4d quaternion_2) {
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
    std::array<double, 3> input_LVLH_vec,
    std::array<double, 3> input_position_vec,
    std::array<double, 3> input_velocity_vec) {
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

std::array<double, 3> convert_ECI_to_LVLH_manual(
    std::array<double, 3> input_ECI_vec,
    std::array<double, 3> input_position_vec,
    std::array<double, 3> input_velocity_vec) {
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

  double v_x_LVLH = input_ECI_vec.at(0) * LVLH_x_unit_vec.dot(ECI_unit_vec_x) +
                    input_ECI_vec.at(1) * LVLH_x_unit_vec.dot(ECI_unit_vec_y) +
                    input_ECI_vec.at(2) * LVLH_x_unit_vec.dot(ECI_unit_vec_z);

  double v_y_LVLH = input_ECI_vec.at(0) * LVLH_y_unit_vec.dot(ECI_unit_vec_x) +
                    input_ECI_vec.at(1) * LVLH_y_unit_vec.dot(ECI_unit_vec_y) +
                    input_ECI_vec.at(2) * LVLH_y_unit_vec.dot(ECI_unit_vec_z);

  double v_z_LVLH = input_ECI_vec.at(0) * LVLH_z_unit_vec.dot(ECI_unit_vec_x) +
                    input_ECI_vec.at(1) * LVLH_z_unit_vec.dot(ECI_unit_vec_y) +
                    input_ECI_vec.at(2) * LVLH_z_unit_vec.dot(ECI_unit_vec_z);

  std::array<double, 3> output_LVLH_arr = {v_x_LVLH, v_y_LVLH, v_z_LVLH};
  return output_LVLH_arr;
}

// baselining cartesian coordinates in ECI frame

std::array<double, 3> calculate_orbital_acceleration(
    const std::array<double, 3> input_r_vec, const double input_spacecraft_mass,
    std::vector<std::array<double, 3>> input_vec_of_force_vectors_in_ECI) {
  // Note: this is the version used in the RK4 solver
  // orbital acceleration = -G m_Earth/distance^3 * r_vec (just based on
  // rearranging F=ma with a the acceleration due to gravitational attraction
  // between satellite and Earth
  // https://en.wikipedia.org/wiki/Newton%27s_law_of_universal_gravitation)
  // going to be assuming Earth's position doesn't change for now
  // also assuming Earth is spherical, can loosen this assumption in the future
  // note: this is in ECI frame

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

  // now add effects from externally-applied forces, e.g., thrusters, if any
  std::array<double, 3> acceleration_vec = acceleration_vec_due_to_gravity;

  for (std::array<double, 3> external_force_vec_in_ECI :
       input_vec_of_force_vectors_in_ECI) {
    acceleration_vec.at(0) +=
        (external_force_vec_in_ECI.at(0) / input_spacecraft_mass);
    acceleration_vec.at(1) +=
        (external_force_vec_in_ECI.at(1) / input_spacecraft_mass);
    acceleration_vec.at(2) +=
        (external_force_vec_in_ECI.at(2) / input_spacecraft_mass);
  }
  return acceleration_vec;
}

std::array<double, 3> convert_cylindrical_to_cartesian(double input_r_comp,
                                                       double input_theta_comp,
                                                       double input_z_comp,
                                                       double input_theta) {
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
    std::vector<ThrustProfileLVLH> input_list_of_thrust_profiles_LVLH,
    double input_evaluation_time,
    const std::array<double, 3> input_velocity_vec, double input_inclination,
    double input_arg_of_periapsis, double input_true_anomaly,
    bool perturbation) {
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

  for (ThrustProfileLVLH thrust_profile : input_list_of_thrust_profiles_LVLH) {
    if ((input_evaluation_time >= thrust_profile.t_start_) &&
        (input_evaluation_time <= thrust_profile.t_end_)) {
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
    double mu = G * mass_Earth;
    double C =
        3 * mu * J2 * radius_Earth * radius_Earth / (2 * pow(distance, 4));
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

  return acceleration_vec;
}

std::array<double, 6> RK4_deriv_function_orbit_position_and_velocity(
    std::array<double, 6> input_position_and_velocity,
    const double input_spacecraft_mass,
    std::vector<std::array<double, 3>> input_vec_of_force_vectors_in_ECI) {
  std::array<double, 6> derivative_of_input_y = {};
  std::array<double, 3> position_array = {};

  for (size_t ind = 0; ind < 3; ind++) {
    derivative_of_input_y.at(ind) = input_position_and_velocity.at(ind + 3);
    position_array.at(ind) = input_position_and_velocity.at(ind);
  }

  std::array<double, 3> calculated_orbital_acceleration =
      calculate_orbital_acceleration(position_array, input_spacecraft_mass,
                                     input_vec_of_force_vectors_in_ECI);

  for (size_t ind = 3; ind < 6; ind++) {
    derivative_of_input_y.at(ind) = calculated_orbital_acceleration.at(ind - 3);
  }

  return derivative_of_input_y;
}

std::array<double, 6> RK45_deriv_function_orbit_position_and_velocity(
    std::array<double, 6> input_position_and_velocity,
    const double input_spacecraft_mass,
    std::vector<ThrustProfileLVLH> input_list_of_thrust_profiles_LVLH,
    double input_evaluation_time, double input_inclination,
    double input_arg_of_periapsis, double input_true_anomaly,
    bool perturbation) {
  std::array<double, 6> derivative_of_input_y = {};
  std::array<double, 3> position_array = {};
  std::array<double, 3> velocity_array = {};

  for (size_t ind = 0; ind < 3; ind++) {
    derivative_of_input_y.at(ind) = input_position_and_velocity.at(ind + 3);
    velocity_array.at(ind) = input_position_and_velocity.at(ind + 3);
    position_array.at(ind) = input_position_and_velocity.at(ind);
  }

  std::array<double, 3> calculated_orbital_acceleration =
      calculate_orbital_acceleration(position_array, input_spacecraft_mass,
                                     input_list_of_thrust_profiles_LVLH,
                                     input_evaluation_time, velocity_array,
                                     input_inclination, input_arg_of_periapsis,
                                     input_true_anomaly, perturbation);

  for (size_t ind = 3; ind < 6; ind++) {
    derivative_of_input_y.at(ind) = calculated_orbital_acceleration.at(ind - 3);
  }

  return derivative_of_input_y;
}

void sim_and_draw_orbit_gnuplot(std::vector<Satellite> input_satellite_vector,
                                double input_timestep,
                                double input_total_sim_time,
                                double input_epsilon, bool perturbation) {
  if (input_satellite_vector.size() < 1) {
    std::cout << "No input Satellite objects\n";
    return;
  }

  // first, open "pipe" to gnuplot
  FILE *gnuplot_pipe = popen("gnuplot -persist", "w");
  // if it exists
  if (gnuplot_pipe) {
    fprintf(gnuplot_pipe, "set terminal qt size 900,700 font 'Helvetica,14'\n");
    // formatting
    fprintf(gnuplot_pipe, "set xlabel 'x [m]' offset 0,-2\n");
    fprintf(gnuplot_pipe, "set ylabel 'y [m]' offset -2,0\n");
    fprintf(gnuplot_pipe, "set zlabel 'z [m]'\n");
    fprintf(gnuplot_pipe,
            "set title 'Simulated orbits up to time %.2f s' offset 0,-7.5\n",
            input_total_sim_time);
    // fprintf(gnuplot_pipe,"set view 70,1,1,1\n");
    fprintf(gnuplot_pipe, "set view equal xyz\n");
    fprintf(gnuplot_pipe, "set xtics offset 0,-1\n");
    fprintf(gnuplot_pipe, "set ytics offset -1,0\n");
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
    for (size_t satellite_index = 0;
         satellite_index < input_satellite_vector.size(); satellite_index++) {
      Satellite current_satellite = input_satellite_vector.at(satellite_index);
      std::array<double, 3> initial_position =
          current_satellite.get_ECI_position();
      fprintf(gnuplot_pipe, "%.17g %.17g %.17g\n", initial_position.at(0),
              initial_position.at(1), initial_position.at(2));

      std::array<double, 3> evolved_position = {};

      double timestep_to_use = input_timestep;
      double current_satellite_time =
          current_satellite.get_instantaneous_time();
      while (current_satellite_time < input_total_sim_time) {
        std::pair<double, int> new_timestep_and_error_code =
            current_satellite.evolve_RK45(input_epsilon, timestep_to_use,
                                          perturbation);
        double new_timestep = new_timestep_and_error_code.first;
        int error_code = new_timestep_and_error_code.second;
        if (error_code != 0) {
          std::cout << "Error detected, halting visualization\n";
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
    fprintf(gnuplot_pipe, "pause mouse keypress\n");

    fprintf(gnuplot_pipe, "exit \n");
    pclose(gnuplot_pipe);

  } else {
    std::cout << "gnuplot not found";
  }

  return;
}

void sim_and_plot_orbital_elem_gnuplot(
    std::vector<Satellite> input_satellite_vector, double input_timestep,
    double input_total_sim_time, double input_epsilon,
    std::string input_orbital_element_name, bool perturbation) {
  if (input_satellite_vector.size() < 1) {
    std::cout << "No input Satellite objects\n";
    return;
  }

  // first, open "pipe" to gnuplot
  FILE *gnuplot_pipe = popen("gnuplot", "w");
  // if it exists
  if (gnuplot_pipe) {
    fprintf(gnuplot_pipe, "set terminal qt size 600,400 font 'Helvetica,14'\n");
    // formatting
    fprintf(gnuplot_pipe, "set xlabel 'Time [s]'\n");
    if (input_orbital_element_name == "Semimajor Axis") {
      fprintf(gnuplot_pipe, "set ylabel '%s [m]'\n",
              input_orbital_element_name.c_str());
    } else if (input_orbital_element_name == "Eccentricity") {
      fprintf(gnuplot_pipe, "set ylabel '%s'\n",
              input_orbital_element_name.c_str());
    } else {
      fprintf(gnuplot_pipe, "set ylabel '%s [rad]'\n",
              input_orbital_element_name.c_str());
    }
    fprintf(gnuplot_pipe, "set title '%s simulated up to time %.2f s'\n",
            input_orbital_element_name.c_str(), input_total_sim_time);
    fprintf(gnuplot_pipe, "set key right bottom\n");

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
                "pplot '-' using 1:2 with lines lw 1 title '%s' \n",
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
    for (size_t satellite_index = 0;
         satellite_index < input_satellite_vector.size(); satellite_index++) {
      Satellite current_satellite = input_satellite_vector.at(satellite_index);
      double val =
          current_satellite.get_orbital_element(input_orbital_element_name);
      double current_satellite_time =
          current_satellite.get_instantaneous_time();
      fprintf(gnuplot_pipe, "%.17g %.17g\n", current_satellite_time, val);

      double evolved_val = {0};

      double timestep_to_use = input_timestep;
      current_satellite_time = current_satellite.get_instantaneous_time();
      while (current_satellite_time < input_total_sim_time) {
        std::pair<double, int> new_timestep_and_error_code =
            current_satellite.evolve_RK45(input_epsilon, timestep_to_use,
                                          perturbation);
        double new_timestep = new_timestep_and_error_code.first;
        int error_code = new_timestep_and_error_code.second;

        if (error_code != 0) {
          std::cout << "Error detected, halting visualization\n";
          fprintf(gnuplot_pipe, "e\n");
          fprintf(gnuplot_pipe, "exit \n");
          pclose(gnuplot_pipe);
          return;
        }
        timestep_to_use = new_timestep;
        evolved_val =
            current_satellite.get_orbital_element(input_orbital_element_name);
        current_satellite_time = current_satellite.get_instantaneous_time();
        fprintf(gnuplot_pipe, "%.17g %.17g\n", current_satellite_time,
                evolved_val);
      }
      fprintf(gnuplot_pipe, "e\n");
    }
    fprintf(gnuplot_pipe, "pause mouse keypress\n");

    fprintf(gnuplot_pipe, "exit \n");
    pclose(gnuplot_pipe);

  } else {
    std::cout << "gnuplot not found";
  }

  return;
}

Matrix3d z_rot_matrix(double input_angle) {
  Matrix3d z_rotation_matrix;
  z_rotation_matrix << cos(input_angle), -sin(input_angle), 0, sin(input_angle),
      cos(input_angle), 0, 0, 0, 1;

  return z_rotation_matrix;
}

Matrix3d y_rot_matrix(double input_angle) {
  Matrix3d y_rotation_matrix;
  y_rotation_matrix << cos(input_angle), 0, sin(input_angle), 0, 1, 0,
      -sin(input_angle), 0, cos(input_angle);

  return y_rotation_matrix;
}

Matrix3d x_rot_matrix(double input_angle) {
  Matrix3d x_rotation_matrix;
  x_rotation_matrix << 1, 0, 0, 0, cos(input_angle), -sin(input_angle), 0,
      sin(input_angle), cos(input_angle);

  return x_rotation_matrix;
}

std::array<double, 4> bodyframe_quaternion_deriv(
    std::array<double, 4> input_bodyframe_quaternion, double input_w_1,
    double input_w_2, double input_w_3) {
  // Assumes quaternions are input as [q_0, q_1, q_2, q_3] and omega is the
  // spacecraft "body rate", represented in the body frame, w.r.t. a chosen
  // reference frame Ref:
  // https://ntrs.nasa.gov/api/citations/20240009554/downloads/Space%20Attitude%20Development%20Control.pdf

  // Let's do this with Eigen so we can use matrix multiplication easily
  Vector4d input_body_quaternion;
  input_body_quaternion << input_bodyframe_quaternion.at(0),
      input_bodyframe_quaternion.at(1), input_bodyframe_quaternion.at(2),
      input_bodyframe_quaternion.at(3);

  Matrix4d coeff_mat;
  coeff_mat << 0, -input_w_1, -input_w_2, -input_w_3, input_w_1, 0, input_w_3,
      -input_w_2, input_w_2, -input_w_3, 0, input_w_1, input_w_3, input_w_2,
      -input_w_1, 0;

  Vector4d quaternion_deriv_vec = 0.5 * coeff_mat * input_body_quaternion;
  std::array<double, 4> quaternion_deriv_array;
  for (size_t ind = 0; ind < quaternion_deriv_array.size(); ind++) {
    quaternion_deriv_array.at(ind) = quaternion_deriv_vec(ind);
  }
  return quaternion_deriv_array;
}

void sim_and_plot_attitude_evolution_gnuplot(
    std::vector<Satellite> input_satellite_vector, double input_timestep,
    double input_total_sim_time, double input_epsilon,
    std::string input_plotted_val_name, bool perturbation) {
  if (input_satellite_vector.size() < 1) {
    std::cout << "No input Satellite objects\n";
    return;
  }

  // first, open "pipe" to gnuplot
  FILE *gnuplot_pipe = popen("gnuplot", "w");
  // if it exists
  if (gnuplot_pipe) {
    fprintf(gnuplot_pipe, "set terminal qt size 600,400 font 'Helvetica,14'\n");
    // formatting
    fprintf(gnuplot_pipe, "set xlabel 'Time [s]'\n");
    if ((input_plotted_val_name == "omega_x") ||
        (input_plotted_val_name == "omega_y") ||
        (input_plotted_val_name == "omega_z")) {
      fprintf(gnuplot_pipe, "set ylabel '%s [rad/s]'\n",
              input_plotted_val_name.c_str());
    } else if ((input_plotted_val_name == "Roll") ||
               (input_plotted_val_name == "Pitch") ||
               (input_plotted_val_name == "Yaw")) {
      fprintf(gnuplot_pipe, "set ylabel '%s [rad]'\n",
              input_plotted_val_name.c_str());
    } else {
      fprintf(gnuplot_pipe, "set ylabel '%s'\n",
              input_plotted_val_name.c_str());
    }
    fprintf(gnuplot_pipe, "set title '%s simulated up to time %.2f s'\n",
            input_plotted_val_name.c_str(), input_total_sim_time);
    fprintf(gnuplot_pipe, "set key right bottom\n");

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
                "pplot '-' using 1:2 with lines lw 1 title '%s' \n",
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
    for (size_t satellite_index = 0;
         satellite_index < input_satellite_vector.size(); satellite_index++) {
      Satellite current_satellite = input_satellite_vector.at(satellite_index);
      double val = current_satellite.get_attitude_val(input_plotted_val_name);
      double current_satellite_time =
          current_satellite.get_instantaneous_time();
      fprintf(gnuplot_pipe, "%.17g %.17g\n", current_satellite_time, val);

      double evolved_val = {0};

      double timestep_to_use = input_timestep;
      current_satellite_time = current_satellite.get_instantaneous_time();
      while (current_satellite_time < input_total_sim_time) {
        std::pair<double, int> new_timestep_and_error_code =
            current_satellite.evolve_RK45(input_epsilon, timestep_to_use,
                                          perturbation);
        double new_timestep = new_timestep_and_error_code.first;
        int error_code = new_timestep_and_error_code.second;

        if (error_code != 0) {
          std::cout << "Error detected, halting visualization\n";
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
    fprintf(gnuplot_pipe, "pause mouse keypress\n");

    fprintf(gnuplot_pipe, "exit \n");
    pclose(gnuplot_pipe);

  } else {
    std::cout << "gnuplot not found";
  }

  return;
}

Matrix3d rollyawpitch_bodyframe_to_LVLH_matrix(double input_roll,
                                               double input_pitch,
                                               double input_yaw) {
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

std::array<double, 4> rollyawpitch_angles_to_quaternion(double input_roll,
                                                        double input_pitch,
                                                        double input_yaw) {
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
    std::array<double, 4> input_bodyframe_quaternion_relative_to_LVLH) {
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

std::array<double, 3> calculate_spacecraft_bodyframe_angular_acceleration(
    Matrix3d J_matrix,
    std::vector<BodyframeTorqueProfile> input_bodyframe_torque_profile_list,
    Vector3d input_omega_I, double input_orbital_angular_acceleration,
    Vector3d input_omega_bodyframe_wrt_LVLH_in_body_frame,
    Matrix3d input_LVLH_to_bodyframe_transformation_matrix,
    Vector3d input_omega_LVLH_wrt_inertial_in_LVLH,
    double input_evaluation_time) {
  // Using approach from
  // https://ntrs.nasa.gov/api/citations/20240009554/downloads/Space%20Attitude%20Development%20Control.pdf
  // , Ch. 4 especially Objective: calculate omega_dot in spacecraft body frame
  // w.r.t. the LVLH frame, expressed in the spacecraft body frame Disturbance
  // torques, if eventually added in, should just be more torque profiles in the
  // satellite object's list of torque profiles So the input_torques vector
  // includes both disturbance and control torques
  Vector3d bodyframe_torque_vec = {0, 0, 0};

  for (BodyframeTorqueProfile bodyframe_torque_profile :
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

Vector4d quaternion_kinematics_equation(
    Vector4d quaternion_of_bodyframe_relative_to_ref_frame,
    Vector3d angular_velocity_vec_wrt_ref_frame_in_body_frame) {
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

std::array<double, 7> RK45_satellite_body_angular_deriv_function(
    std::array<double, 7> combined_bodyframe_angular_array, Matrix3d J_matrix,
    std::vector<BodyframeTorqueProfile> input_bodyframe_torque_profile_list,
    Vector3d input_omega_I, double input_orbital_angular_acceleration,
    Matrix3d input_LVLH_to_bodyframe_transformation_matrix,
    Vector3d input_omega_LVLH_wrt_inertial_in_LVLH,
    double input_evaluation_time) {
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

Vector3d calculate_omega_I(
    Vector3d input_bodyframe_ang_vel_vector_wrt_lvlh,
    Matrix3d input_LVLH_to_bodyframe_transformation_matrix,
    double input_orbital_rate) {
  // Eq. 4.17 in
  // https://ntrs.nasa.gov/api/citations/20240009554/downloads/Space%20Attitude%20Development%20Control.pdf
  Vector3d omega_LVLH_wrt_inertial_in_LVLH = {0, -input_orbital_rate, 0};
  Vector3d omega_I = input_bodyframe_ang_vel_vector_wrt_lvlh +
                     input_LVLH_to_bodyframe_transformation_matrix *
                         omega_LVLH_wrt_inertial_in_LVLH;
  return omega_I;
}

Matrix3d construct_J_matrix(double input_Jxx, double input_Jyy,
                            double input_Jzz) {
  Matrix3d J_matrix = Matrix3d::Zero();
  J_matrix(0, 0) = input_Jxx;
  J_matrix(1, 1) = input_Jyy;
  J_matrix(2, 2) = input_Jzz;
  return J_matrix;
}

std::array<double, 13>
RK45_combined_orbit_position_velocity_attitude_deriv_function(
    std::array<double, 13> combined_position_velocity_bodyframe_angular_array,
    Matrix3d J_matrix,
    std::vector<BodyframeTorqueProfile> input_bodyframe_torque_profile_list,
    Vector3d input_omega_I, double input_orbital_angular_acceleration,
    Matrix3d input_LVLH_to_bodyframe_transformation_matrix,
    Vector3d input_omega_LVLH_wrt_inertial_in_LVLH,
    double input_spacecraft_mass,
    std::vector<ThrustProfileLVLH> input_list_of_thrust_profiles_LVLH,
    double input_evaluation_time, double input_inclination,
    double input_arg_of_periapsis, double input_true_anomaly,
    bool perturbation) {
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
          perturbation);

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
    std::array<double, 4> input_quaternion) {
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
    std::array<double, 3> input_LVLH_frame_array, double input_roll,
    double input_yaw, double input_pitch) {
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