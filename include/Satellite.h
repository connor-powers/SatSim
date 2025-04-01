#ifndef SATELLITE_HEADER
#define SATELLITE_HEADER

#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <stdexcept>

using json = nlohmann::json;

// Define constants
const double G =
    6.674 *
    pow(10, -11);  // https://en.wikipedia.org/wiki/Gravitational_constant
const double mass_Earth =
    5.9722 * pow(10, 24);  // https://en.wikipedia.org/wiki/Earth_mass
const double radius_Earth =
    6378137;  // https://en.wikipedia.org/wiki/Earth_radius

class ThrustProfileLVLH {
  // Note: for now, thrust forces are assumed to act through center of mass of
  // satellite.
 public:
  double t_start_ = {0};
  double t_end_ = {0};
  std::array<double, 3> LVLH_force_vec_ = {0, 0, 0};
  ThrustProfileLVLH(const double t_start, const double t_end,
                    const std::array<double, 3> LVLH_force_vec) {
    t_start_ = t_start;
    t_end_ = t_end;
    LVLH_force_vec_ = LVLH_force_vec;
  }
  ThrustProfileLVLH(
      const double t_start, const double t_end,
      const std::array<double, 3> LVLH_normalized_force_direction_vec,
      const double input_force_magnitude) {
    t_start_ = t_start;
    t_end_ = t_end;
    for (size_t ind = 0; ind < 3; ind++) {
      LVLH_force_vec_.at(ind) =
          input_force_magnitude * LVLH_normalized_force_direction_vec.at(ind);
    }
  }

  bool operator == (const ThrustProfileLVLH& input_profile){
    return ((t_start_ == input_profile.t_start_) && (t_end_ == input_profile.t_end_) && (std::equal(LVLH_force_vec_.begin(), LVLH_force_vec_.end(),input_profile.LVLH_force_vec_.begin())));
  }
};

class BodyframeTorqueProfile {
 public:
  double t_start_ = {0};
  double t_end_ = {0};
  std::array<double, 3> bodyframe_torque_list = {0, 0, 0};
  BodyframeTorqueProfile(const double t_start, const double t_end,
                         const std::array<double, 3> bodyframe_torque_vec) {
    t_start_ = t_start;
    t_end_ = t_end;
    bodyframe_torque_list = bodyframe_torque_vec;
  }
  BodyframeTorqueProfile(
      const double t_start, const double t_end,
      const std::array<double, 3> bodyframe_normalized_torque_axis_vec,
      const double input_torque_magnitude) {
    t_start_ = t_start;
    t_end_ = t_end;
    for (size_t ind = 0; ind < 3; ind++) {
      bodyframe_torque_list.at(ind) =
          input_torque_magnitude * bodyframe_normalized_torque_axis_vec.at(ind);
    }
  }
  bool operator == (const BodyframeTorqueProfile& input_profile){
    return ((t_start_ == input_profile.t_start_) && (t_end_ == input_profile.t_end_) && (std::equal(bodyframe_torque_list.begin(), bodyframe_torque_list.end(),input_profile.bodyframe_torque_list.begin())));
  }
};

class Satellite {
 private:
  double inclination_ = {0};
  double raan_ = {
      0};  // Assuming RAAN can be used interchangeably with longitude of
           // ascending node for the Earth-orbiting satellites simulated here
  double arg_of_periapsis_ = {0};
  double eccentricity_ = {0};
  double a_ = {0};
  double true_anomaly_ = {0};
  double orbital_period_ = {0};
  double m_ = {1};  // default value to prevent infinities in acceleration
                    // calculations from a=F/m
  // double I_={1}; //moment of inertia, taken to be same for all 3 principal
  // axes, set to default value for same reasons as mass
  double t_ = {0};

  double orbital_rate_ = {0};
  double orbital_angular_acceleration_ = {0};  // Time derivative of orbital
                                               // rate

  // Now body-frame attributes
  // Assuming diagonal J matrix
  double J_11_ = {1};
  double J_22_ = {1};
  double J_33_ = {1};
  // The following angles are angles of the satellite body frame with respect to
  // the LVLH frame, represented in the body frame
  double pitch_angle_ = {0};
  double roll_angle_ = {0};
  double yaw_angle_ = {0};

  // body-frame angular velocities relative to the LVLH frame, represented in
  // the body frame
  std::array<double, 3> body_angular_velocity_vec_wrt_LVLH_in_body_frame_ = {
      0, 0, 0};

  // quaternion representing attitude of satellite body frame with respect to
  // the LVLH frame
  std::array<double, 4> quaternion_satellite_bodyframe_wrt_LVLH_;

  std::string name_ = "";

  std::array<double, 3> perifocal_position_ = {0, 0, 0};
  std::array<double, 3> perifocal_velocity_ = {0, 0, 0};

  std::array<double, 3> ECI_position_ = {0, 0, 0};
  std::array<double, 3> ECI_velocity_ = {0, 0, 0};

  std::vector<ThrustProfileLVLH> thrust_profile_list_ = {};
  std::vector<BodyframeTorqueProfile> bodyframe_torque_profile_list_ = {};

  std::vector<std::array<double, 3>> list_of_LVLH_forces_at_this_time_ = {};
  std::vector<std::array<double, 3>> list_of_ECI_forces_at_this_time_ = {};
  std::vector<std::array<double, 3>> list_of_body_frame_torques_at_this_time_ =
      {};
  // std::vector<std::array<double,3>>
  // list_of_body_frame_torques_at_this_time_={};

  std::pair<double, double> calculate_eccentric_anomaly(
      const double input_eccentricity, const double input_true_anomaly,
      const double input_semimajor_axis);
  void initialize_body_angular_velocity_vec_wrt_LVLH_in_body_frame();

 public:
  std::string plotting_color_ = "";
  Satellite(const std::string input_file_name) {
    // baselining JSON input file format specifying initial orbital parameters
    // of satellite semimajor axis is read in units of km angles are read in
    // units of degrees, then internally translated to radians
    std::ifstream input_filestream(input_file_name);
    json input_data = json::parse(input_filestream);

    inclination_ = input_data.at("Inclination");
    // convert to radians
    inclination_ *= (M_PI / 180.0);
    if (inclination_ == 0) {
      throw std::invalid_argument(
          "Zero inclination orbits are not currently supported");
    }

    raan_ = input_data.at("RAAN");
    // convert to radians
    raan_ *= (M_PI / 180.0);

    arg_of_periapsis_ = input_data.at("Argument of Periapsis");
    // convert to radians
    arg_of_periapsis_ *= (M_PI / 180.0);

    eccentricity_ = input_data.at("Eccentricity");
    // If circular orbit, arg of periapsis is undefined, using convention of
    // setting it to 0 in this case
    if (eccentricity_ == 0) {
      arg_of_periapsis_ = 0;
    }

    a_ = input_data.at("Semimajor Axis");
    a_ *= 1000.0;  // converting from km to m

    true_anomaly_ = input_data.at("True Anomaly");
    // convert to radians
    true_anomaly_ *= (M_PI / 180.0);

    // making initial pitch angle an optional parameter
    if (input_data.find("Initial Pitch Angle") != input_data.end()) {
      pitch_angle_ = input_data.at("Initial Pitch Angle");
      // convert to radians
      pitch_angle_ *= (M_PI / 180.0);
    }
    // making initial roll angle an optional parameter
    if (input_data.find("Initial Roll Angle") != input_data.end()) {
      roll_angle_ = input_data.at("Initial Roll Angle");
      // convert to radians
      roll_angle_ *= (M_PI / 180.0);
    }
    // making initial yaw angle an optional parameter
    if (input_data.find("Initial Yaw Angle") != input_data.end()) {
      yaw_angle_ = input_data.at("Initial Yaw Angle");
      // convert to radians
      yaw_angle_ *= (M_PI / 180.0);
    }

    initialize_and_normalize_body_quaternion(roll_angle_, pitch_angle_,
                                             yaw_angle_);

    m_ = input_data.at("Mass");
    name_ = input_data.at("Name");

    // making plotting color an optional parameter
    if (input_data.find("Plotting Color") != input_data.end()) {
      plotting_color_ = input_data.at("Plotting Color");
    }

    t_ = 0;  // for now, assuming satellites are initialized at time t=0;

    orbital_period_ = calculate_orbital_period();

    // updated workflow
    perifocal_position_ = calculate_perifocal_position();
    perifocal_velocity_ = calculate_perifocal_velocity();

    ECI_position_ = convert_perifocal_to_ECI(perifocal_position_);
    ECI_velocity_ = convert_perifocal_to_ECI(perifocal_velocity_);

    orbital_rate_ = calculate_instantaneous_orbit_rate();
    initialize_body_angular_velocity_vec_wrt_LVLH_in_body_frame();
    // making initial omega_x an optional parameter
    if (input_data.find("Initial omega_x") != input_data.end()) {
      double initial_omega_x_wrt_LVLH_in_body_frame =
          input_data.at("Initial omega_x");
      body_angular_velocity_vec_wrt_LVLH_in_body_frame_.at(0) +=
          initial_omega_x_wrt_LVLH_in_body_frame;
    }
    // making initial omega_y an optional parameter
    if (input_data.find("Initial omega_y") != input_data.end()) {
      double initial_omega_y_wrt_LVLH_in_body_frame =
          input_data.at("Initial omega_y");
      body_angular_velocity_vec_wrt_LVLH_in_body_frame_.at(1) +=
          initial_omega_y_wrt_LVLH_in_body_frame;
    }

    // making initial omega_z an optional parameter
    if (input_data.find("Initial omega_z") != input_data.end()) {
      double initial_omega_z_wrt_LVLH_in_body_frame =
          input_data.at("Initial omega_z");
      body_angular_velocity_vec_wrt_LVLH_in_body_frame_.at(2) +=
          initial_omega_z_wrt_LVLH_in_body_frame;
    }

    orbital_angular_acceleration_ =
        calculate_instantaneous_orbit_angular_acceleration();
  }

  std::array<double, 3> get_ECI_position() { return ECI_position_; }
  std::array<double, 3> get_ECI_velocity() { return ECI_velocity_; }
  double get_speed() {
    // shouldn't matter which frame I use, might as well use perifocal coords
    // since it's fewer operations (no W-direction component so can omit that
    // term, whereas there's x,y,z components in ECI)
    return sqrt(pow(perifocal_velocity_.at(0), 2) +
                pow(perifocal_velocity_.at(1), 2));
  }
  double get_speed_ECI() {
    return sqrt(pow(ECI_velocity_.at(0), 2) +
                pow(ECI_velocity_.at(1), 2) +
                pow(ECI_velocity_.at(2),2));
  }
  double get_radius() {
    // shouldn't matter which frame I use, might as well use perifocal coords
    // since it's fewer operations (no W-direction component so can omit that
    // term, whereas there's x,y,z components in ECI)
    return sqrt(pow(perifocal_position_.at(0), 2) +
                pow(perifocal_position_.at(1), 2));
  }
  double get_radius_ECI() {
    return sqrt(pow(ECI_position_.at(0), 2) +
                pow(ECI_position_.at(1), 2) +
                pow(ECI_position_.at(2),2));
  }
  double get_total_energy() {
    double orbital_radius = get_radius();
    double gravitational_potential_energy =
        -G * mass_Earth * m_ / orbital_radius;

    double orbital_speed = get_speed();
    double kinetic_energy = (1.0 / 2.0) * m_ * (orbital_speed * orbital_speed);

    return (gravitational_potential_energy + kinetic_energy);
  }

  double get_instantaneous_time() { return t_; }
  std::string get_name() { return name_; }
  void evolve_RK4(const double input_timestep);

  std::array<double, 3> body_frame_to_ECI(
      const std::array<double, 3> input_vector);

  std::array<double, 3> ECI_to_body_frame(
      const std::array<double, 3> input_vector);

  std::array<double, 3> calculate_perifocal_position();

  std::array<double, 3> calculate_perifocal_velocity();

  std::array<double, 3> convert_perifocal_to_ECI(
      const std::array<double, 3> input_perifocal_vec);
  std::array<double, 3> convert_ECI_to_perifocal(
      const std::array<double, 3> input_ECI_vec);

  // std::array<double,3> convert_LVLH_to_ECI(std::array<double,3>
  // input_LVLH_vec);

  void add_LVLH_thrust_profile(
      const std::array<double, 3> input_LVLH_normalized_thrust_direction,
      const double input_LVLH_thrust_magnitude,
      const double input_thrust_start_time, const double input_thrust_end_time);
  void add_LVLH_thrust_profile(
      const std::array<double, 3> input_LVLH_thrust_vector,
      const double input_thrust_start_time, const double input_thrust_end_time);

  void add_bodyframe_torque_profile(
      const std::array<double, 3> input_bodyframe_direction_unit_vec,
      const double input_bodyframe_torque_magnitude,
      const double input_torque_start_time, const double input_torque_end_time);
  void add_bodyframe_torque_profile(
      const std::array<double, 3> input_bodyframe_torque_vector,
      const double input_torque_start_time, const double input_torque_end_time);

  int update_orbital_elements_from_position_and_velocity();
  std::array<double, 6> get_orbital_elements();

  std::pair<double, int> evolve_RK45(const double input_epsilon,
                                     const double input_initial_timestep,
                                     const bool perturbation = true);

  double get_orbital_element(const std::string orbital_element_name);
  double calculate_instantaneous_orbit_rate();
  double calculate_instantaneous_orbit_angular_acceleration();
  void initialize_and_normalize_body_quaternion(const double roll_angle,
                                                const double pitch_angle,
                                                const double yaw_angle);
  double get_attitude_val(const std::string input_attitude_val_name);
  double calculate_orbital_period();

};

#endif