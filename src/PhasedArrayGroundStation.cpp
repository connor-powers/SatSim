#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include "PhasedArrayGroundStation.h"
#include "utils.h"

using Eigen::Vector3d;

void PhasedArrayGroundStation::generate_ECEF_position() {
    ECEF_position_ = convert_lat_long_to_ECEF(latitude, longitude, height);
}

Vector3d PhasedArrayGroundStation::get_ECI_position(double input_time) {
    Vector3d ECI_position_vec = convert_ECEF_to_ECI(ECEF_position_, input_time);
    return ECI_position_vec;
  }

double PhasedArrayGroundStation::distance_to_satellite(Satellite input_satellite) {
    // Computes Euclidean distance to a satellite
    std::array<double,3> satellite_position_ECI = input_satellite.get_ECI_position();
    double current_time = input_satellite.get_instantaneous_time();
    Vector3d groundstation_position_ECI = get_ECI_position(current_time);

    // Convert to Eigen Vector3d
    Vector3d satellite_position_vec;
    for (size_t ind=0; ind<satellite_position_ECI.size();ind++) {
        satellite_position_vec(ind) = satellite_position_ECI.at(ind);
    }
    // Order doesn't matter here
    Vector3d difference_vec = (satellite_position_vec - groundstation_position_ECI);
    return difference_vec.norm();
}

double PhasedArrayGroundStation::angle_to_satellite_from_normal(Satellite input_satellite) {

    // Computes angle between vector pointing from ground station to satellite
    // and the vector normal to the ground station's antenna plane
    // Note: for now, this antenna plane is assumed to be tangential to the Earth
    // at the antenna's location.

    // Ref: https://archive.aoe.vt.edu/lutze/AOE4134/12SatelliteLookAngle.pdf
    // but I don't think you need to use that topocentric-horizontal coordinate system
    // to just get angle between antenna plane normal and station-satellite vector
    double current_time = input_satellite.get_instantaneous_time();
    Vector3d groundstation_position_ECI = get_ECI_position(current_time);

    std::array<double,3> satellite_position_ECI = input_satellite.get_ECI_position();
    // Convert to Eigen Vector3d
    Vector3d satellite_position_vec;
    for (size_t ind=0; ind<satellite_position_ECI.size();ind++) {
        satellite_position_vec(ind) = satellite_position_ECI.at(ind);
    }

    Vector3d station_to_satellite_vec = satellite_position_vec - groundstation_position_ECI;
    // Normalize to avoid needing magnitudes in angle calculation
    // So this is just a unit vector now
    station_to_satellite_vec.normalize();

    // Antenna plane normal is (given current assumptions) in the same
    // direction as the vector from the center of the Earth to the ground station
    Vector3d antenna_plane_normal_direction = groundstation_position_ECI;
    antenna_plane_normal_direction.normalize();

    // See, e.g., https://www.wikihow.com/Find-the-Angle-Between-Two-Vectors
    double angle_radians = acos(station_to_satellite_vec.dot(antenna_plane_normal_direction));
    // Make sure returned angle is in degrees
    return (angle_radians * (180/M_PI));
}

int PhasedArrayGroundStation::num_sats_connected_at_this_time(double input_time) {
    // Goal: check the number of satellites whose time ranges of being linked with this ground station
    // include the current time

    // Method: two-pointer search along each element of linked_sats_list to find 
    // first and last ranges that overlap the input time
    int connected_sats = 0;

    for (auto key_val_pair : linked_sats_map_) {
        size_t satellite_ind = key_val_pair.first;
        std::vector<std::pair<double,double>> satellite_comms_range_list = key_val_pair.second;
        // For a given satellite
        // I just need to identify whether there exists a single range that encompasses the input time
        // I can do this by just finding the last range where the t_start <= input_time and 
        // see if t_end >= input_time
        if (satellite_comms_range_list.size() == 1) {
            if ((satellite_comms_range_list.at(0).first <= input_time) && (satellite_comms_range_list.at(0).second >= input_time)) {
                connected_sats += 1;
                continue;
            }
        }
        int right_ind = satellite_comms_range_list.size() - 1;
        int last_matching_pair_ind;
        while (right_ind >= 0) {
            std::pair<double,double> right_range = satellite_comms_range_list.at(right_ind);
            if (right_range.second < input_time) {
                break;
            }
            if (right_range.first > input_time) {
                right_ind -= 1;
            }
            else {
                if (right_range.second >= input_time) {
                    connected_sats += 1;
                }
                break;
            }
        }
    }
    return connected_sats;
}

void PhasedArrayGroundStation::update_linked_sats_map(const size_t satellite_index, const double connection_time, const double previous_time) {
    // If you're just extending an existing range, update the t_end of existing pair,
    // otherwise create a new pair with the current time as the t_start
    if (linked_sats_map_.count(satellite_index) == 0) {
        std::pair<double,double> connection_time_range = {connection_time,connection_time};
        std::vector<std::pair<double,double>> tmp_vector = {connection_time_range};
        linked_sats_map_.emplace(satellite_index,tmp_vector);
    }
    else {
        
        // Look at most recently added range
        std::vector<std::pair<double,double>> existing_connection_time_ranges = linked_sats_map_[satellite_index];
        std::pair<double,double> most_recent_range = existing_connection_time_ranges.at(existing_connection_time_ranges.size() - 1);
        if (most_recent_range.second == previous_time) {
            // Just extending existing range
            most_recent_range.second = connection_time;
            existing_connection_time_ranges.at(existing_connection_time_ranges.size() - 1) = most_recent_range;
            linked_sats_map_[satellite_index] = existing_connection_time_ranges;
        }
        else {
            // New connection range
            std::vector<std::pair<double,double>> existing_connection_time_ranges = linked_sats_map_[satellite_index];
            std::pair<double,double> new_connection_time_range = {connection_time,connection_time};
            existing_connection_time_ranges.push_back(new_connection_time_range);
            linked_sats_map_[satellite_index] = existing_connection_time_ranges;
        }
    }
}
