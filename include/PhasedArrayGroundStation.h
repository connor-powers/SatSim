#ifndef GROUNDSTATION_HEADER
#define GROUNDSTATION_HEADER

#include <iostream>
#include <Eigen/Dense>
#include <stdexcept>
#include "Satellite.h"

using Eigen::Vector3d;

class PhasedArrayGroundStation {
    // Phased array ground station class
    public:
      double latitude = {0};
      double longitude = {0};
      double height = {0};
      Vector3d ECEF_position_;
      double max_beam_angle_from_normal_ = {0}; //deg
      int num_beams_ = {0};
      // To enforce the maximum number of satellites a ground station can communicate with at any one time,
      // going to create a map where each key is the satellite index of a given satellite (i.e., index in a satellite vector being plotted)
      // Each key has a val which is a vector of pairs, where the first element is the start of a time range where the
      // ground station is communicating with that satellite, and the second element is the end time for that range
      // So communications with satellites are characterized by vectors of time ranges

      std::map<size_t,std::vector<std::pair<double,double>>> linked_sats_map_;
  
      PhasedArrayGroundStation(double latitude, double longitude, double height, double max_beam_angle_from_normal, int num_beams = 1) {
        generate_ECEF_position();
        num_beams_ = num_beams;
        max_beam_angle_from_normal_ = max_beam_angle_from_normal;
      }
      void generate_ECEF_position();
      Vector3d get_ECI_position(double input_time);
      double distance_to_satellite(Satellite input_satellite);
      double angle_to_satellite_from_normal(Satellite input_satellite);
      int num_sats_connected_at_this_time(double input_time);
      void update_linked_sats_map(const size_t satellite_index, const double connection_time, const double previous_time);
  };

#endif