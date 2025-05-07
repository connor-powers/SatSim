#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <iostream>

#include "PhasedArrayGroundStation.h"
#include "utils.h"

using Eigen::Vector3d;

TEST(PhasedArrayGroundStationTests, InitializationTest) {
  double gs_latitude = 2;
  double gs_longitude = 5;
  double gs_altitude = 10;
  double max_beam_angle_from_normal = 65;  // deg
  int num_beams = 2;
  PhasedArrayGroundStation example_ground_station(
      gs_latitude, gs_longitude, gs_altitude, max_beam_angle_from_normal,
      num_beams);
}

TEST(PhasedArrayGroundStationTests, ECI_test_0) {
  double gs_latitude = 2;
  double gs_longitude = 5;
  double gs_altitude = 10;
  double max_beam_angle_from_normal = 65;  // deg
  int num_beams = 2;
  PhasedArrayGroundStation example_ground_station(
      gs_latitude, gs_longitude, gs_altitude, max_beam_angle_from_normal,
      num_beams);
  Vector3d sample_ECI_position = example_ground_station.get_ECI_position(0);
}

TEST(PhasedArrayGroundStationTests, Evolved_ECI_test) {
  double gs_latitude = 2;
  double gs_longitude = 5;
  double gs_altitude = 10;
  double max_beam_angle_from_normal = 65;  // deg
  int num_beams = 2;
  PhasedArrayGroundStation example_ground_station(
      gs_latitude, gs_longitude, gs_altitude, max_beam_angle_from_normal,
      num_beams);

  // Let's look at the difference over 10 seconds
  Vector3d initial_ECI_position = example_ground_station.get_ECI_position(0);
  Vector3d evolved_ECI_position = example_ground_station.get_ECI_position(10);

  // Earth rotates to the east, and at J2000 the ECI y axis points East
  // but at J2000 there's a rotation angle of ~280 degrees between ECEF and ECI
  // so ECI x and y components should be increasing
  EXPECT_TRUE(evolved_ECI_position(0) > initial_ECI_position(0))
      << "Diff: " << evolved_ECI_position(0) - initial_ECI_position(0) << "\n";
  EXPECT_TRUE(evolved_ECI_position(1) > initial_ECI_position(1))
      << "Diff: " << evolved_ECI_position(1) - initial_ECI_position(1) << "\n";
  EXPECT_TRUE(evolved_ECI_position(2) == initial_ECI_position(2));
}
