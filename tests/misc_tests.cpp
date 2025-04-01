#include <gtest/gtest.h>

#include <iostream>

#include "Satellite.h"
#include "utils.h"

TEST(MiscTests, ThrustProfileInitializationTest1) {
  double t_start = 1.0;
  double t_end = 9.0;
  std::array<double,3> thrust_vec = {1.0,101.2,-0.4};
  double magnitude = sqrt(pow(thrust_vec.at(0),2) + pow(thrust_vec.at(1),2) + pow(thrust_vec.at(2),2));
  std::array<double,3> thrust_vec_direction = {0.0,0.0,0.0};
  for (size_t ind=0;ind<thrust_vec.size();ind++) {
    thrust_vec_direction.at(ind) = thrust_vec.at(ind)/magnitude;
  }

  ThrustProfileLVLH thrust_profile_1(t_start,t_end,thrust_vec);
  ThrustProfileLVLH thrust_profile_2(t_start,t_end,thrust_vec_direction,magnitude);

  EXPECT_TRUE(thrust_profile_1 == thrust_profile_2)
      << "Thrust profiles initialized differently didn't agree.\n";
}

TEST(MiscTests, TorqueProfileInitializationTest1) {
  double t_start = 1.0;
  double t_end = 9.0;
  std::array<double,3> torque_vec = {1.0,101.2,-0.4};
  double magnitude = sqrt(pow(torque_vec.at(0),2) + pow(torque_vec.at(1),2) + pow(torque_vec.at(2),2));
  std::array<double,3> torque_vec_direction = {0.0,0.0,0.0};
  for (size_t ind=0;ind<torque_vec.size();ind++) {
    torque_vec_direction.at(ind) = torque_vec.at(ind)/magnitude;
  }

  BodyframeTorqueProfile torque_profile_1(t_start,t_end,torque_vec);
  BodyframeTorqueProfile torque_profile_2(t_start,t_end,torque_vec_direction,magnitude);

  EXPECT_TRUE(torque_profile_1 == torque_profile_2)
      << "Torque profiles initialized differently didn't agree.\n";
}