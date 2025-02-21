#include <iostream>
#include <gtest/gtest.h>
#include "Satellite.h"

//Testing calculated orbital speed based on input orbital parameters at e=0 against known formula for circular orbital speed
//https://en.wikipedia.org/wiki/Circular_orbit#Velocity

TEST(CircularOrbitTests,OrbitalSpeed1){

    Satellite test_satellite("../tests/circular_orbit_test_1_input.json");
    double calculated_radius=test_satellite.get_radius();
    double calculated_speed=test_satellite.get_speed();
    double expected_circular_orbital_speed=sqrt(G*mass_Earth/calculated_radius);
    EXPECT_TRUE(abs(calculated_speed-expected_circular_orbital_speed)<pow(10,-10)) << "Calculated orbital speed did not match expected value within tolerance. Difference: " << calculated_speed-expected_circular_orbital_speed << "\n";
}

TEST(CircularOrbitTests,OrbitalSpeed2){

    Satellite test_satellite("../tests/circular_orbit_test_2_input.json");
    double calculated_radius=test_satellite.get_radius();
    double calculated_speed=test_satellite.get_speed();
    double expected_circular_orbital_speed=sqrt(G*mass_Earth/calculated_radius);
    EXPECT_TRUE(abs(calculated_speed-expected_circular_orbital_speed)<pow(10,-10)) << "Calculated orbital speed did not match expected value within tolerance. Difference: " << calculated_speed-expected_circular_orbital_speed << "\n";
}