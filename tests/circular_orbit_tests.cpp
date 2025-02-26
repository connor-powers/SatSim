#include <iostream>
#include <gtest/gtest.h>
#include "Satellite.h"
#include "utils.h"

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

TEST(CircularOrbitTests,TotalEnergyTimestep1){

    Satellite test_satellite("../tests/circular_orbit_test_2_input.json");
    double initial_energy=test_satellite.get_total_energy();
    double test_timestep=1; //s
    int num_timesteps=1;
    test_satellite.evolve_RK4(test_timestep,num_timesteps);
    double evolved_energy=test_satellite.get_total_energy();
    EXPECT_TRUE(abs(initial_energy-evolved_energy)<pow(10,-10)) << "Total energy not preserved within tolerance. Difference: " << initial_energy-evolved_energy << "\n";
}

TEST(CircularOrbitTests,EvolvedOrbitalRadius1){

    Satellite test_satellite("../tests/circular_orbit_test_2_input.json");
    double calculated_initial_radius=test_satellite.get_radius();
    double test_timestep=1; //s
    int num_timesteps=1;
    test_satellite.evolve_RK4(test_timestep,num_timesteps);
    double calculated_evolved_radius=test_satellite.get_radius();

    EXPECT_TRUE(abs(calculated_initial_radius-calculated_evolved_radius)<pow(10,-10)) << "Orbital radius not constant within tolerance. Difference: " << calculated_initial_radius-calculated_evolved_radius << "\n";
}

TEST(CircularOrbitTests,EvolvedOrbitalSpeed1){

    Satellite test_satellite("../tests/circular_orbit_test_2_input.json");
    double calculated_initial_speed=test_satellite.get_speed();
    double test_timestep=1; //s
    int num_timesteps=1;
    test_satellite.evolve_RK4(test_timestep,num_timesteps);
    double calculated_evolved_speed=test_satellite.get_speed();

    EXPECT_TRUE(abs(calculated_initial_speed-calculated_evolved_speed)<pow(10,-10)) << "Orbital speed not constant within tolerance. Difference: " << calculated_initial_speed-calculated_evolved_speed << "\n";
}