#include <iostream>
#include <gtest/gtest.h>
#include "Satellite.h"
#include "utils.h"


TEST(EllipticalOrbitTests,EvolvedOrbitalSpeed1){
    //Starting at true anomaly=0 means it's starting at perigee, which is where its orbital speed should be maximum
    Satellite test_satellite("../tests/elliptical_orbit_test_1.json");
    double calculated_initial_speed=test_satellite.get_speed();
    double test_timestep=1; //s
    test_satellite.evolve_RK4(test_timestep);
    double calculated_evolved_speed=test_satellite.get_speed();

    EXPECT_TRUE(calculated_initial_speed>calculated_evolved_speed) << "Perigee speed not larger than calculated evolved speed. Difference: " << calculated_initial_speed-calculated_evolved_speed << "\n";
}

TEST(EllipticalOrbitTests,EvolvedOrbitalSpeed2){
    //Starting at true anomaly=180 means it's starting at apogee, which is where its orbital speed should be minimum
    Satellite test_satellite("../tests/elliptical_orbit_test_2.json");
    double calculated_initial_speed=test_satellite.get_speed();
    double test_timestep=1; //s
    test_satellite.evolve_RK4(test_timestep);
    double calculated_evolved_speed=test_satellite.get_speed();

    EXPECT_TRUE(calculated_initial_speed<calculated_evolved_speed) << "Apogee speed not smaller than calculated evolved speed. Difference: " << calculated_initial_speed-calculated_evolved_speed << "\n";
}