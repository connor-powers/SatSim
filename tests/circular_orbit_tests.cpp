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
    test_satellite.evolve_RK4(test_timestep);
    double evolved_energy=test_satellite.get_total_energy();
    EXPECT_TRUE(abs(initial_energy-evolved_energy)<pow(10,-10)) << "Total energy not preserved within tolerance. Difference: " << initial_energy-evolved_energy << "\n";
}

TEST(CircularOrbitTests,EvolvedOrbitalRadius1){

    Satellite test_satellite("../tests/circular_orbit_test_2_input.json");
    double calculated_initial_radius=test_satellite.get_radius();
    double test_timestep=1; //s
    test_satellite.evolve_RK4(test_timestep);
    double calculated_evolved_radius=test_satellite.get_radius();

    EXPECT_TRUE(abs(calculated_initial_radius-calculated_evolved_radius)<pow(10,-10)) << "Orbital radius not constant within tolerance. Difference: " << calculated_initial_radius-calculated_evolved_radius << "\n";
}

TEST(CircularOrbitTests,EvolvedOrbitalSpeed1){

    Satellite test_satellite("../tests/circular_orbit_test_2_input.json");
    double calculated_initial_speed=test_satellite.get_speed();
    double test_timestep=1; //s
    test_satellite.evolve_RK4(test_timestep);
    double calculated_evolved_speed=test_satellite.get_speed();

    EXPECT_TRUE(abs(calculated_initial_speed-calculated_evolved_speed)<pow(10,-10)) << "Orbital speed not constant within tolerance. Difference: " << calculated_initial_speed-calculated_evolved_speed << "\n";
}



TEST(CircularOrbitTests,BasicOrbitalElementsTest){

    Satellite test_satellite("../tests/circular_orbit_test_2_input.json");
    std::array<double,6> initial_orbit_elements=test_satellite.get_orbital_elements();

    test_satellite.update_orbital_elements_from_position_and_velocity();
    std::array<double,6> recalculated_orbit_elements=test_satellite.get_orbital_elements();

    std::array<std::string,6> orbital_element_name_array;
    orbital_element_name_array.at(0)="Semimajor Axis";
    orbital_element_name_array.at(1)="Eccentricity";
    orbital_element_name_array.at(2)="Inclination";
    orbital_element_name_array.at(3)="RAAN";
    orbital_element_name_array.at(4)="Argument of Periapsis";
    orbital_element_name_array.at(5)="True Anomaly";


    for (size_t orbital_elem_index=0; orbital_elem_index<6;orbital_elem_index++){
        if (orbital_elem_index==0){
            //Setting a different tolerance for semimajor axis than the other orbital parameters since its scale is so different
            EXPECT_TRUE(abs(initial_orbit_elements.at(orbital_elem_index)-recalculated_orbit_elements.at(orbital_elem_index))<pow(10,-7)) << orbital_element_name_array.at(orbital_elem_index) << " was not constant within tolerance. Diff:" <<initial_orbit_elements.at(orbital_elem_index)-recalculated_orbit_elements.at(orbital_elem_index) << "\n";
        }
        else {
            EXPECT_TRUE(abs(initial_orbit_elements.at(orbital_elem_index)-recalculated_orbit_elements.at(orbital_elem_index))<pow(10,-14)) << orbital_element_name_array.at(orbital_elem_index) << " was not constant within tolerance. Diff:" <<initial_orbit_elements.at(orbital_elem_index)-recalculated_orbit_elements.at(orbital_elem_index) << "\n";
        }
    }
}

TEST(CircularOrbitTests,ConstantEvolvedOrbitalElementsTest){

    Satellite test_satellite("../tests/circular_orbit_test_2_input.json");
    std::array<double,6> initial_orbit_elements=test_satellite.get_orbital_elements();

    double test_timestep=1; //s
    test_satellite.evolve_RK4(test_timestep);
    std::array<double,6> evolved_orbit_elements=test_satellite.get_orbital_elements();

    std::array<std::string,6> orbital_element_name_array;
    orbital_element_name_array.at(0)="Semimajor Axis";
    orbital_element_name_array.at(1)="Eccentricity";
    orbital_element_name_array.at(2)="Inclination";
    orbital_element_name_array.at(3)="RAAN";
    orbital_element_name_array.at(4)="Argument of Periapsis";
    orbital_element_name_array.at(5)="True Anomaly";


    for (size_t orbital_elem_index=0; orbital_elem_index<5;orbital_elem_index++){
        //True anomaly shouldn't be constant over evolution
        if (orbital_elem_index==0){
            //Setting a different tolerance for semimajor axis than the other orbital parameters since its scale is so different
            EXPECT_TRUE(abs(initial_orbit_elements.at(orbital_elem_index)-evolved_orbit_elements.at(orbital_elem_index))<pow(10,-7)) << orbital_element_name_array.at(orbital_elem_index) << " was not constant within tolerance. Diff:" <<initial_orbit_elements.at(orbital_elem_index)-evolved_orbit_elements.at(orbital_elem_index) << "\n";
        }
        else {
            EXPECT_TRUE(abs(initial_orbit_elements.at(orbital_elem_index)-evolved_orbit_elements.at(orbital_elem_index))<pow(10,-14)) << orbital_element_name_array.at(orbital_elem_index) << " was not constant within tolerance. Diff:" <<initial_orbit_elements.at(orbital_elem_index)-evolved_orbit_elements.at(orbital_elem_index) << "\n";
        }
    }
}




TEST(CircularOrbitTests,Thruster_Eccentricity_Change){

    Satellite test_satellite("../tests/circular_orbit_test_2_input.json");
    std::array<double,3> LVLH_thrust_direction={1,0,0};
    double thrust_magnitude=100; //N
    double t_thrust_start=1;
    double t_thrust_end=100;

    test_satellite.add_LVLH_thrust_profile(LVLH_thrust_direction,thrust_magnitude,t_thrust_start,t_thrust_end);
    double test_timestep=1; //s
    for (int timestep_ind=0;timestep_ind<120;timestep_ind++){
        test_satellite.evolve_RK4(test_timestep);
    }
    std::array<double,6> evolved_orbit_elements=test_satellite.get_orbital_elements();
    double resulting_eccentricity=evolved_orbit_elements.at(1);

    EXPECT_TRUE(resulting_eccentricity>0) << "Resulting eccentricity was not greater than 0. Calculated value: " << resulting_eccentricity << "\n";
}