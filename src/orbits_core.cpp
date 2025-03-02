#include <iostream>
#include "Satellite.h"
#include "utils.h"

int main()
{
    Satellite test_sat_1("../input.json");
    Satellite test_sat_2("../input_2.json");

    std::vector<Satellite> satellite_vector={test_sat_1,test_sat_2};

    double timestep=2;
    double total_sim_time=10000;
    sim_and_draw_orbit_gnuplot(satellite_vector,timestep,total_sim_time);

    return 0;
}