#include <iostream>
#include "Satellite.h"
#include "utils.h"

int main()
{
    Satellite test_sat("../input.json");
    double timestep=2;
    double total_sim_time=10000;
    sim_and_draw_orbit_gnuplot(test_sat,timestep,total_sim_time);

    return 0;
}