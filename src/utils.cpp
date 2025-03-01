#include <iostream>
#include <cmath>
#include "Satellite.h"
//baselining cartesian coordinates in ECI frame

std::array<double,3> calculate_orbital_acceleration(const std::array<double,3> input_r_vec)
{
    //orbital acceleration = -G m_Earth/distance^3 * r_vec (just based on rearranging F=ma with a the acceleration due to gravitational attraction between satellite and Earth https://en.wikipedia.org/wiki/Newton%27s_law_of_universal_gravitation)
    //going to be assuming Earth's position doesn't change for now
    //also assuming Earth is spherical, can loosen this assumption in the future
    std::array<double,3> acceleration_vec=input_r_vec; //shouldn't need to explicitly call copy function because input_r_vec is passed by value not ref
    


    const double distance=sqrt(input_r_vec.at(0)*input_r_vec.at(0) + input_r_vec.at(1)*input_r_vec.at(1) + input_r_vec.at(2)*input_r_vec.at(2));
    const double overall_factor= -G*mass_Earth/pow(distance,3);



    for (size_t ind=0;ind<input_r_vec.size();ind++){
        acceleration_vec.at(ind)*=overall_factor;
    }

    return acceleration_vec;
}

std::array<double,6> RK4_deriv_function_orbit_position_and_velocity(std::array<double,6> input_position_and_velocity){
    std::array<double,6> derivative_of_input_y={};
    std::array<double,3> position_array={};

    for (size_t ind=0;ind<3;ind++){
        derivative_of_input_y.at(ind)=input_position_and_velocity.at(ind+3);
        position_array.at(ind)=input_position_and_velocity.at(ind);
    }
    
    std::array<double,3> calculated_orbital_acceleration=calculate_orbital_acceleration(position_array);

    for (size_t ind=3;ind<6;ind++){
        derivative_of_input_y.at(ind)=calculated_orbital_acceleration.at(ind-3);
    }

    return derivative_of_input_y;
}






void sim_and_draw_orbit_gnuplot(Satellite input_satellite,double input_timestep, double input_total_sim_time){
    std::string satellite_name=input_satellite.get_name();
    //first, open "pipe" to gnuplot
    FILE *gnuplot_pipe = popen("gnuplot -persist", "w");
    //if it exists
    if (gnuplot_pipe){


        //formatting
        fprintf(gnuplot_pipe,"set xlabel 'x'\n");
        fprintf(gnuplot_pipe,"set ylabel 'y'\n");
        fprintf(gnuplot_pipe,"set zlabel 'z'\n");
        fprintf(gnuplot_pipe,"set title '%s orbit up to time %.2f'\n",satellite_name.c_str(),input_total_sim_time);
        fprintf(gnuplot_pipe,"set view 45,45\n");





        //plotting
        //first let's set the stage for plotting the Earth
        fprintf(gnuplot_pipe,"R_Earth=%f\n",radius_Earth);
        fprintf(gnuplot_pipe,"set isosamples 50,50\n");
        fprintf(gnuplot_pipe,"set parametric\n");
        fprintf(gnuplot_pipe,"set urange [-pi/2:pi/2]\n");
        fprintf(gnuplot_pipe,"set vrange [0:2*pi]\n");




        fprintf(gnuplot_pipe,"splot R_Earth*cos(u)*cos(v),R_Earth*cos(u)*sin(v),R_Earth*sin(u) lw 1 lc rgb 'blue' notitle,'-' lw 1 lc rgb 'violet' notitle\n");
        //now the orbit data, inline
        std::array<double,3> initial_position=input_satellite.get_position();
        fprintf(gnuplot_pipe,"%f %f %f\n",initial_position.at(0),initial_position.at(1),initial_position.at(2));



        std::array<double,3> evolved_position={};

        int num_timesteps=std::ceil(input_total_sim_time/input_timestep);

        for (int timestep=0;timestep<num_timesteps;timestep++){
            input_satellite.evolve_RK4(input_timestep);
            evolved_position=input_satellite.get_position();
            fprintf(gnuplot_pipe,"%f %f %f\n",evolved_position.at(0),evolved_position.at(1),evolved_position.at(2));
        }
        fprintf(gnuplot_pipe,"e\n");


    
        fprintf(gnuplot_pipe,"exit \n");
        pclose(gnuplot_pipe);



    }
    else {
        std::cout << "gnuplot not found";
    }

    return;
}