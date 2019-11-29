#ifndef EVOLVING_ROBOT_RANDOM_SEARCH_HPP
#define EVOLVING_ROBOT_RANDOM_SEARCH_HPP

#include "simulate.hpp"

// program constants
#define NUM_OF_TRIALS 1
#define NUM_OF_EVALS 30000
#define OPENGL_TXT "cubes_rand.txt"

// global macros
#define PI 3.14159
#define OMEGA (4*PI)
#define INITIAL_HEIGHT 0.000000001

// spring macros
#define MIN_K_SPRING 500
#define MAX_K_SPRING 2000  // 1000 is a good value, so do times .5 and times 2 for min and max
#define MIN_A L0_SIDE //0
#define MAX_A L0_SIDE //(2*L0_SIDE)
#define MIN_B 0 // a + b*sin(wt+c)
#define MAX_B (L0_SIDE/2)
#define MIN_C 0
#define MAX_C (PI/2)

/*
struct Cube {
    std::vector<Mass> mass;
    std::vector<Spring> spring;
    double fitness;
};*/

// where all the GA stuff happens
void random_search();
Cube initialize_cube();

#endif //EVOLVING_ROBOT_GENETIC_ALGORITHM_HPP
