#ifndef EVOLVING_ROBOT_GENETIC_ALGORITHM_HPP
#define EVOLVING_ROBOT_GENETIC_ALGORITHM_HPP

#include "simulate.hpp"

// program constants
#define NUM_OF_TRIALS 1
#define POP_SIZE 30 // needs to be even number
#define NUM_OF_EVALS 1000

// global macros
#define PI 3.14159
#define OMEGA (4*PI)
#define INITIAL_HEIGHT 0.000000001
#define PROB_OF_MUT 0.3
#define PROB_PER_PARAM 0.333333
#define MIN_SWING 0.8
#define MAX_SWING 1.2

// spring macros
#define MIN_K_SPRING 500
#define MAX_K_SPRING 2000  // 1000 is a good value, so do times .5 and times 2 for min and max
#define MIN_A L0_SIDE //0
#define MAX_A L0_SIDE //(2*L0_SIDE)
#define MIN_B 0 // a + b*sin(wt+c)
#define MAX_B (L0_SIDE/2)
#define MIN_C 0
#define MAX_C (PI/2)

// where all the GA stuff happens
void loop();
Cube initialize_cube();
std::vector<int> randomize_array_of_springs();
void crossover(Cube &, Cube &);
void mutation(Cube &);
void tournament_selection(std::vector<Cube> &, std::vector<Cube> &, std::vector<Cube> &);


#endif //EVOLVING_ROBOT_GENETIC_ALGORITHM_HPP
