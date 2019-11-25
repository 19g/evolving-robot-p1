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
#define INITIAL_HEIGHT 0.02
#define PROB_OF_MUT 0.1
#define MIN_SWING 0.8
#define MAX_SWING 1.2

// spring macros
#define MIN_K_SPRING 100
#define MAX_K_SPRING 1000000
#define MIN_A 0
#define MAX_A (2*L0_SIDE)
#define MIN_B 0
#define MAX_B 5
#define MIN_C 0
#define MAX_C (2*PI)
#define MIN_D 0
#define MAX_D 1
#define MIN_E 0
#define MAX_E 1


/*
struct Cube {
    std::vector<Mass> mass;
    std::vector<Spring> spring;
    double fitness;
};*/

// where all the GA stuff happens
void loop();
Cube initialize_cube();
std::vector<int> randomize_array_of_springs();
void crossover(Cube &, Cube &);
void mutation(Cube &);
void tournament_selection(std::vector<Cube> &, std::vector<Cube> &, std::vector<Cube> &);


#endif //EVOLVING_ROBOT_GENETIC_ALGORITHM_HPP
