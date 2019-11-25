#ifndef PHYSICS_SIM_GENETIC_ALGORITHM_HPP
#define PHYSICS_SIM_GENETIC_ALGORITHM_HPP

using namespace std;

#include "simulate.hpp"

// program constants
#define NUM_OF_TRIALS 10
#define POP_SIZE 30
#define NUM_OF_EVALS 10000

// spring macros
#define MIN_K_SPRING 100
#define MAX_K_SPRING 1000000
#define MIN_A 0
#define MAX_A 2*L0_SIDE
#define MIN_B 0
#define MAX_B 1
#define MIN_C 0
#define MAX_C 1
#define MIN_D 0
#define MAX_D 1
#define MIN_E 0
#define MAX_E 1

// global macros
#define OMEGA 3.14159
#define INITIAL_HEIGHT 0.02
#define PROB_OF_MUT 0.1
#define MIN_SWING 0.8
#define MAX_SWING 1.2

/*
struct Cube {
    vector<Mass> mass;
    vector<Spring> spring;
    double fitness;
};*/

// where all the GA stuff happens
void loop();
Cube initialize_cube();
vector<int> randomize_array_of_springs();
void crossover(Cube &, Cube &);
void mutation(Cube &);
void tournament_selection(vector<Cube> &, vector<Cube> &, vector<Cube> &);


#endif //PHYSICS_SIM_GENETICALGORITHM_HPP
