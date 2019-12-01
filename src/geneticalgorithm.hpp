#ifndef EVOLVING_ROBOT_GENETIC_ALGORITHM_HPP
#define EVOLVING_ROBOT_GENETIC_ALGORITHM_HPP

#include "simulate.hpp"

// program constants
#define NUM_OF_TRIALS 5
#define POP_SIZE 60 // needs to be even number
#define NUM_OF_EVALS 1000
#define DIVERSITY_TXT "diversity.txt"
#define OPENGL_TXT "cubes_ga.txt"
#define LEARNING_TXT "learning_ga.txt"

// global macros
#define PI 3.14159
#define INITIAL_HEIGHT 0.0000001
#define PROB_OF_MUT 0.6
#define PROB_PER_PARAM 0.333333
#define MIN_SWING 0.9
#define MAX_SWING 1.1

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
void loop(int);
Cube initialize_cube();
std::vector<int> randomize_array_of_springs();
void crossover(Cube &, Cube &);
void mutation(Cube &);
void tournament_selection(std::vector<Cube> &, std::vector<Cube> &, std::vector<Cube> &);
void calculate_diversity(std::vector<Cube> &, std::ofstream &);


#endif //EVOLVING_ROBOT_GENETIC_ALGORITHM_HPP
