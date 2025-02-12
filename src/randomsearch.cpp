#include <iostream>
#include <cmath>
#include <thread>
#include <random>
#include <vector>
#include <fstream>
#include "simulate.hpp"
#include "randomsearch.hpp"

using namespace std;

int main() {
    srand(time(nullptr));

    // use threads for performance
    thread t[NUM_OF_TRIALS];
    for (int i = 0; i < NUM_OF_TRIALS; i++) {
        t[i] = thread(random_search, i);
    }
    for (int i = 0; i < NUM_OF_TRIALS; i++) {
        t[i].join();
    }

    return 0;
}

void random_search(int thread_num) {
    // begin timer
    clock_t begin = clock();

    // initialize file
    ofstream learning_file;
    learning_file.open(to_string(thread_num) + LEARNING_TXT);

    // initialize parent population randomly
    Cube parent = initialize_cube();
    simulation_loop(parent, thread_num, false);

    // random search loop
    for (int eval = 0; eval <= NUM_OF_EVALS; eval++) {
        // initialize offspring
        Cube child = initialize_cube();
        simulation_loop(child, thread_num, false);

        if (parent.fitness < child.fitness) {
            parent = child;
        }

        // write to file
        learning_file << parent.distance;
        if (eval != NUM_OF_EVALS - 1) {
            learning_file << ",";
        }

        if (eval % 1 == 0) {
            cout << eval << ": " << parent.distance << "\n";
        }
    }

    // end timer
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    double iters_per_sec = NUM_OF_ITERATIONS / elapsed_secs;
    cout << "iter/sec: " << iters_per_sec << "\n";

    // output to file for opengl
    simulation_loop(parent, thread_num, true);

    // close file
    learning_file.close();
}

Cube initialize_cube() {
    // random variable generation
    random_device rd;
    mt19937 mt(rd());
//    uniform_real_distribution<double> a_val(MIN_A, MAX_A);
    uniform_real_distribution<double> b_val(MIN_B, MAX_B);
    uniform_real_distribution<double> c_val(MIN_C, MAX_C);
    uniform_real_distribution<double> k_spring_val(MIN_K_SPRING, MAX_K_SPRING);

    // temp objects
    vector<Mass> mass;
    vector<Spring> spring;

    // create masses
    Mass temp_mass = {WEIGHT_PER_MASS,
                      {0.0, 0.0, 0.0},
                      {0.0, 0.0, 0.0},
                      {0.0, 0.0, 0.0}};
    temp_mass.p = {0.0, 0.0, INITIAL_HEIGHT}; // 0
    mass.emplace_back(temp_mass);
    temp_mass.p = {L0_SIDE, 0.0, INITIAL_HEIGHT}; // 1
    mass.emplace_back(temp_mass);
    temp_mass.v = {0.0, 0.0, 0.0}; // 0
    temp_mass.p = {L0_SIDE, L0_SIDE, INITIAL_HEIGHT}; // 2
    mass.emplace_back(temp_mass);
    temp_mass.p = {0.0, L0_SIDE, INITIAL_HEIGHT}; // 3
    mass.emplace_back(temp_mass);
    temp_mass.p = {0.0, 0.0, INITIAL_HEIGHT + L0_SIDE}; // 4
    mass.emplace_back(temp_mass);
    temp_mass.p = {L0_SIDE, 0.0, INITIAL_HEIGHT + L0_SIDE}; // 5
    mass.emplace_back(temp_mass);
    temp_mass.p = {L0_SIDE, L0_SIDE, INITIAL_HEIGHT + L0_SIDE}; // 6
    mass.emplace_back(temp_mass);
    temp_mass.p = {0.0, L0_SIDE, INITIAL_HEIGHT + L0_SIDE}; // 7
    mass.emplace_back(temp_mass);

    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[0].p, mass[1].p), 0, 1, dist(mass[0].p, mass[1].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[0].p, mass[2].p), 0, 2, dist(mass[0].p, mass[2].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[0].p, mass[3].p), 0, 3, dist(mass[0].p, mass[3].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[0].p, mass[4].p), 0, 4, dist(mass[0].p, mass[4].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[0].p, mass[5].p), 0, 5, dist(mass[0].p, mass[5].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[0].p, mass[6].p), 0, 6, dist(mass[0].p, mass[6].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[0].p, mass[7].p), 0, 7, dist(mass[0].p, mass[7].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[1].p, mass[2].p), 1, 2, dist(mass[1].p, mass[2].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[1].p, mass[3].p), 1, 3, dist(mass[1].p, mass[3].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[1].p, mass[4].p), 1, 4, dist(mass[1].p, mass[4].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[1].p, mass[5].p), 1, 5, dist(mass[1].p, mass[5].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[1].p, mass[6].p), 1, 6, dist(mass[1].p, mass[6].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[1].p, mass[7].p), 1, 7, dist(mass[1].p, mass[7].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[2].p, mass[3].p), 2, 3, dist(mass[2].p, mass[3].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[2].p, mass[4].p), 2, 4, dist(mass[2].p, mass[4].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[2].p, mass[5].p), 2, 5, dist(mass[2].p, mass[5].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[2].p, mass[6].p), 2, 6, dist(mass[2].p, mass[6].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[2].p, mass[7].p), 2, 7, dist(mass[2].p, mass[7].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[3].p, mass[4].p), 3, 4, dist(mass[3].p, mass[4].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[3].p, mass[5].p), 3, 5, dist(mass[3].p, mass[5].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[3].p, mass[6].p), 3, 6, dist(mass[3].p, mass[6].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[3].p, mass[7].p), 3, 7, dist(mass[3].p, mass[7].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[4].p, mass[5].p), 4, 5, dist(mass[4].p, mass[5].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[4].p, mass[6].p), 4, 6, dist(mass[4].p, mass[6].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[4].p, mass[7].p), 4, 7, dist(mass[4].p, mass[7].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[5].p, mass[6].p), 5, 6, dist(mass[5].p, mass[6].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[5].p, mass[7].p), 5, 7, dist(mass[5].p, mass[7].p), b_val(mt), c_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[6].p, mass[7].p), 6, 7, dist(mass[6].p, mass[7].p), b_val(mt), c_val(mt)});

    Cube cube = {mass, spring};

    return cube;
}
