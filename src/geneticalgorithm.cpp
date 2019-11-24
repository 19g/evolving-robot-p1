#include <iostream>
#include <cmath>
#include <thread>
#include <random>
#include <vector>
#include "simulate.hpp"
#include "geneticalgorithm.hpp"

using namespace std;

int main() {
    // use threads for performance
    thread t[NUM_OF_TRIALS];
    for (int i = 0; i < NUM_OF_TRIALS; i++) {
        t[i] = thread(loop);
    }
    for (int i = 0; i < NUM_OF_TRIALS; i++) {
        t[i].join();
    }

    return 0;
}

void loop() {
    // begin timer
    clock_t begin = clock();

    // initialize population randomly
    vector<Cube> population(POP_SIZE);
    for (int i = 0; i < POP_SIZE; i++) {
        Cube temp = initialize_cube();
        population[i] = temp;
    }

    for (int i = 0; i < POP_SIZE; i++) {
        print_mass(population[i].mass[0]);
    }
    // evolutionary loop
    for (int iteration = 0; iteration < NUM_OF_ITERATIONS; iteration++) {



    }

    // end timer
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    double iters_per_sec = NUM_OF_ITERATIONS / elapsed_secs;
    cout << "iter/sec: " << iters_per_sec << "\n";

    // output to file for opengl
}

Cube initialize_cube() {
    // random variable generation
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> a_val(MIN_A, MAX_A);
    uniform_real_distribution<double> b_val(MIN_B, MAX_B);
    uniform_real_distribution<double> c_val(MIN_C, MAX_C);
    uniform_real_distribution<double> d_val(MIN_D, MAX_D);
    uniform_real_distribution<double> e_val(MIN_E, MAX_E);
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

    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[0].p, mass[1].p), 0, 1, a_val(mt), b_val(mt), c_val(mt), d_val(mt), e_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[0].p, mass[2].p), 0, 2, a_val(mt), b_val(mt), c_val(mt), d_val(mt), e_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[0].p, mass[3].p), 0, 3, a_val(mt), b_val(mt), c_val(mt), d_val(mt), e_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[0].p, mass[4].p), 0, 4, a_val(mt), b_val(mt), c_val(mt), d_val(mt), e_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[0].p, mass[5].p), 0, 5, a_val(mt), b_val(mt), c_val(mt), d_val(mt), e_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[0].p, mass[6].p), 0, 6, a_val(mt), b_val(mt), c_val(mt), d_val(mt), e_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[0].p, mass[7].p), 0, 7, a_val(mt), b_val(mt), c_val(mt), d_val(mt), e_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[1].p, mass[2].p), 1, 2, a_val(mt), b_val(mt), c_val(mt), d_val(mt), e_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[1].p, mass[3].p), 1, 3, a_val(mt), b_val(mt), c_val(mt), d_val(mt), e_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[1].p, mass[4].p), 1, 4, a_val(mt), b_val(mt), c_val(mt), d_val(mt), e_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[1].p, mass[5].p), 1, 5, a_val(mt), b_val(mt), c_val(mt), d_val(mt), e_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[1].p, mass[6].p), 1, 6, a_val(mt), b_val(mt), c_val(mt), d_val(mt), e_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[1].p, mass[7].p), 1, 7, a_val(mt), b_val(mt), c_val(mt), d_val(mt), e_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[2].p, mass[3].p), 2, 3, a_val(mt), b_val(mt), c_val(mt), d_val(mt), e_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[2].p, mass[4].p), 2, 4, a_val(mt), b_val(mt), c_val(mt), d_val(mt), e_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[2].p, mass[5].p), 2, 5, a_val(mt), b_val(mt), c_val(mt), d_val(mt), e_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[2].p, mass[6].p), 2, 6, a_val(mt), b_val(mt), c_val(mt), d_val(mt), e_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[2].p, mass[7].p), 2, 7, a_val(mt), b_val(mt), c_val(mt), d_val(mt), e_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[3].p, mass[4].p), 3, 4, a_val(mt), b_val(mt), c_val(mt), d_val(mt), e_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[3].p, mass[5].p), 3, 5, a_val(mt), b_val(mt), c_val(mt), d_val(mt), e_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[3].p, mass[6].p), 3, 6, a_val(mt), b_val(mt), c_val(mt), d_val(mt), e_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[3].p, mass[7].p), 3, 7, a_val(mt), b_val(mt), c_val(mt), d_val(mt), e_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[4].p, mass[5].p), 4, 5, a_val(mt), b_val(mt), c_val(mt), d_val(mt), e_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[4].p, mass[6].p), 4, 6, a_val(mt), b_val(mt), c_val(mt), d_val(mt), e_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[4].p, mass[7].p), 4, 7, a_val(mt), b_val(mt), c_val(mt), d_val(mt), e_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[5].p, mass[6].p), 5, 6, a_val(mt), b_val(mt), c_val(mt), d_val(mt), e_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[5].p, mass[7].p), 5, 7, a_val(mt), b_val(mt), c_val(mt), d_val(mt), e_val(mt)});
    spring.emplace_back(Spring {k_spring_val(mt), dist(mass[6].p, mass[7].p), 6, 7, a_val(mt), b_val(mt), c_val(mt), d_val(mt), e_val(mt)});

    Cube cube = {mass, spring};

    return cube;
}