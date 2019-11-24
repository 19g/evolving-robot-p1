#include <iostream>
#include <cmath>
#include <thread>
#include <random>
#include <vector>
#include "simulate.hpp"
#include "geneticalgorithm.hpp"

using namespace std;

int main() {
    srand(time(nullptr));

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

    // initialize parent population randomly
    vector<Cube> parents(POP_SIZE);
    for (int i = 0; i < POP_SIZE; i++) {
        Cube temp = initialize_cube();
        parents[i] = temp;
    }

    for (int i = 0; i < POP_SIZE; i++) {
        print_mass(parents[i].mass[0]);
    }
    // evolutionary loop
    for (int eval = 0; eval < NUM_OF_EVALS; eval++) {
        // get random order of individuals for crossover
        vector<int> order = randomize_array_of_springs();
        // initialize offspring
        vector<Cube> children(parents);

        // crossover
        for (int i = 0; i < POP_SIZE; i++) {
            crossover(parents[order[i]], parents[order[i+1]]);
        }

        // mutation
        for (int i = 0; i < POP_SIZE; i++) {
            mutation(children[i]);
        }

        // get fitness of population
        for (int i = 0; i < POP_SIZE; i++) {
            simulation_loop(parents[i].mass, parents[i].spring);
            simulation_loop(children[i].mass, children[i].spring);
        }

        // selection
        tournament_selection(parents, children);


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

vector<int> randomize_array_of_springs() {
    // create vector with order of parents to be crossed over
    vector<int> order(NUM_OF_SPRINGS);
    for (int i = 0; i < NUM_OF_SPRINGS; i++) {
        order[i] = i;
    }
    for (int i = order.size() - 1; i > 0; i--) {
        int j = (int)(rand() % (i+1));
        int temp = order[i];
        order[i] = order[j];
        order[j] = temp;
    }

    return order;
}

void crossover(Cube &A, Cube &B) {
    // choose random points for crossover
    int start = rand() % NUM_OF_SPRINGS;
    int end = rand() % NUM_OF_SPRINGS;
    while (end == start) {
        end = rand() % NUM_OF_SPRINGS;
    }

    // swap the spring values at that point
    for (int i = start; i < end; i++) {
        Spring temp = A.spring[i];
        A.spring[i] = B.spring[i];
        B.spring[i] = temp;
    }
}

void mutation(Cube &individual) {
    // random variable generation
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> mut_chance(0, 1);

    for (int i = 0; i < NUM_OF_SPRINGS; i++) {
        if (mut_chance(mt) < PROB_OF_MUT) {
            uniform_real_distribution<int> spring(0, NUM_OF_SPRINGS);
            uniform_real_distribution<double> swing(MIN_SWING, MAX_SWING);
            individual.spring[spring(mt)].a =  individual.spring[spring(mt)].a * swing(mt);
            individual.spring[spring(mt)].b =  individual.spring[spring(mt)].b * swing(mt);
            individual.spring[spring(mt)].c =  individual.spring[spring(mt)].c * swing(mt);
            individual.spring[spring(mt)].d =  individual.spring[spring(mt)].d * swing(mt);
            individual.spring[spring(mt)].e =  individual.spring[spring(mt)].e * swing(mt);
            individual.spring[spring(mt)].k =  individual.spring[spring(mt)].k * swing(mt);
        }
    }
}

void tournament_selection(vector<Cube> &parents, vector<Cube> &children) {
    // take elite children;
    int best_parent_index;
    int best_child_index;
}
















